#include <stdio.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include "incbin.h"
#include <time.h>
#include <math.h>
#include <alsa/asoundlib.h>
#include <signal.h>
#include <pthread.h>
#include <fftw3.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>

// ====== Définitions des constantes ======
#define NUM_COLS 31
#define LR_BAR_Y   192
#define LR_BAR_H   800
#define FFT_SIZE 16384
#define FFT_SIZE_b2 ((FFT_SIZE/2)+1)
#define HOP_SIZE 1024
#define AUDIO_BUFFER_SIZE HOP_SIZE

#define NUM_SCALE 9

const float db_mark[NUM_SCALE] = { +6, 0, -3, -8, -14, -20, -30, -40, -80 };
const int   y_mark[NUM_SCALE]  = { 192, 264, 308, 392, 496, 592, 752, 888, 992 };

const float sample_rate = 48000.0f;
const float freq_per_bin = 48000.0f / FFT_SIZE;

const float band_freqs[NUM_COLS] = {
    20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250,
    315, 400, 500, 630, 800, 1000, 1250, 1600, 2000,
    2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500, 16000, 20000
};

const float band_edges[NUM_COLS+1] = {
    16.00, 22.36, 28.07, 35.48, 44.72, 56.12, 71.05, 89.44,
    111.80, 141.42, 179.60, 223.61, 279.51, 354.81, 447.21, 561.23,
    707.11, 894.43, 1118.03, 1409.49, 1788.85, 2236.07, 2795.08,
    3548.13, 4472.14, 5612.49, 7071.07, 8944.27, 11180.34,
    14142.14, 17888.54, 22400.00
};

// ====== Données de calibration ======
// Référence FFT par bande (à mesurer)
static float ref_val[NUM_COLS] = {
    53.03, 37.27, 41.88, 34.97, 55.34, 39.77, 19.11, 35.20, 40.41, 38.98,
    38.21, 46.98, 45.60, 41.15, 39.22, 39.04, 36.95, 40.23, 40.43, 37.16,
    41.96, 40.30, 39.98, 39.72, 38.29, 39.82, 40.04, 39.42, 41.26, 40.13, 41.48
};

bool calibrate_spectro = 0;  // Flag calibration spectrogramme

float ref_vu_L = 1.0f;  // Calibration vumètre gauche
float ref_vu_R = 1.0f;  // Calibration vumètre droit
bool calibrate_vu = 0;  // Flag calibration vumètre

// ====== Variables globales pour audio et FFT ======
float hann_window[FFT_SIZE];
fftwf_plan fft_plan;
float FFT_in[FFT_SIZE];
fftwf_complex FFT_out[FFT_SIZE_b2];

volatile sig_atomic_t running = 1;

short audio_ringbuffer[FFT_SIZE * 2]; // Buffer circulaire stéréo interleaved
int audio_ringbuffer_pos = 0;

short fft_samples[FFT_SIZE * 2];      // Buffer pour FFT stéréo (L,R,...)
pthread_mutex_t fft_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t fft_cond = PTHREAD_COND_INITIALIZER;
volatile int new_block_ready = 0;

float col_values[NUM_COLS] = {0.0f};  // Valeurs niveaux par bande FFT
float spectro_peak_hold[NUM_COLS] = {0};

// Constantes de temps en secondes
float tau_peak_hold = 2.0f; //persistance des pics
float tau_rms = 0.1f;          // 100 ms RMS
float tau_peak_release = 0.3f; // 300 ms Peak relâchement

// Positions et dimensions des barres du spectrogramme
const int bar_x[NUM_COLS] = {
    120,160,200,240,280,320,360,400,440,480,520,560,600,640,680,720,
    760,800,840,880,920,960,1000,1040,1080,1120,1160,1200,1240,1280,1320
};
const int bar_y = 192;
const int bar_w = 22;
const int bar_h = 800;

// Positions et dimensions des barres L, S, R
const int lr_bar_x[3] = { 1482, 1527, 1558 };
const int lr_bar_w[3] = { 36, 22, 36 };
float lr_values[3] = {0.0f};       // Valeurs RMS L, S, R
float lr_peak_values[3] = {0.0f};  // Valeurs peak L, S, R
float balance_value = 0.0f;

const char *lr_labels[3] = {"L", "S", "R"};

// Dégradés pour affichage des barres
typedef struct { Uint8 r, g, b; } Color;

Color gradientStops[8] = {
    {0, 0, 0},      {0, 0, 128},    {0, 0, 255},    {0, 128, 0},
    {0, 255, 0},    {255, 255, 0},  {255, 128, 0},  {255, 0, 0}
};
static SDL_Texture *gradTex = NULL;

Color gradientLUT[800]; // LUT de dégradé verticale pour barres (hauteur = 800)

const int circle_cy = 172;
const int circle_radius = 11;

// ====== Prototypes ======
void handle_sigint(int sig);
SDL_Texture *loadPngFromMemory(SDL_Renderer *renderer);
Color getGradientColor(float v);
int circle_cx(int i);
void getCircleColor(float v, Uint8 *r, Uint8 *g, Uint8 *b);
void drawFilledCircle(SDL_Renderer *renderer, int x0, int y0, int radius, Uint8 r, Uint8 g, Uint8 b, Uint8 a);
void getGradLeftColor(float t, Uint8 *r, Uint8 *g, Uint8 *b);
void getGradRightColor(float t, Uint8 *r, Uint8 *g, Uint8 *b);
void draw_spectrogram_bars(SDL_Renderer *ren, const float *col_values);
void draw_lr_bars(SDL_Renderer *ren, const float *rms_values, const float *peak_values);
void draw_balance_bar(SDL_Renderer *ren, float value);
snd_pcm_t *init_alsa_capture(const char *device);
void update_audio_levels(snd_pcm_t *handle, float *lr_values, float *balance_value);
void* fft_thread_func(void *arg);
void save_calibration_file();
int load_calibration_file();
void init_hann_window();
void init_gradientLUT();
int calc_filled_height(float val);

// ====== Fonctions ======

// -- Handler signal Ctrl-C --
void handle_sigint(int sig) {
    running = 0;
}

// -- Chargement PNG en mémoire avec INCBIN --
INCBIN(background_png, "backgroundVibeScope.png");

SDL_Texture *loadPngFromMemory(SDL_Renderer *renderer) {
    const unsigned char *data = _binary_background_png_start;
    size_t len = _binary_background_png_end - _binary_background_png_start;
    SDL_RWops *rw = SDL_RWFromConstMem(data, len);
    if (!rw) return NULL;
    SDL_Surface *surface = IMG_Load_RW(rw, 1);
    if (!surface) return NULL;
    SDL_Texture *tex = SDL_CreateTextureFromSurface(renderer, surface);
    SDL_FreeSurface(surface);
    return tex;
}

// -- Calcul du dégradé vertical sur les barres spectro --
Color getGradientColor(float v) {
    if (v <= 0) return gradientStops[0];
    if (v >= 1) return gradientStops[7];
    float pos = v * 7.0f;
    int idx = (int)pos;
    float t = pos - idx;
    Color a = gradientStops[idx];
    Color b = gradientStops[idx + 1];
    Color c;
    c.r = a.r + (b.r - a.r) * t;
    c.g = a.g + (b.g - a.g) * t;
    c.b = a.b + (b.b - a.b) * t;
    return c;
}

// -- Coordonnées en X des cercles indicateurs --
int circle_cx(int i) { return bar_x[i] + bar_w / 2; }

// -- Couleur cercle indicateur selon niveau (vert-orange-rouge) --
void getCircleColor(float v, Uint8 *r, Uint8 *g, Uint8 *b) {
    float level_db = -60.0f + v * 66.0f;
    if (level_db < 0.0f) {
        *r = 0; *g = 255; *b = 0;
    } else if (level_db < 3.0f) {
        *r = 255; *g = 140; *b = 0;
    } else {
        *r = 255; *g = 0; *b = 0;
    }
}

// -- Dessine un cercle plein (pour indicateurs) --
void drawFilledCircle(SDL_Renderer *renderer, int x0, int y0, int radius, Uint8 r, Uint8 g, Uint8 b, Uint8 a) {
    SDL_SetRenderDrawColor(renderer, r, g, b, a);
    for (int y = -radius; y <= radius; y++) {
        for (int x = -radius; x <= radius; x++) {
            if (x*x + y*y <= radius*radius)
                SDL_RenderDrawPoint(renderer, x0 + x, y0 + y);
        }
    }
}

// -- Dégradé pour barre balance gauche --
void getGradLeftColor(float t, Uint8 *r, Uint8 *g, Uint8 *b) {
    if (t < 0.5f) {
        float k = t / 0.5f;
        *r = 255;
        *g = (Uint8)(255 * k);
        *b = 0;
    } else {
        float k = (t - 0.5f) / 0.5f;
        *r = (Uint8)(255 * (1 - k));
        *g = 255;
        *b = 0;
    }
}

// -- Dégradé pour barre balance droite --
void getGradRightColor(float t, Uint8 *r, Uint8 *g, Uint8 *b) {
    if (t < 0.5f) {
        float k = t / 0.5f;
        *r = (Uint8)(255 * k);
        *g = 255;
        *b = 0;
    } else {
        float k = (t - 0.5f) / 0.5f;
        *r = 255;
        *g = (Uint8)(255 * (1 - k));
        *b = 0;
    }
}

// -- Conversion niveau dB en coordonnée verticale --
int db_to_y(float db) {
    if (db >= db_mark[0]) return y_mark[0];
    if (db <= db_mark[NUM_SCALE-1]) return y_mark[NUM_SCALE-1];
    for (int i = 0; i < NUM_SCALE-1; ++i) {
        if (db <= db_mark[i] && db > db_mark[i+1]) {
            float t = (db - db_mark[i+1]) / (db_mark[i] - db_mark[i+1]);
            return (int)(y_mark[i+1] + t * (y_mark[i] - y_mark[i+1]));
        }
    }
    return y_mark[NUM_SCALE-1];
}

// -- Dessine les barres du spectrogramme --
void draw_spectrogram_bars(SDL_Renderer *ren, const float *col_values) {
    // Static pour garder les valeurs de peak hold du spectro
    static float draw_spectro_peak_hold[NUM_COLS] = {0};
    static int first_call = 1;

    // On met à jour draw_spectro_peak_hold avec col_values ici
    // (Attention, dans fft_thread_func on a aussi un peak hold, donc ici on synchronise avec celui-là)
    // Pour la simplicité, on met à jour ici aussi avec un lissage léger (ou copier le tableau depuis fft_thread_func si possible)

    if (first_call) {
        // Initialisation à la première passe
        for (int i = 0; i < NUM_COLS; i++) {
            draw_spectro_peak_hold[i] = col_values[i];
        }
        first_call = 0;
    } else {
        float cur_time = SDL_GetTicks() * 0.001f;
        static float last_update_time = 0.0f;
        float dt = cur_time - last_update_time;
        if (dt <= 0.0f) dt = 0.001f;
        float alpha_decay = 1.0f - expf(-dt / tau_peak_hold);

        for (int i = 0; i < NUM_COLS; i++) {
            if (col_values[i] > draw_spectro_peak_hold[i]) {
                draw_spectro_peak_hold[i] = col_values[i];
            } else {
                draw_spectro_peak_hold[i] = (1 - alpha_decay) * draw_spectro_peak_hold[i] + alpha_decay * col_values[i];
                if (draw_spectro_peak_hold[i] < col_values[i])
                    draw_spectro_peak_hold[i] = col_values[i];
            }
        }
        last_update_time = cur_time;
    }

    for (int i = 0; i < NUM_COLS; ++i) {
        SDL_Rect rc = {bar_x[i], bar_y, bar_w, bar_h};
        SDL_SetRenderDrawColor(ren, 255, 255, 255, 255);
        SDL_RenderDrawRect(ren, &rc);

        int y_head = db_to_y(col_values[i]);
        int y_bottom = bar_y + bar_h - 1;

        // Remplissage du bas vers le haut jusqu'à y_head
        int filled_h = (bar_y + bar_h) - y_head;     // hauteur à remplir
        if (filled_h > 0) {
		SDL_Rect src = {0, 0, 1, filled_h};                   // on prend le haut (0…filled_h-1)
		SDL_Rect dst = {bar_x[i] + 1, y_head, bar_w - 2, filled_h};
		SDL_RenderCopyEx(ren,          // même coût qu’un RenderCopy
						 gradTex,
						 &src,
						 &dst,
						 0.0,          // pas de rotation
						 NULL,
						 SDL_FLIP_VERTICAL);   // ← inverse l’axe Y
		}

        // Ligne du peak hold (petit trait blanc semi-transparent)
        int y_peak = db_to_y(draw_spectro_peak_hold[i]);
        if (y_peak < bar_y) y_peak = bar_y;
        if (y_peak > y_bottom) y_peak = y_bottom;

        SDL_SetRenderDrawColor(ren, 255, 255, 255, 180);
        SDL_Rect peak_rect = {bar_x[i] + 2, y_peak, bar_w - 4, 2};
        SDL_RenderFillRect(ren, &peak_rect);

        // Cercles indicateurs de niveau
        Uint8 r, g, b;
        float v = (col_values[i] + 60.0f) / 66.0f;
        getCircleColor(v, &r, &g, &b);
        int cx = circle_cx(i);
        int cy = circle_cy;
        drawFilledCircle(ren, cx, cy, circle_radius, r, g, b, 255);
    }
}


// -- Calcul hauteur remplie pour barres RMS ou Peak (en pixels) --
int calc_filled_height(float val) {
    float val_db = 20.0f * log10f(fmaxf(val, 1e-5f));
    int y_pos = db_to_y(val_db);
    int filled_h = 0;
    if (y_pos < LR_BAR_Y) {
        filled_h = LR_BAR_H;
    } else if (y_pos > LR_BAR_Y + LR_BAR_H) {
        filled_h = 0;
    } else {
        filled_h = (LR_BAR_Y + LR_BAR_H) - y_pos;
    }
    if (filled_h < 0) filled_h = 0;
    if (filled_h > LR_BAR_H) filled_h = LR_BAR_H;
    return filled_h;
}

// -- Dessine les barres RMS (opaque) et Peak (semi-transparent) pour L, S, R --
/* === draw_lr_bars : version identique au spectro === */
void draw_lr_bars(SDL_Renderer *ren,
                  const float *rms_values,
                  const float *peak_values)
{
    for (int i = 0; i < 3; ++i) {
        /* Cadre blanc */
        SDL_Rect rc = { lr_bar_x[i], LR_BAR_Y, lr_bar_w[i], LR_BAR_H };
        SDL_SetRenderDrawColor(ren, 255, 255, 255, 255);
        SDL_RenderDrawRect(ren, &rc);

        /* ------------------ RMS plein ------------------ */
        int filled = calc_filled_height(rms_values[i]);
        if (filled) {
            SDL_Rect src = { 0, 0, 1, filled };           /* bas de la texture */
            SDL_Rect dst = { lr_bar_x[i] + 1,
                             LR_BAR_Y + LR_BAR_H - filled,
                             lr_bar_w[i] - 2,
                             filled };
            SDL_SetTextureAlphaMod(gradTex, 255);          /* opaque */
            SDL_RenderCopyEx(ren, gradTex, &src, &dst,
                             0.0, NULL, SDL_FLIP_VERTICAL); /* noir en bas */
        }

        /* ------------------ Peak semi-transparent ------------------ */
        int filled_peak = calc_filled_height(peak_values[i]);
        if (filled_peak) {
            SDL_Rect src = { 0, 0, 1, filled_peak };
            SDL_Rect dst = { lr_bar_x[i] + 1,
                             LR_BAR_Y + LR_BAR_H - filled_peak,
                             lr_bar_w[i] - 2,
                             filled_peak };
            SDL_SetTextureAlphaMod(gradTex, 128);          /* 50 % */
            SDL_RenderCopyEx(ren, gradTex, &src, &dst,
                             0.0, NULL, SDL_FLIP_VERTICAL);
        }
    }
    /* Remise de l’opacité pour la suite (spectrogramme) */
    SDL_SetTextureAlphaMod(gradTex, 255);
}


// -- Dessine la barre de balance L/R avec dégradé coloré --
void draw_balance_bar(SDL_Renderer *ren, float value) {
    if (value < -1.0f) value = -1.0f;
    if (value >  1.0f) value =  1.0f;

    int left_x  = 1482, left_w  = 56;
    int right_x = 1538, right_w = 56;
    int bar_y   = 164,  bar_h   = 22;

    if (value < 0.0f) {
        int fill = (int)(left_w * (-value));
        for (int dx = 0; dx < fill; ++dx) {
            float t = (left_w - 1 - dx) / (float)(left_w - 1);
            Uint8 r, g, b;
            getGradLeftColor(t, &r, &g, &b);
            SDL_SetRenderDrawColor(ren, r, g, b, 255);
            SDL_RenderDrawLine(ren, left_x + (left_w - 1 - dx), bar_y, left_x + (left_w - 1 - dx), bar_y + bar_h - 1);
        }
    }
    else if (value > 0.0f) {
        int fill = (int)(right_w * value);
        for (int dx = 0; dx < fill; ++dx) {
            float t = dx / (float)(right_w - 1);
            Uint8 r, g, b;
            getGradRightColor(t, &r, &g, &b);
            SDL_SetRenderDrawColor(ren, r, g, b, 255);
            SDL_RenderDrawLine(ren, right_x + dx, bar_y, right_x + dx, bar_y + bar_h - 1);
        }
    }
}

// -- Initialisation de la fenêtre de Hann pour FFT --
void init_hann_window() {
    for (int i = 0; i < FFT_SIZE; i++)
        hann_window[i] = 0.5f * (1.0f - cosf(2.0f * M_PI * i / (FFT_SIZE - 1)));
}

// -- Initialisation du LUT dégradé vertical --
// --- Création du LUT EN MÉMOIRE + texture GPU en même temps ---
void init_gradient_texture(SDL_Renderer *ren)
{
    // 1. Surface RGBA 1 x bar_h
    SDL_Surface *surf = SDL_CreateRGBSurfaceWithFormat(
        0, 1, bar_h, 32, SDL_PIXELFORMAT_RGBA8888);
    if (!surf) {
        fprintf(stderr, "SDL_CreateRGBSurface error: %s\n", SDL_GetError());
        exit(1);
    }

    Uint32 *pix = (Uint32 *)surf->pixels;
    for (int dy = 0; dy < bar_h; ++dy) {
        float val = (float)dy / (bar_h - 1);
        Color c = getGradientColor(val);          // ta table 8 couleurs
        pix[dy] = SDL_MapRGBA(surf->format, c.r, c.g, c.b, 255);
        gradientLUT[dy] = c;                      // ← garde aussi le LUT si d’autres fonctions l’utilisent
    }

    // 2. Texture GPU
    gradTex = SDL_CreateTextureFromSurface(ren, surf);
    SDL_FreeSurface(surf);
    if (!gradTex) {
        fprintf(stderr, "SDL_CreateTextureFromSurface error: %s\n", SDL_GetError());
        exit(1);
    }
    SDL_SetTextureBlendMode(gradTex, SDL_BLENDMODE_BLEND);  // autorise alpha pour le peak
}


// -- Initialisation ALSA capture (input audio) --
snd_pcm_t *init_alsa_capture(const char *device) {
    snd_pcm_t *pcm_handle = NULL;
    snd_pcm_hw_params_t *hw_params = NULL;
    int err;

    if ((err = snd_pcm_open(&pcm_handle, device, SND_PCM_STREAM_CAPTURE, 0)) < 0) {
        fprintf(stderr, "ALSA open error (%s): %s\n", device, snd_strerror(err));
        return NULL;
    }

    snd_pcm_hw_params_alloca(&hw_params);
    if ((err = snd_pcm_hw_params_any(pcm_handle, hw_params)) < 0) {
        fprintf(stderr, "ALSA: Impossible d'obtenir les paramètres hardware: %s\n", snd_strerror(err));
        goto fail;
    }
    if ((err = snd_pcm_hw_params_set_format(pcm_handle, hw_params, SND_PCM_FORMAT_S16_LE)) < 0) {
        fprintf(stderr, "ALSA: Impossible de définir le format S16_LE: %s\n", snd_strerror(err));
        goto fail;
    }
    if ((err = snd_pcm_hw_params_set_channels(pcm_handle, hw_params, 2)) < 0) {
        fprintf(stderr, "ALSA: Impossible de définir 2 canaux: %s\n", snd_strerror(err));
        goto fail;
    }
    unsigned int rate = (unsigned int)sample_rate;
    if ((err = snd_pcm_hw_params_set_rate_near(pcm_handle, hw_params, &rate, 0)) < 0) {
        fprintf(stderr, "ALSA: Impossible de définir la fréquence d'échantillonnage: %s\n", snd_strerror(err));
        goto fail;
    }
    if (rate != 48000)
        fprintf(stderr, "ALSA: Fréquence obtenue: %u Hz (au lieu de 48000 Hz)\n", rate);

    if ((err = snd_pcm_hw_params_set_access(pcm_handle, hw_params, SND_PCM_ACCESS_RW_INTERLEAVED)) < 0) {
        fprintf(stderr, "ALSA: Impossible de définir l'accès interleaved: %s\n", snd_strerror(err));
        goto fail;
    }

    snd_pcm_uframes_t buffer_size = AUDIO_BUFFER_SIZE;
    snd_pcm_hw_params_set_period_size_near(pcm_handle, hw_params, &buffer_size, 0);

    if ((err = snd_pcm_hw_params(pcm_handle, hw_params)) < 0) {
        fprintf(stderr, "ALSA: Impossible d'appliquer les paramètres hardware: %s\n", snd_strerror(err));
        goto fail;
    }

    if ((err = snd_pcm_prepare(pcm_handle)) < 0) {
        fprintf(stderr, "ALSA: snd_pcm_prepare() a échoué: %s\n", snd_strerror(err));
        goto fail;
    }

    fprintf(stderr, "ALSA: Capture initialisée sur %s\n", device);
    return pcm_handle;

fail:
    if (pcm_handle)
        snd_pcm_close(pcm_handle);
    return NULL;
}

// -- Mise à jour des niveaux audio RMS, Peak et balance --
void update_audio_levels(snd_pcm_t *handle, float *lr_values, float *balance_value) {
    static float smoothed_rms_l = 0.0f;
    static float smoothed_rms_r = 0.0f;
    static float smoothed_peak_l = 0.0f;
    static float smoothed_peak_r = 0.0f;
    static float last_time = 0.0f;
    static int waiting_for_stability = 0;
    static float stable_rms_l = 0.0f;
    static float stable_rms_r = 0.0f;
    static float stable_since = 0.0f;

    short buffer[HOP_SIZE * 2];
    int nframes;

    while (1) {
        nframes = snd_pcm_readi(handle, buffer, HOP_SIZE);
        if (nframes == -EAGAIN) { SDL_Delay(1); continue; }
        if (nframes == -EPIPE) { snd_pcm_prepare(handle); continue; }
        if (nframes < 0) {
            fprintf(stderr, "Erreur de lecture ALSA: %s\n", snd_strerror(nframes));
            lr_values[0] = lr_values[1] = lr_values[2] = 0.0f;
            lr_peak_values[0] = lr_peak_values[1] = lr_peak_values[2] = 0.0f;
            *balance_value = 0.0f;
            return;
        }
        break;
    }

    // Calcul RMS et peak bruts
    double sum_l = 0, sum_r = 0;
    float peak_l = 0.0f;
    float peak_r = 0.0f;

    for (int i = 0; i < nframes; ++i) {
        float sample_l = buffer[2*i] / 32768.0f;
        float sample_r = buffer[2*i+1] / 32768.0f;

        sum_l += sample_l * sample_l;
        sum_r += sample_r * sample_r;

        float abs_l = fabsf(sample_l);
        float abs_r = fabsf(sample_r);
        if (abs_l > peak_l) peak_l = abs_l;
        if (abs_r > peak_r) peak_r = abs_r;
    }
    float rms_l = sqrt(sum_l / nframes);
    float rms_r = sqrt(sum_r / nframes);

    // Gestion du temps
    float cur_time = SDL_GetTicks() * 0.001f; // secondes
    float dt = cur_time - last_time;
    if (dt <= 0.0f) dt = 0.001f; // éviter div0 ou négatif

    // Calibration vumètre (idem que dans ta version)

    if (calibrate_vu) {
        float threshold = 0.01f;
        if (rms_l > threshold && rms_r > threshold) {
            if (!waiting_for_stability) {
                stable_since = cur_time;
                stable_rms_l = rms_l;
                stable_rms_r = rms_r;
                waiting_for_stability = 1;
                printf("Calibration VU: attente stabilité...\n");
            } else {
                if (cur_time - stable_since >= 1.0f) {
                    ref_vu_L = stable_rms_l;
                    ref_vu_R = stable_rms_r;
                    printf("[CALIB VU] ref_vu_L = %.6f, ref_vu_R = %.6f\n", ref_vu_L, ref_vu_R);
                    calibrate_vu = 0;
                    save_calibration_file();
                    waiting_for_stability = 0;
                }
            }
        } else {
            waiting_for_stability = 0;
        }
    }

    // Calcul alpha
    float alpha_rms = 1.0f - expf(-dt / tau_rms);

    // Lissage RMS
    smoothed_rms_l = alpha_rms * rms_l + (1.0f - alpha_rms) * smoothed_rms_l;
    smoothed_rms_r = alpha_rms * rms_r + (1.0f - alpha_rms) * smoothed_rms_r;

    // Lissage Peak, attaque instantanée, relâchement lissé
    if (peak_l > smoothed_peak_l) {
        smoothed_peak_l = peak_l;
    } else {
        float alpha_peak = 1.0f - expf(-dt / tau_peak_release);
        smoothed_peak_l = alpha_peak * peak_l + (1.0f - alpha_peak) * smoothed_peak_l;
    }

    if (peak_r > smoothed_peak_r) {
        smoothed_peak_r = peak_r;
    } else {
        float alpha_peak = 1.0f - expf(-dt / tau_peak_release);
        smoothed_peak_r = alpha_peak * peak_r + (1.0f - alpha_peak) * smoothed_peak_r;
    }

    last_time = cur_time;

    // Normalisation RMS corrigée avec calibration
    float norm_l = (smoothed_rms_l / ref_vu_L) * 0.70710678f;
    float norm_r = (smoothed_rms_r / ref_vu_R) * 0.70710678f;
    float norm_s = (norm_l + norm_r) * 0.5f;

    if (norm_l > 1.0f) norm_l = 1.0f;
    if (norm_r > 1.0f) norm_r = 1.0f;
    if (norm_s > 1.0f) norm_s = 1.0f;

    lr_values[0] = norm_l;
    lr_values[2] = norm_r;
    lr_values[1] = norm_s;

    // Clamp peaks
    if (smoothed_peak_l > 1.0f) smoothed_peak_l = 1.0f;
    if (smoothed_peak_r > 1.0f) smoothed_peak_r = 1.0f;
    float peak_s = (smoothed_peak_l + smoothed_peak_r) * 0.5f;
    if (peak_s > 1.0f) peak_s = 1.0f;

    lr_peak_values[0] = smoothed_peak_l;
    lr_peak_values[2] = smoothed_peak_r;
    lr_peak_values[1] = peak_s;

    // Calcul balance
    float total = rms_l + rms_r;
    if (total > 1e-6)
        *balance_value = (rms_r - rms_l) / total;
    else
        *balance_value = 0.0f;

    // Écriture dans buffer circulaire et signal FFT (idem ta version)

    pthread_mutex_lock(&fft_mutex);

    int frames_to_write = HOP_SIZE;
    int pos = audio_ringbuffer_pos;

    for (int i = 0; i < frames_to_write; ++i) {
        int write_pos = (pos + i) % FFT_SIZE;
        audio_ringbuffer[write_pos*2 + 0] = buffer[i*2 + 0];
        audio_ringbuffer[write_pos*2 + 1] = buffer[i*2 + 1];
    }

    audio_ringbuffer_pos = (audio_ringbuffer_pos + frames_to_write) % FFT_SIZE;

    for (int i = 0; i < FFT_SIZE; ++i) {
        int src_pos = (audio_ringbuffer_pos + i) % FFT_SIZE;
        fft_samples[i*2 + 0] = audio_ringbuffer[src_pos*2 + 0];
        fft_samples[i*2 + 1] = audio_ringbuffer[src_pos*2 + 1];
    }

    new_block_ready = 1;
    pthread_cond_signal(&fft_cond);
    pthread_mutex_unlock(&fft_mutex);
}


// -- Thread FFT : calcul de la FFT et calibration spectre --
void* fft_thread_func(void *arg) {
    short block_stereo[FFT_SIZE * 2];

    static int last_band = -1;
    static float last_time = 0.0f;
    static int already_calibrated = -1;

    // Tableau pour stocker le peak hold (persistance des pics)
    static float spectro_peak_hold[NUM_COLS] = {0};

    // Temps de décroissance du peak hold en secondes (ex : 1 seconde)
    
    static float last_peak_update_time = 0.0f;

    while (running) {
        pthread_mutex_lock(&fft_mutex);
        while (!new_block_ready && running)
            pthread_cond_wait(&fft_cond, &fft_mutex);
        if (!running) {
            pthread_mutex_unlock(&fft_mutex);
            break;
        }
        memcpy(block_stereo, fft_samples, sizeof(short) * 2 * FFT_SIZE);
        new_block_ready = 0;
        pthread_mutex_unlock(&fft_mutex);

        // Conversion stéréo vers mono et application fenêtre Hann
        for (int i = 0; i < FFT_SIZE; i++) {
            float L = block_stereo[2*i] / 32768.0f;
            float R = block_stereo[2*i+1] / 32768.0f;
            FFT_in[i] = 0.5f * (L + R) * hann_window[i];
        }

        fftwf_execute_dft_r2c(fft_plan, FFT_in, FFT_out);

        static float mag[FFT_SIZE_b2];
        for (int i = 0; i < FFT_SIZE_b2; i++)mag[i] = hypotf(FFT_out[i][0], FFT_out[i][1]);

        // Calcul niveaux moyens par bande de fréquence
        for (int c = 0; c < NUM_COLS; c++) {
            float f_low  = band_edges[c];
            float f_high = band_edges[c+1];

            int bin_start = (int)ceilf(f_low / freq_per_bin);
            int bin_end   = (int)floorf(f_high / freq_per_bin);

            if (bin_start < 0) bin_start = 0;
            if (bin_end > FFT_SIZE_b2) bin_end = FFT_SIZE_b2;

            float sum = 0.0f;
			int n_bins = bin_end - bin_start + 1;
			for (int i = bin_start; i <= bin_end; i++) sum += mag[i];
            float avg = (n_bins > 0) ? (sum / n_bins) : 1e-20f;

            float avg_norm = avg / ref_val[c];
            float db = 20.0f * log10f(avg_norm + 1e-20f);
            col_values[c] = db;
        }

        // Mise à jour peak hold (persistance pics)
        float cur_time = SDL_GetTicks() * 0.001f; // secondes
        float dt = cur_time - last_peak_update_time;
        if (dt <= 0.0f) dt = 0.001f;
        float alpha_decay = 1.0f - expf(-dt / tau_peak_hold);

        for (int c = 0; c < NUM_COLS; c++) {
            if (col_values[c] > spectro_peak_hold[c]) {
                spectro_peak_hold[c] = col_values[c]; // nouveau pic max
            } else {
                spectro_peak_hold[c] = (1 - alpha_decay) * spectro_peak_hold[c] + alpha_decay * col_values[c];
                // éviter que le pic descende en dessous du niveau actuel
                if (spectro_peak_hold[c] < col_values[c])
                    spectro_peak_hold[c] = col_values[c];
            }
        }
        last_peak_update_time = cur_time;

        // Calibration spectrogramme simplifiée (1s stable sur la même bande)
        if (calibrate_spectro) {
            // Réutiliser cur_time défini ci-dessus

            // Trouver la bande la plus forte
            int best_band = -1;
            float max_db = -1000.0f;
            for (int c = 0; c < NUM_COLS; ++c) {
                if (col_values[c] > max_db) {
                    max_db = col_values[c];
                    best_band = c;
                }
            }

            float threshold_db = -20.0f;

            if (best_band != last_band || max_db < threshold_db) {
                last_band = best_band;
                last_time = cur_time;
                already_calibrated = -1;
            }

            if (best_band == last_band && max_db > threshold_db &&
                (cur_time - last_time > 1.0) && already_calibrated != best_band) {

                float f_low  = band_edges[best_band];
                float f_high = band_edges[best_band+1];
                int bin_start = (int)ceilf(f_low / freq_per_bin);
                int bin_end   = (int)floorf(f_high / freq_per_bin);
                float sum = 0.0f;
                int n_bins = bin_end - bin_start + 1;
                for (int i = bin_start; i <= bin_end; i++) sum += mag[i];
                float avg = (n_bins > 0) ? (sum / n_bins) : 1e-20f;

                printf("[CALIB] Bande %d (%.1f Hz) : ref_val[%d] = %.2f (avant: %.2f)\n",
                       best_band, band_freqs[best_band], best_band, avg, ref_val[best_band]);

                ref_val[best_band] = avg;
                if (best_band + 1 < NUM_COLS) {
                    ref_val[best_band + 1] = avg;
                }
                save_calibration_file();
                already_calibrated = best_band;
            }

            if (max_db < -20.0f) {
                already_calibrated = -1;
            }
        }

        usleep(1000);
    }
    return NULL;
}


// -- Sauvegarde des calibrations dans cal.msr --
void save_calibration_file() {
    FILE *f = fopen("cal.msr", "w");
    if (!f) {
        fprintf(stderr, "Erreur ouverture fichier cal.msr\n");
        return;
    }

    for (int i = 0; i < NUM_COLS; i++) {
        fprintf(f, "%.2f", ref_val[i]);
        if (i < NUM_COLS - 1)
            fprintf(f, ";");
    }
    fprintf(f, "\n");

    fprintf(f, "%.6f;%.6f\n", ref_vu_L, ref_vu_R);

    fclose(f);
}

fftwf_plan create_fftwf_plan_with_wisdom(int size, float* FFT_in, fftwf_complex* FFT_out) {
    fftwf_plan plan = NULL;
    const char* wisdom_filename = "wisdom.fftw";
    int wisdom_loaded = 0;
    // Tente de charger le wisdom
    FILE *wisdom_file = fopen(wisdom_filename, "r");
    if (wisdom_file) {
        fclose(wisdom_file);
        wisdom_loaded = fftwf_import_wisdom_from_filename(wisdom_filename);
        if (wisdom_loaded) {
            printf("[FFTW] Wisdom chargé depuis '%s'.\n", wisdom_filename);
        } else {
            printf("[FFTW] Fichier wisdom présent, mais échec du chargement. Il sera recalculé.\n");
        }
    } else {
        printf("[FFTW] Pas de wisdom détecté, création initiale...\n");
    }
    // Création du plan (toujours avec MEASURE pour garantir l'optimum)
    plan = fftwf_plan_dft_r2c_1d(size, FFT_in, FFT_out, FFTW_MEASURE);
    if (!plan) {
        fprintf(stderr, "[FFTW] Erreur : impossible de créer le plan FFTW !\n");
        exit(1);
    }
    // Si wisdom pas chargé, le sauver pour la prochaine fois
    if (!wisdom_loaded) {
        if (fftwf_export_wisdom_to_filename(wisdom_filename)) {
            printf("[FFTW] Wisdom sauvegardé dans '%s'.\n", wisdom_filename);
        } else {
            fprintf(stderr, "[FFTW] Erreur à la sauvegarde du wisdom !\n");
        }
    }
    return plan;
}


// -- Chargement des calibrations depuis cal.msr --
int load_calibration_file() {
    FILE *f = fopen("cal.msr", "r");
    if (!f) {
        fprintf(stderr, "Impossible d'ouvrir cal.msr\n");
        return -1;
    }

    char line[1024];

    if (!fgets(line, sizeof(line), f)) {
        fprintf(stderr, "Erreur lecture première ligne de cal.msr\n");
        fclose(f);
        return -1;
    }

    char *token = strtok(line, ";");
    int i = 0;
    while (token && i < NUM_COLS) {
        ref_val[i] = strtof(token, NULL);
        token = strtok(NULL, ";");
        i++;
    }

    if (i != NUM_COLS) {
        fprintf(stderr, "Nombre de valeurs dans cal.msr incorrect\n");
        fclose(f);
        return -1;
    }

    if (!fgets(line, sizeof(line), f)) {
        fprintf(stderr, "Erreur lecture deuxième ligne de cal.msr\n");
        fclose(f);
        return -1;
    }

    char *token_L = strtok(line, ";");
    char *token_R = strtok(NULL, ";");

    if (!token_L || !token_R) {
        fprintf(stderr, "Format calibration vumètre incorrect\n");
        fclose(f);
        return -1;
    }

    ref_vu_L = strtof(token_L, NULL);
    ref_vu_R = strtof(token_R, NULL);

    fclose(f);
    return 0;
}

void print_help(const char *progname) {
    printf("Usage: %s [options]\n", progname);
    printf("Options:\n");
    printf("  -D <device>           : ALSA audio capture device (default: 'default')\n");
    printf("  -CS                   : Active calibration spectrogram\nUse:\nfreqs=(20 25 31.5 40 50 63 80 100 125 160 200 250 315 400 500 630 800 1000 1250 1600 2000 2500 3150 4000 5000 6300 8000 10000 12500 16000 20000)\nfor f in \"${freqs[@]}\"; do\nsox -n -c 2 -t alsa plughw:CARD=loopTest1,DEV=1 synth 5 sine $f vol 0dB\ndone");
    printf("  -CV                   : Enable VU meter calibration\nUse:\nsox -n -t alsa -b 16 -e signed-integer -c 2 plughw:CARD=loopTest1,DEV=1 synth 60 sine 1000 vol 0dB\n");
    printf("  -PT <tau_peak_hold>   : Spectrogram peak persistence time (seconds, default: 2.0)\n");
    printf("  -TR <tau_rms>         : RMS smoothing time constant (seconds, default: 0.1)\n");
    printf("  -TP <tau_peak_release>: Peak release time constant (seconds, default: 0.3)\n");
    printf("  -h, --help            : Show this help\n");
}

// ====== Fonction principale ======
int main(int argc, char *argv[]) {

    const char *alsa_device = "default";

    // Lecture arguments ligne de commande
	for (int i = 1; i < argc; ++i) {
		if (strcmp(argv[i], "-D") == 0 && i + 1 < argc) {
			alsa_device = argv[i+1];
			i++;
		}
		else if (strcmp(argv[i], "-CS") == 0) {
			calibrate_spectro = 1;
		}
		else if (strcmp(argv[i], "-CV") == 0) {
			calibrate_vu = 1;
		}
		else if (strcmp(argv[i], "-PT") == 0 && i + 1 < argc) {
			tau_peak_hold = atof(argv[i+1]);
			if (tau_peak_hold < 0.1f) tau_peak_hold = 0.1f; // sécurité
			i++;
		}
	}

	
    fft_plan = create_fftwf_plan_with_wisdom(FFT_SIZE, FFT_in, FFT_out);
    init_hann_window();
    

    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        fprintf(stderr, "SDL_Init error: %s\n", SDL_GetError());
        return 1;
    }

    if (!(IMG_Init(IMG_INIT_PNG) & IMG_INIT_PNG)) {
        fprintf(stderr, "IMG_Init error: %s\n", IMG_GetError());
        SDL_Quit();
        return 1;
    }

    SDL_Window *win = SDL_CreateWindow("Audio Analyser",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
        1920, 1080, SDL_WINDOW_SHOWN);
    SDL_SetWindowFullscreen(win, SDL_WINDOW_FULLSCREEN_DESKTOP);
    SDL_Renderer *ren = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED);
	init_gradient_texture(ren);

    SDL_Texture *bg = loadPngFromMemory(ren);
    if (!bg) {
        fprintf(stderr, "Load memory PNG error: %s\n", SDL_GetError());
        SDL_DestroyRenderer(ren);
        SDL_DestroyWindow(win);
        IMG_Quit();
        SDL_Quit();
        return 1;
    }

    snd_pcm_t *pcm_handle = init_alsa_capture(alsa_device);
    if (!pcm_handle) {
        SDL_DestroyRenderer(ren);
        SDL_DestroyWindow(win);
        IMG_Quit();
        SDL_Quit();
        return 1;
    }

    if (load_calibration_file() == 0) {
        printf("Calibration file loaded successfully.\n");
    }

    pthread_t fft_thread;
    pthread_create(&fft_thread, NULL, fft_thread_func, NULL);

    signal(SIGINT, handle_sigint);

    while (running) {
        update_audio_levels(pcm_handle, lr_values, &balance_value);

        SDL_RenderClear(ren);
        SDL_RenderCopy(ren, bg, NULL, NULL);

        draw_spectrogram_bars(ren, col_values);
        draw_lr_bars(ren, lr_values, lr_peak_values);
        draw_balance_bar(ren, balance_value);

        SDL_RenderPresent(ren);
        SDL_Delay(1);
    }

    running = 0;
    pthread_cond_signal(&fft_cond);
    pthread_join(fft_thread, NULL);
    fftwf_destroy_plan(fft_plan);

    snd_pcm_close(pcm_handle);
    SDL_DestroyTexture(bg);
    SDL_DestroyRenderer(ren);
    SDL_DestroyWindow(win);
    IMG_Quit();
    SDL_Quit();

    return 0;
}
