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

// ==================== FFT =============================

// Taille de la FFT (puissance de 2 recommandée pour la perf.)
// Utilisée pour la taille des buffers d'entrée/sortie FFT.
#define FFT_SIZE           16384

// Taille de la sortie FFT réelle→complexe (1/2 N + 1)
// Pour le buffer de magnitude FFT.
#define FFT_SIZE_b2        ((FFT_SIZE/2)+1)

// Fréquence représentée par chaque bin FFT (pour conversion bin ↔ Hz)
const float fft_freq_per_bin = 48000.0f / FFT_SIZE;

// Fenêtre de Hann (fenêtrage spectral pour la FFT)
float fft_hann_window[FFT_SIZE];

// Plan FFTW (objet optimisé pour la transformation, créé au lancement)
fftwf_plan fft_plan;

// Buffer d'entrée (mono) pour la FFT
float fft_input[FFT_SIZE];

// Buffer de sortie complexe (spectre fréquence)
fftwf_complex fft_output[FFT_SIZE_b2];

// Mutex & condition : synchronisation thread audio/FFT
pthread_mutex_t fft_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t  fft_cond  = PTHREAD_COND_INITIALIZER;
volatile int    fft_block_ready = 0;   // Flag : nouveau bloc FFT dispo

// ==================== AUDIO / ALSA ====================

// Taille du buffer audio lu par ALSA à chaque itération
#define AUDIO_BUFFER_SIZE   HOP_SIZE
#define HOP_SIZE           1024      // Nombre d'échantillons par "hop"

// Taux d'échantillonnage (Hz) : utilisé pour les conversions temporelles/fréquences
const float audio_sample_rate = 48000.0f;

// Buffer circulaire (stéréo intercalé) : stockage temporaire du flux audio
short audio_ringbuffer[FFT_SIZE * 2];
int   audio_ringbuffer_pos = 0;

// Buffer temporaire pour le thread FFT (copie du ringbuffer)
short audio_fft_samples[FFT_SIZE * 2];

// ==================== SPECTROGRAMME ====================

// Nombre de bandes de fréquence pour l'affichage du spectrogramme
#define NUM_COLS           31

// Fréquences centrales de chaque bande (Hz)
const float spectro_band_freqs[NUM_COLS] = {    20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250,
    315, 400, 500, 630, 800, 1000, 1250, 1600, 2000,
    2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500, 16000, 20000
 };

// Bornes de chaque bande (Hz): pour calculer l’énergie par bande
const float spectro_band_edges[NUM_COLS+1] = {     
	16.00, 22.36, 28.07, 35.48, 44.72, 56.12, 71.05, 89.44,
    111.80, 141.42, 179.60, 223.61, 279.51, 354.81, 447.21, 561.23,
    707.11, 894.43, 1118.03, 1409.49, 1788.85, 2236.07, 2795.08,
    3548.13, 4472.14, 5612.49, 7071.07, 8944.27, 11180.34,
    14142.14, 17888.54, 22400.00
 };

// Calibration "spectrale" (modifiable à l'exécution, sauvegardée)
// Normalise la réponse par bande (utilisé pour compenser la fenêtre, le micro, etc.)
static float spectro_ref_val[NUM_COLS] = {     53.03, 37.27, 41.88, 34.97, 55.34, 39.77, 19.11, 35.20, 40.41, 38.98,
    38.21, 46.98, 45.60, 41.15, 39.22, 39.04, 36.95, 40.23, 40.43, 37.16,
    41.96, 40.30, 39.98, 39.72, 38.29, 39.82, 40.04, 39.42, 41.26, 40.13, 41.48
 };

// Flag mode calibration, et buffers d’état pour le spectrogramme (mis à jour en temps réel)
bool   spectro_calibrate       = 0;
float  spectro_band_db[NUM_COLS]      = {0.0f};  // Valeur RMS/dB par bande
float  spectro_peak_hold[NUM_COLS]    = {0};     // Pic visuel "peak hold"
float  spectro_peak_time[NUM_COLS]    = {0};     // Dernier timestamp du peak
float  spectro_tau_peak_hold   = 0.5f; // Temps de persistance du peak hold (s)
float  spectro_tau_peak_fall   = 2.0f; // Temps de décroissance du peak hold (s)

// Coordonnées graphiques pour l’affichage des barres du spectrogramme
const int spectro_bar_x[NUM_COLS] = {     120,160,200,240,280,320,360,400,440,480,520,560,600,640,680,720,
    760,800,840,880,920,960,1000,1040,1080,1120,1160,1200,1240,1280,1320
 };
const int spectro_bar_y = 192;
const int spectro_bar_w = 22;
const int spectro_bar_h = 800;

// ==================== LSR (LEFT/SIDE/RIGHT) ================

// Coordonnées et tailles pour affichage des barres L/S/R (vu-mètre stéréo + médiane)
const int lsr_bar_x[3]   = { 1482, 1527, 1558 };
const int lsr_bar_w[3]   = { 36, 22, 36 };
const int lsr_bar_y = 192;
const int lsr_bar_h = 800;

// Buffers de niveaux audio (L, S, R) : valeurs normalisées pour affichage
float lsr_rms[3]         = {0.0f};      // RMS gauche, side, droite
float lsr_peak[3]        = {0.0f};      // Peak instantané L, S, R
float lsr_balance        = 0.0f;        // Balance gauche/droite
static float lsr_phase_rho_smoothed = 0.0f;   /* valeur lissée pour l’affichage */
const  float lsr_phase_tau          = 0.25f;  /* τ = 250 ms : réactif mais lisible */

// Libellés d'affichage et calibration vumètre (modifiable à l'exécution)
const char *lsr_labels[3]= {"L", "S", "R"};
float lsr_ref_vu_L       = 1.0f;
float lsr_ref_vu_R       = 1.0f;
bool  lsr_calibrate      = 0;           // Mode calibration vumètre

// Paramètres de lissage/relâchement pour l’affichage
float lsr_tau_rms        = 0.1f;        // Lissage RMS (s)
float lsr_tau_peak_release = 0.3f;      // Relâchement peak (s)

// Coordonnées et tailles pour affichage des barres diff L/R
const int lsr_diffLR_left_x  = 1482, lsr_diffLR_left_w  = 56;
const int lsr_diffLR_right_x = 1538, lsr_diffLR_right_w = 56;
const int lsr_diffLR_bar_y   = 175,  lsr_diffLR_bar_h   = 11;

// Coordonnées et tailles pour affichage des barres phase L/R
const int lsr_phaseLR_bar_x = 1482, lsr_phaseLR_bar_y = 164;
const int lsr_phaseLR_bar_w = 112,  lsr_phaseLR_bar_h = 11;

/* ===============   Hardware  ============================ */
static float attenuator_db = 0;

/* ===============   OSCILLO L / R  ============================ */
const int osc_L_x = 670, osc_L_y = 30;     /* cadre haut (L)  – SVG */
const int osc_R_x = 670, osc_R_y = 75;     /* cadre bas  (R)  – SVG */
const int osc_w  = 672,  osc_h  = 45;

static float osc_buf[HOP_SIZE*2];
static int   osc_pos = 0;                /* index circulaire courant */

// ==================== UI / AFFICHAGE =======================

// Nombre de repères dB affichés verticalement (graduation dB sur les barres)
#define ui_num_scale      9

// Valeurs dB et positions y correspondantes pour l'affichage
const float ui_db_mark[ui_num_scale] = { +6, 0, -3, -8, -14, -20, -30, -40, -80 };
const int   ui_y_mark[ui_num_scale]  = { 192, 264, 308, 392, 496, 592, 752, 888, 992 };

// Définition d'une couleur RGB (pour gradient/barres)
typedef struct { Uint8 r, g, b; } Color;

// Points d'arrêt du gradient de couleur pour le niveau audio
Color ui_gradient_stops[8] = {
    {0, 0, 0},      {0, 0, 128},    {0, 0, 255},    {0, 128, 0},
    {0, 255, 0},    {255, 255, 0},  {255, 128, 0},  {255, 0, 0}
};
static SDL_Texture *ui_gradient_texture = NULL; // Texture verticale du gradient
Color ui_gradient_LUT[800];                    // Lookup-table gradient en RAM

// Position/tailles graphiques pour cercles indicateurs (niveau, clipping)
const int ui_circle_cy     = 172;
const int ui_circle_radius = 11;

// ==================== APP =============================

// Flag d’exécution principale (set à 0 pour arrêter l’app proprement)
volatile sig_atomic_t app_running = 1;



// ====== Prototypes =============================================================================================
void handle_sigint(int sig);
SDL_Texture *loadPngFromMemory(SDL_Renderer *renderer);
Color getGradientColor(float v);
int circle_cx(int i);
void getCircleColor(float v, Uint8 *r, Uint8 *g, Uint8 *b);
void drawFilledCircle(SDL_Renderer *renderer, int x0, int y0, int radius, Uint8 r, Uint8 g, Uint8 b, Uint8 a);
void getGradLeftColor(float t, Uint8 *r, Uint8 *g, Uint8 *b);
void getGradRightColor(float t, Uint8 *r, Uint8 *g, Uint8 *b);
void draw_spectrogram_bars(SDL_Renderer *ren, const float *spectro_band_db);
void draw_lr_bars(SDL_Renderer *ren, const float *rms_values, const float *peak_values);
void draw_balance_bar(SDL_Renderer *ren, float value);
snd_pcm_t *init_alsa_capture(const char *device);
void audio_update_levels(snd_pcm_t *handle, float *lsr_rms, float *lsr_balance);
void* fft_thread_func(void *arg);
void save_calibration_file();
int load_calibration_file();
void init_hann_window();
void init_gradientLUT();
int calc_filled_height(float val);

// ====== Fonctions =============================================================================================

// -- Handler signal Ctrl-C --
void handle_sigint(int sig) {
    app_running = 0;
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
    if (v <= 0) return ui_gradient_stops[0];
    if (v >= 1) return ui_gradient_stops[7];
    float pos = v * 7.0f;
    int idx = (int)pos;
    float t = pos - idx;
    Color a = ui_gradient_stops[idx];
    Color b = ui_gradient_stops[idx + 1];
    Color c;
    c.r = a.r + (b.r - a.r) * t;
    c.g = a.g + (b.g - a.g) * t;
    c.b = a.b + (b.b - a.b) * t;
    return c;
}

// -- Coordonnées en X des cercles indicateurs --
int circle_cx(int i) { return spectro_bar_x[i] + spectro_bar_w / 2; }

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
    if (db >= ui_db_mark[0]) return ui_y_mark[0];
    if (db <= ui_db_mark[ui_num_scale-1]) return ui_y_mark[ui_num_scale-1];
    for (int i = 0; i < ui_num_scale-1; ++i) {
        if (db <= ui_db_mark[i] && db > ui_db_mark[i+1]) {
            float t = (db - ui_db_mark[i+1]) / (ui_db_mark[i] - ui_db_mark[i+1]);
            return (int)(ui_y_mark[i+1] + t * (ui_y_mark[i] - ui_y_mark[i+1]));
        }
    }
    return ui_y_mark[ui_num_scale-1];
}

// -- Dessine les barres du spectrogramme --
void draw_spectrogram_bars(SDL_Renderer *ren, const float *spectro_band_db) {
    for (int i = 0; i < NUM_COLS; ++i) {
        SDL_Rect rc = {spectro_bar_x[i], spectro_bar_y, spectro_bar_w, spectro_bar_h};
        SDL_SetRenderDrawColor(ren, 255, 255, 255, 255);
        SDL_RenderDrawRect(ren, &rc);

        int y_head = db_to_y(spectro_band_db[i]);
        int y_bottom = spectro_bar_y + spectro_bar_h - 1;
        // Remplissage du bas vers le haut jusqu'à y_head
        int filled_h = (spectro_bar_y + spectro_bar_h) - y_head;
        if (filled_h > 0) {
            SDL_Rect src = {0, 0, 1, filled_h}; // on prend le haut (0…filled_h-1)
            SDL_Rect dst = {spectro_bar_x[i] + 1, y_head, spectro_bar_w - 2, filled_h};
            SDL_RenderCopyEx(ren, ui_gradient_texture, &src, &dst, 0.0, NULL, SDL_FLIP_VERTICAL);
        }

        // Ligne du peak hold (petit trait blanc semi-transparent)
        int y_peak = db_to_y(spectro_peak_hold[i]);
        if (y_peak < spectro_bar_y) y_peak = spectro_bar_y;
        if (y_peak > y_bottom) y_peak = y_bottom;
        SDL_SetRenderDrawColor(ren, 255, 255, 255, 180);
        SDL_Rect peak_rect = {spectro_bar_x[i] + 2, y_peak, spectro_bar_w - 4, 2};
        SDL_RenderFillRect(ren, &peak_rect);

        // Cercles indicateurs de niveau
        Uint8 r, g, b;
        float v = (spectro_band_db[i] + 60.0f) / 66.0f;
        getCircleColor(v, &r, &g, &b);
        int cx = circle_cx(i);
        int cy = ui_circle_cy;
        drawFilledCircle(ren, cx, cy, ui_circle_radius, r, g, b, 255);
    }
}


// -- Calcul hauteur remplie pour barres RMS ou Peak (en pixels) --
int calc_filled_height(float val) {
    float val_db = 20.0f * log10f(fmaxf(val, 1e-5f)) + attenuator_db;
    int y_pos = db_to_y(val_db);
    int filled_h = 0;
    if (y_pos < lsr_bar_y) {
        filled_h = lsr_bar_h;
    } else if (y_pos > lsr_bar_y + lsr_bar_h) {
        filled_h = 0;
    } else {
        filled_h = (lsr_bar_y + lsr_bar_h) - y_pos;
    }
    if (filled_h < 0) filled_h = 0;
    if (filled_h > lsr_bar_h) filled_h = lsr_bar_h;
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
        SDL_Rect rc = { lsr_bar_x[i], lsr_bar_y, lsr_bar_w[i], lsr_bar_h };
        SDL_SetRenderDrawColor(ren, 255, 255, 255, 255);
        SDL_RenderDrawRect(ren, &rc);

        /* ------------------ RMS plein ------------------ */
        int filled = calc_filled_height(rms_values[i]);
        if (filled) {
            SDL_Rect src = { 0, 0, 1, filled };           /* bas de la texture */
            SDL_Rect dst = { lsr_bar_x[i] + 1,
                             lsr_bar_y + lsr_bar_h - filled,
                             lsr_bar_w[i] - 2,
                             filled };
            SDL_SetTextureAlphaMod(ui_gradient_texture, 255);          /* opaque */
            SDL_RenderCopyEx(ren, ui_gradient_texture, &src, &dst,
                             0.0, NULL, SDL_FLIP_VERTICAL); /* noir en bas */
        }

        /* ------------------ Peak semi-transparent ------------------ */
        int filled_peak = calc_filled_height(peak_values[i]);
        if (filled_peak) {
            SDL_Rect src = { 0, 0, 1, filled_peak };
            SDL_Rect dst = { lsr_bar_x[i] + 1,
                             lsr_bar_y + lsr_bar_h - filled_peak,
                             lsr_bar_w[i] - 2,
                             filled_peak };
            SDL_SetTextureAlphaMod(ui_gradient_texture, 128);          /* 50 % */
            SDL_RenderCopyEx(ren, ui_gradient_texture, &src, &dst,
                             0.0, NULL, SDL_FLIP_VERTICAL);
        }
    }
    /* Remise de l’opacité pour la suite (spectrogramme) */
    SDL_SetTextureAlphaMod(ui_gradient_texture, 255);
}


// -- Dessine la barre de balance L/R avec dégradé coloré --
void draw_balance_bar(SDL_Renderer *ren, float value) {
    if (value < -1.0f) value = -1.0f;
    if (value >  1.0f) value =  1.0f;

    if (value < 0.0f) {
        int fill = (int)(lsr_diffLR_left_w * (-value));
        for (int dx = 0; dx < fill; ++dx) {
            float t = (lsr_diffLR_left_w - 1 - dx) / (float)(lsr_diffLR_left_w - 1);
            Uint8 r, g, b;
            getGradLeftColor(t, &r, &g, &b);
            SDL_SetRenderDrawColor(ren, r, g, b, 255);
            SDL_RenderDrawLine(ren, lsr_diffLR_left_x + (lsr_diffLR_left_w - 1 - dx), lsr_diffLR_bar_y, lsr_diffLR_left_x + (lsr_diffLR_left_w - 1 - dx), lsr_diffLR_bar_y + lsr_diffLR_bar_h - 1);
        }
    }
    else if (value > 0.0f) {
        int fill = (int)(lsr_diffLR_right_w * value);
        for (int dx = 0; dx < fill; ++dx) {
            float t = dx / (float)(lsr_diffLR_right_w - 1);
            Uint8 r, g, b;
            getGradRightColor(t, &r, &g, &b);
            SDL_SetRenderDrawColor(ren, r, g, b, 255);
            SDL_RenderDrawLine(ren, lsr_diffLR_right_x + dx, lsr_diffLR_bar_y, lsr_diffLR_right_x + dx, lsr_diffLR_bar_y + lsr_diffLR_bar_h - 1);
        }
    }
}

/* === Barre de corrélation –1 … +1 ====================== */
void draw_phase_bar(SDL_Renderer *ren, float rho)
{
    /* curseur couleur */
    Uint8 r,g,b;
    if (rho < 0.0f)      { r=255; g=0;   b=0;   }
    else if (rho < 0.5f) { r=255; g=140; b=0;   }
    else                 { r=0;   g=255; b=0;   }

    int cx = lsr_phaseLR_bar_x + (int)((rho + 1.0f) * 0.5f * lsr_phaseLR_bar_w);
    SDL_SetRenderDrawColor(ren, r, g, b, 255);
    SDL_Rect cur = {cx-2, lsr_phaseLR_bar_y-2, 4, lsr_phaseLR_bar_h+4};
    SDL_RenderFillRect(ren, &cur);
}

static void draw_oscilloscopes(SDL_Renderer *ren, const float *buf)
{
    float y_gain = 1.0f;
    int trigger = -1;

    // Cherche le front montant pour L
    for (int i = 1; i < HOP_SIZE; ++i) {
        float prev = buf[((osc_pos + i - 1) % HOP_SIZE) * 2];
        float curr = buf[((osc_pos + i    ) % HOP_SIZE) * 2];
        if (prev < 0.0f && curr >= 0.0f) {
            trigger = i;
            break;
        }
    }
    if (trigger == -1) trigger = 0;

    // ---- Oscilloscope L (vert) ----
    SDL_SetRenderDrawColor(ren, 80, 255, 80, 255); // vert
    for (int i = 1; i < HOP_SIZE; ++i) {
        int idx0 = (osc_pos + trigger + i - 1) % HOP_SIZE;
        int idx1 = (osc_pos + trigger + i    ) % HOP_SIZE;
        float v0 = buf[idx0 * 2];
        float v1 = buf[idx1 * 2];
        int px0 = osc_L_x + (i-1) * osc_w / (HOP_SIZE-1);
        int py0 = osc_L_y + osc_h/2 - (int)(v0 * (osc_h/2) * y_gain);
        int px1 = osc_L_x + (i  ) * osc_w / (HOP_SIZE-1);
        int py1 = osc_L_y + osc_h/2 - (int)(v1 * (osc_h/2) * y_gain);
        SDL_RenderDrawLine(ren, px0, py0, px1, py1);
    }

    // ---- Oscilloscope R (bleu) ----
    SDL_SetRenderDrawColor(ren, 80, 200, 255, 255); // bleu
    for (int i = 1; i < HOP_SIZE; ++i) {
        int idx0 = (osc_pos + trigger + i - 1) % HOP_SIZE;
        int idx1 = (osc_pos + trigger + i    ) % HOP_SIZE;
        float v0 = buf[idx0 * 2 + 1];
        float v1 = buf[idx1 * 2 + 1];
        int px0 = osc_R_x + (i-1) * osc_w / (HOP_SIZE-1);
        int py0 = osc_R_y + osc_h/2 - (int)(v0 * (osc_h/2) * y_gain);
        int px1 = osc_R_x + (i  ) * osc_w / (HOP_SIZE-1);
        int py1 = osc_R_y + osc_h/2 - (int)(v1 * (osc_h/2) * y_gain);
        SDL_RenderDrawLine(ren, px0, py0, px1, py1);
    }
}


 



// -- Initialisation de la fenêtre de Hann pour FFT --
void init_hann_window() {
    for (int i = 0; i < FFT_SIZE; i++)
        fft_hann_window[i] = 0.5f * (1.0f - cosf(2.0f * M_PI * i / (FFT_SIZE - 1)));
}

// -- Initialisation du LUT dégradé vertical --
// --- Création du LUT EN MÉMOIRE + texture GPU en même temps ---
void init_gradient_texture(SDL_Renderer *ren)
{
    // 1. Surface RGBA 1 x spectro_bar_h
    SDL_Surface *surf = SDL_CreateRGBSurfaceWithFormat(
        0, 1, spectro_bar_h, 32, SDL_PIXELFORMAT_RGBA8888);
    if (!surf) {
        fprintf(stderr, "SDL_CreateRGBSurface error: %s\n", SDL_GetError());
        exit(1);
    }

    Uint32 *pix = (Uint32 *)surf->pixels;
    for (int dy = 0; dy < spectro_bar_h; ++dy) {
        float val = (float)dy / (spectro_bar_h - 1);
        Color c = getGradientColor(val);          // ta table 8 couleurs
        pix[dy] = SDL_MapRGBA(surf->format, c.r, c.g, c.b, 255);
        ui_gradient_LUT[dy] = c;                      // ← garde aussi le LUT si d’autres fonctions l’utilisent
    }

    // 2. Texture GPU
    ui_gradient_texture = SDL_CreateTextureFromSurface(ren, surf);
    SDL_FreeSurface(surf);
    if (!ui_gradient_texture) {
        fprintf(stderr, "SDL_CreateTextureFromSurface error: %s\n", SDL_GetError());
        exit(1);
    }
    SDL_SetTextureBlendMode(ui_gradient_texture, SDL_BLENDMODE_BLEND);  // autorise alpha pour le peak
}


// -- Initialisation ALSA capture (input audio) --
snd_pcm_t *init_alsa_capture(const char *device) {
    snd_pcm_t *alsa_handle = NULL;
    snd_pcm_hw_params_t *hw_params = NULL;
    int err;

    if ((err = snd_pcm_open(&alsa_handle, device, SND_PCM_STREAM_CAPTURE, 0)) < 0) {
        fprintf(stderr, "ALSA open error (%s): %s\n", device, snd_strerror(err));
        return NULL;
    }

    snd_pcm_hw_params_alloca(&hw_params);
    if ((err = snd_pcm_hw_params_any(alsa_handle, hw_params)) < 0) {
        fprintf(stderr, "ALSA: Impossible d'obtenir les paramètres hardware: %s\n", snd_strerror(err));
        goto fail;
    }
    if ((err = snd_pcm_hw_params_set_format(alsa_handle, hw_params, SND_PCM_FORMAT_S16_LE)) < 0) {
        fprintf(stderr, "ALSA: Impossible de définir le format S16_LE: %s\n", snd_strerror(err));
        goto fail;
    }
    if ((err = snd_pcm_hw_params_set_channels(alsa_handle, hw_params, 2)) < 0) {
        fprintf(stderr, "ALSA: Impossible de définir 2 canaux: %s\n", snd_strerror(err));
        goto fail;
    }
    unsigned int rate = (unsigned int)audio_sample_rate;
    if ((err = snd_pcm_hw_params_set_rate_near(alsa_handle, hw_params, &rate, 0)) < 0) {
        fprintf(stderr, "ALSA: Impossible de définir la fréquence d'échantillonnage: %s\n", snd_strerror(err));
        goto fail;
    }
    if (rate != 48000)
        fprintf(stderr, "ALSA: Fréquence obtenue: %u Hz (au lieu de 48000 Hz)\n", rate);

    if ((err = snd_pcm_hw_params_set_access(alsa_handle, hw_params, SND_PCM_ACCESS_RW_INTERLEAVED)) < 0) {
        fprintf(stderr, "ALSA: Impossible de définir l'accès interleaved: %s\n", snd_strerror(err));
        goto fail;
    }

    snd_pcm_uframes_t buffer_size = AUDIO_BUFFER_SIZE;
    snd_pcm_hw_params_set_period_size_near(alsa_handle, hw_params, &buffer_size, 0);

    if ((err = snd_pcm_hw_params(alsa_handle, hw_params)) < 0) {
        fprintf(stderr, "ALSA: Impossible d'appliquer les paramètres hardware: %s\n", snd_strerror(err));
        goto fail;
    }

    if ((err = snd_pcm_prepare(alsa_handle)) < 0) {
        fprintf(stderr, "ALSA: snd_pcm_prepare() a échoué: %s\n", snd_strerror(err));
        goto fail;
    }

    fprintf(stderr, "ALSA: Capture initialisée sur %s\n", device);
    return alsa_handle;

fail:
    if (alsa_handle)
        snd_pcm_close(alsa_handle);
    return NULL;
}

// -- Mise à jour des niveaux audio RMS, Peak et balance --
void audio_update_levels(snd_pcm_t *handle, float *lsr_rms, float *lsr_balance) {
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
            lsr_rms[0] = lsr_rms[1] = lsr_rms[2] = 0.0f;
            lsr_peak[0] = lsr_peak[1] = lsr_peak[2] = 0.0f;
            *lsr_balance = 0.0f;
            return;
        }
        break;
    }

    // Calcul RMS et peak bruts
	double dot = 0.0;
	double sum_l = 0.0, sum_r = 0.0; // sert pour RMS ET corrélation !
	float peak_l = 0.0f, peak_r = 0.0f;

	for (int i = 0; i < nframes; ++i) {
		float l = buffer[2*i]   / 32768.0f;
		float r = buffer[2*i+1] / 32768.0f;
		// Corrélation
		dot    += l * r;
		sum_l  += l * l;
		sum_r  += r * r;
		// Peak
		float abs_l = fabsf(l);
		float abs_r = fabsf(r);
		if (abs_l > peak_l) peak_l = abs_l;
		if (abs_r > peak_r) peak_r = abs_r;
		// Oscilloscope buffer
		osc_buf[osc_pos*2]   = l;
		osc_buf[osc_pos*2+1] = r;
		osc_pos              = (osc_pos + 1) & (HOP_SIZE-1); // si HOP_SIZE = puissance de 2
	}

	float rho = 0.0f;
	if (sum_l > 1e-9 && sum_r > 1e-9) // on remplace energyL/energyR
		rho = dot / sqrtf(sum_l * sum_r);

	float rms_l = sqrt(sum_l / nframes);
	float rms_r = sqrt(sum_r / nframes);


    // Gestion du temps
    float cur_time = SDL_GetTicks() * 0.001f; // secondes
    float dt = cur_time - last_time;
    if (dt <= 0.0f) dt = 0.001f; // éviter div0 ou négatif

    float alpha_corr = 1.0f - expf(-dt / lsr_phase_tau);
    lsr_phase_rho_smoothed = alpha_corr * rho + (1.0f - alpha_corr) * lsr_phase_rho_smoothed;

    // Calibration vumètre (idem que dans ta version)

    if (lsr_calibrate) {
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
                    lsr_ref_vu_L = stable_rms_l;
                    lsr_ref_vu_R = stable_rms_r;
                    printf("[CALIB VU] lsr_ref_vu_L = %.6f, lsr_ref_vu_R = %.6f\n", lsr_ref_vu_L, lsr_ref_vu_R);
                    lsr_calibrate = 0;
                    save_calibration_file();
                    waiting_for_stability = 0;
                }
            }
        } else {
            waiting_for_stability = 0;
        }
    }

    // Calcul alpha
    float alpha_rms = 1.0f - expf(-dt / lsr_tau_rms);

    // Lissage RMS
    smoothed_rms_l = alpha_rms * rms_l + (1.0f - alpha_rms) * smoothed_rms_l;
    smoothed_rms_r = alpha_rms * rms_r + (1.0f - alpha_rms) * smoothed_rms_r;

    // Lissage Peak, attaque instantanée, relâchement lissé
    if (peak_l > smoothed_peak_l) {
        smoothed_peak_l = peak_l;
    } else {
        float alpha_peak = 1.0f - expf(-dt / lsr_tau_peak_release);
        smoothed_peak_l = alpha_peak * peak_l + (1.0f - alpha_peak) * smoothed_peak_l;
    }

    if (peak_r > smoothed_peak_r) {
        smoothed_peak_r = peak_r;
    } else {
        float alpha_peak = 1.0f - expf(-dt / lsr_tau_peak_release);
        smoothed_peak_r = alpha_peak * peak_r + (1.0f - alpha_peak) * smoothed_peak_r;
    }

    last_time = cur_time;

    // Normalisation RMS corrigée avec calibration
    float norm_l = (smoothed_rms_l / lsr_ref_vu_L) * 0.70710678f;
    float norm_r = (smoothed_rms_r / lsr_ref_vu_R) * 0.70710678f;
    float norm_s = (norm_l + norm_r) * 0.5f;

    if (norm_l > 1.0f) norm_l = 1.0f;
    if (norm_r > 1.0f) norm_r = 1.0f;
    if (norm_s > 1.0f) norm_s = 1.0f;

    lsr_rms[0] = norm_l;
    lsr_rms[2] = norm_r;
    lsr_rms[1] = norm_s;

    // Clamp peaks
    if (smoothed_peak_l > 1.0f) smoothed_peak_l = 1.0f;
    if (smoothed_peak_r > 1.0f) smoothed_peak_r = 1.0f;
    float peak_s = (smoothed_peak_l + smoothed_peak_r) * 0.5f;
    if (peak_s > 1.0f) peak_s = 1.0f;

    lsr_peak[0] = smoothed_peak_l;
    lsr_peak[2] = smoothed_peak_r;
    lsr_peak[1] = peak_s;

    // Calcul balance
    float total = rms_l + rms_r;
    if (total > 1e-6)
        *lsr_balance = (rms_r - rms_l) / total;
    else
        *lsr_balance = 0.0f;

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
        audio_fft_samples[i*2 + 0] = audio_ringbuffer[src_pos*2 + 0];
        audio_fft_samples[i*2 + 1] = audio_ringbuffer[src_pos*2 + 1];
    }

    fft_block_ready = 1;
    pthread_cond_signal(&fft_cond);
    pthread_mutex_unlock(&fft_mutex);
}


// -- Thread FFT : calcul de la FFT et calibration spectre --
void* fft_thread_func(void *arg) {
    short block_stereo[FFT_SIZE * 2];
    static int last_band = -1;
    static float last_time = 0.0f;
    static int already_calibrated = -1;
    static float last_peak_update_time = 0.0f;

    while (app_running) {
        pthread_mutex_lock(&fft_mutex);
        while (!fft_block_ready && app_running)
            pthread_cond_wait(&fft_cond, &fft_mutex);
        if (!app_running) {
            pthread_mutex_unlock(&fft_mutex);
            break;
        }
        memcpy(block_stereo, audio_fft_samples, sizeof(short) * 2 * FFT_SIZE);
        fft_block_ready = 0;
        pthread_mutex_unlock(&fft_mutex);

        // Conversion stéréo vers mono et application fenêtre Hann
        for (int i = 0; i < FFT_SIZE; i++) {
            float L = block_stereo[2*i] / 32768.0f;
            float R = block_stereo[2*i+1] / 32768.0f;
            fft_input[i] = 0.5f * (L + R) * fft_hann_window[i];
        }
        fftwf_execute_dft_r2c(fft_plan, fft_input, fft_output);
        static float mag[FFT_SIZE_b2];
        for (int i = 0; i < FFT_SIZE_b2; i++) mag[i] = hypotf(fft_output[i][0], fft_output[i][1]);

        // Calcul niveaux moyens par bande de fréquence
        for (int c = 0; c < NUM_COLS; c++) {
            float f_low  = spectro_band_edges[c];
            float f_high = spectro_band_edges[c+1];
            int bin_start = (int)ceilf(f_low / fft_freq_per_bin);
            int bin_end   = (int)floorf(f_high / fft_freq_per_bin);
            if (bin_start < 0) bin_start = 0;
            if (bin_end > FFT_SIZE_b2) bin_end = FFT_SIZE_b2;
            float sum = 0.0f;
            int n_bins = bin_end - bin_start + 1;
            for (int i = bin_start; i <= bin_end; i++) sum += mag[i];
            float avg = (n_bins > 0) ? (sum / n_bins) : 1e-20f;
            float avg_norm = avg / spectro_ref_val[c];
            float db = 20.0f * log10f(avg_norm + 1e-20f) + attenuator_db;
            spectro_band_db[c] = db;
        }

        // === GESTION PEAK HOLD PRO ===
        float cur_time = SDL_GetTicks() * 0.001f; // secondes
        float dt = cur_time - last_peak_update_time;
        if (dt <= 0.0f) dt = 0.001f;

        for (int c = 0; c < NUM_COLS; c++) {
            float current = spectro_band_db[c];
            // Nouveau pic : update peak + timestamp
            if (current > spectro_peak_hold[c]) {
                spectro_peak_hold[c] = current;
                spectro_peak_time[c] = cur_time;
            } else {
                // Si HOLD pas fini : freeze
                if ((cur_time - spectro_peak_time[c]) < spectro_tau_peak_hold) {
                    // Rien à faire, on fige la valeur du peak
                } else {
                    // Après HOLD : décroissance linéaire
                    float db_down_per_sec = 20.0f / spectro_tau_peak_fall;
                    float db_drop = db_down_per_sec * dt;
                    float new_peak = spectro_peak_hold[c] - db_drop;
                    // Ne pas descendre sous le niveau actuel
                    if (new_peak < current)
                        new_peak = current;
                    spectro_peak_hold[c] = new_peak;
                }
            }
        }
        last_peak_update_time = cur_time;

        // Calibration spectrogramme simplifiée (1s stable sur la même bande)
        if (spectro_calibrate) {
            // Réutiliser cur_time défini ci-dessus
            // Trouver la bande la plus forte
            int best_band = -1;
            float max_db = -1000.0f;
            for (int c = 0; c < NUM_COLS; ++c) {
                if (spectro_band_db[c] > max_db) {
                    max_db = spectro_band_db[c];
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
                float f_low  = spectro_band_edges[best_band];
                float f_high = spectro_band_edges[best_band+1];
                int bin_start = (int)ceilf(f_low / fft_freq_per_bin);
                int bin_end   = (int)floorf(f_high / fft_freq_per_bin);
                float sum = 0.0f;
                int n_bins = bin_end - bin_start + 1;
                for (int i = bin_start; i <= bin_end; i++) sum += mag[i];
                float avg = (n_bins > 0) ? (sum / n_bins) : 1e-20f;
                printf("[CALIB] Bande %d (%.1f Hz) : spectro_ref_val[%d] = %.2f (avant: %.2f)\n",
                       best_band, spectro_band_freqs[best_band], best_band, avg, spectro_ref_val[best_band]);
                spectro_ref_val[best_band] = avg;
                if (best_band + 1 < NUM_COLS) {
                    spectro_ref_val[best_band + 1] = avg;
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
        fprintf(f, "%.2f", spectro_ref_val[i]);
        if (i < NUM_COLS - 1)
            fprintf(f, ";");
    }
    fprintf(f, "\n");

    fprintf(f, "%.6f;%.6f\n", lsr_ref_vu_L, lsr_ref_vu_R);

    fclose(f);
}

fftwf_plan create_fftwf_plan_with_wisdom(int size, float* fft_input, fftwf_complex* fft_output) {
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
    plan = fftwf_plan_dft_r2c_1d(size, fft_input, fft_output, FFTW_MEASURE);
    if (!plan) {
        fprintf(stderr, "[FFTW] Erreur : impossible de créer le plan FFTW !\n");
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
        spectro_ref_val[i] = strtof(token, NULL);
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

    lsr_ref_vu_L = strtof(token_L, NULL);
    lsr_ref_vu_R = strtof(token_R, NULL);

    fclose(f);
    return 0;
}

void print_help(const char *progname) {
    printf("Usage: %s [options]\n", progname);
    printf("Options disponibles :\n");
    printf("  -D <device>                 : Périphérique ALSA d'acquisition audio (défaut : 'default')\n");
	printf("  -AD <attenuator_db>         : Compensation atténuateur dB (ex : 9 pour +9 dB)\n");
    printf("  -CS                         : Active la calibration du spectrogramme\n");
    printf("     Utilisation :\n");
    printf("       freqs=(20 25 31.5 40 50 63 80 100 125 160 200 250 315 400 500 630 800 1000 1250 1600 2000 2500 3150 4000 5000 6300 8000 10000 12500 16000 20000)\n");
    printf("       for f in \"${freqs[@]}\"; do\n");
    printf("         sox -n -c 2 -t alsa plughw:CARD=loopTest1,DEV=1 synth 5 sine $f vol 0dB\n");
    printf("       done\n");
    printf("  -CV                         : Active la calibration du vumètre L/R\n");
    printf("     Utilisation :\n");
    printf("       sox -n -t alsa -b 16 -e signed-integer -c 2 plughw:CARD=loopTest1,DEV=1 synth 60 sine 1000 vol 0dB\n");
    printf("  -PT <tau_peak_hold>         : Durée de persistance du peak hold du spectrogramme (en s, défaut : 0.5)\n");
    printf("  -PF <tau_peak_fall>         : Temps de décroissance du peak hold du spectrogramme (en s, défaut : 2.0)\n");
    printf("  -TR <lsr_tau_rms>           : Constante de temps pour le lissage RMS (en s, défaut : 0.1)\n");
    printf("  -TP <lsr_tau_peak_release>  : Constante de temps pour la relâche du peak (en s, défaut : 0.3)\n");
    printf("  -h, --help                  : Affiche ce message d'aide\n");
    printf("\n");
}

// ====== Fonction principale ======
int main(int argc, char *argv[]) {
    const char *alsa_device = "default";
    // Parsing des arguments ligne de commande
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-D") == 0 && i + 1 < argc) {
            alsa_device = argv[++i];
        }
        else if (strcmp(argv[i], "-CS") == 0) {
            spectro_calibrate = 1;
        }
        else if (strcmp(argv[i], "-CV") == 0) {
            lsr_calibrate = 1;
        }
        else if (strcmp(argv[i], "-PT") == 0 && i + 1 < argc) {
            spectro_tau_peak_hold = atof(argv[++i]);
            if (spectro_tau_peak_hold < 0.1f) spectro_tau_peak_hold = 0.1f;
        }
        else if (strcmp(argv[i], "-PF") == 0 && i + 1 < argc) {
            spectro_tau_peak_fall = atof(argv[++i]);
            if (spectro_tau_peak_fall < 0.1f) spectro_tau_peak_fall = 0.1f;
        }
        else if (strcmp(argv[i], "-TR") == 0 && i + 1 < argc) {
            lsr_tau_rms = atof(argv[++i]);
            if (lsr_tau_rms < 0.01f) lsr_tau_rms = 0.01f;
        }
        else if (strcmp(argv[i], "-TP") == 0 && i + 1 < argc) {
            lsr_tau_peak_release = atof(argv[++i]);
            if (lsr_tau_peak_release < 0.01f) lsr_tau_peak_release = 0.01f;
        }
		else if (strcmp(argv[i], "-AD") == 0 && i + 1 < argc) {
			attenuator_db = atof(argv[++i]);
		}
        else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_help(argv[0]);
            return 0;
        }
        else {
            fprintf(stderr, "Option inconnue : %s\n", argv[i]);
            print_help(argv[0]);
            return 1;
        }
    }

	
    fft_plan = create_fftwf_plan_with_wisdom(FFT_SIZE, fft_input, fft_output);
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

    snd_pcm_t *alsa_handle = init_alsa_capture(alsa_device);
    if (!alsa_handle) {
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

    while (app_running) {
        audio_update_levels(alsa_handle, lsr_rms, &lsr_balance);

        SDL_RenderClear(ren);
        SDL_RenderCopy(ren, bg, NULL, NULL);

        draw_spectrogram_bars(ren, spectro_band_db);
        draw_lr_bars(ren, lsr_rms, lsr_peak);
        draw_balance_bar(ren, lsr_balance);
		draw_phase_bar(ren, lsr_phase_rho_smoothed);
		draw_oscilloscopes(ren, osc_buf);

        SDL_RenderPresent(ren);
        SDL_Delay(1);
    }

    app_running = 0;
    pthread_cond_signal(&fft_cond);
    pthread_join(fft_thread, NULL);
    fftwf_destroy_plan(fft_plan);

    snd_pcm_close(alsa_handle);
    SDL_DestroyTexture(bg);
    SDL_DestroyRenderer(ren);
    SDL_DestroyWindow(win);
    IMG_Quit();
    SDL_Quit();

    return 0;
}
