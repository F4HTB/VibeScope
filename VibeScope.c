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

// ====== Defines ======
#define NUM_COLS 31
#define LR_BAR_Y   192
#define LR_BAR_H   800
#define AUDIO_BUFFER_SIZE 1024
#define FFT_SIZE 1024

// ====== Variables globales ======
volatile sig_atomic_t running = 1;

// --- FFT synchronisation/partage ---
short fft_samples[FFT_SIZE * 2]; // stéréo (L,R,...)
pthread_mutex_t fft_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t fft_cond = PTHREAD_COND_INITIALIZER;
volatile int new_block_ready = 0;

// --- Spectre et affichage ---
float col_values[NUM_COLS] = {0.0f};
int selected_col = 0;

const int bar_x[NUM_COLS] = {
    120,160,200,240,280,320,360,400,440,480,520,560,600,640,680,720,
    760,800,840,880,920,960,1000,1040,1080,1120,1160,1200,1240,1280,1320
};
const int bar_y = 192;
const int bar_w = 22;
const int bar_h = 800;

const int lr_bar_x[3] = { 1482, 1527, 1558 };
const int lr_bar_w[3] = { 36, 22, 36 };
float lr_values[3] = {0.0f};
float balance_value = 0.0f;
const char *lr_labels[3] = {"L", "S", "R"};

// --- Couleurs gradients ---
typedef struct { Uint8 r, g, b; } Color;
Color gradientStops[8] = {
    {0, 0, 0},      {0, 0, 128},    {0, 0, 255},    {0, 128, 0},
    {0, 255, 0},    {255, 255, 0},  {255, 128, 0},  {255, 0, 0}
};
Color gradientLUT[800]; //const int bar_h = 800;
int gradientLUT_initialized = 0;
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
void set_column_value(int col, float value);
void draw_spectrogram_bars(SDL_Renderer *ren, const float *col_values);
void draw_lr_bars(SDL_Renderer *ren, const float *lr_values);
void draw_balance_bar(SDL_Renderer *ren, float value);
snd_pcm_t *init_alsa_capture(const char *device);
void update_audio_levels(snd_pcm_t *handle, float *lr_values, float *balance_value);
void* fft_thread_func(void *arg);

// ====== Fonctions ======

// -- Signal handler --
void handle_sigint(int sig) {
    running = 0;
}

// -- INCBIN : Chargement du PNG de fond --
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

// -- Dégradé vertical barres spectro --
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

// -- Coordonnées cercles indicateurs --
int circle_cx(int i) { return bar_x[i] + bar_w / 2; }

// -- Couleur du cercle d’après niveau --
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

// -- Cercle rempli --
void drawFilledCircle(SDL_Renderer *renderer, int x0, int y0, int radius, Uint8 r, Uint8 g, Uint8 b, Uint8 a) {
    SDL_SetRenderDrawColor(renderer, r, g, b, a);
    for (int y = -radius; y <= radius; y++) {
        for (int x = -radius; x <= radius; x++) {
            if (x*x + y*y <= radius*radius)
                SDL_RenderDrawPoint(renderer, x0 + x, y0 + y);
        }
    }
}

// -- Dégradés balance L/R --
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

// -- Outils valeurs spectro --
void set_column_value(int col, float value) {
    if (col < 0 || col >= NUM_COLS) return;
    if (value < 0.0f) value = 0.0f;
    if (value > 1.0f) value = 1.0f;
    col_values[col] = value;
    selected_col = col;
}

// -- Affichage du spectrogramme --
void draw_spectrogram_bars(SDL_Renderer *ren, const float *col_values) {
    for (int i = 0; i < NUM_COLS; ++i) {
        SDL_Rect rc = {bar_x[i], bar_y, bar_w, bar_h};
        SDL_SetRenderDrawColor(ren, 255, 255, 255, 255);
        SDL_RenderDrawRect(ren, &rc);

        int filled_h = (int)(bar_h * col_values[i]);
		if (filled_h < 0) filled_h = 0;
		if (filled_h > bar_h) filled_h = bar_h;

        // --- Ici tu utilises gradientLUT au lieu de recalculer chaque pixel ---
        for (int dy = 0; dy < filled_h; ++dy) {
            // float val = (float)dy / (bar_h - 1);    // Plus besoin
            // Color col = getGradientColor(val);      // Plus besoin
            Color col = gradientLUT[dy];
            SDL_SetRenderDrawColor(ren, col.r, col.g, col.b, 255);
            SDL_Rect seg = {bar_x[i] + 1, bar_y + bar_h - dy - 1, bar_w - 2, 1};
            SDL_RenderFillRect(ren, &seg);
        }

        // Cercles de niveau en haut (inchangé)
        Uint8 r, g, b;
        getCircleColor(col_values[i], &r, &g, &b);
        int cx = circle_cx(i);
        int cy = circle_cy;
        drawFilledCircle(ren, cx, cy, circle_radius, r, g, b, 255);
    }
}


// -- Affichage des barres L, S, R --
void draw_lr_bars(SDL_Renderer *ren, const float *lr_values) {
    for (int i = 0; i < 3; ++i) {
        SDL_Rect rc = {lr_bar_x[i], LR_BAR_Y, lr_bar_w[i], LR_BAR_H};
        SDL_SetRenderDrawColor(ren, 255, 255, 255, 255);
        SDL_RenderDrawRect(ren, &rc);

        int filled_h = (int)(LR_BAR_H * lr_values[i]);
        for (int dy = 0; dy < filled_h; ++dy) {
            float val = (float)dy / (LR_BAR_H - 1);
            Color col = getGradientColor(val);
            SDL_SetRenderDrawColor(ren, col.r, col.g, col.b, 255);
            SDL_Rect seg = {lr_bar_x[i] + 1, LR_BAR_Y + LR_BAR_H - dy - 1, lr_bar_w[i] - 2, 1};
            SDL_RenderFillRect(ren, &seg);
        }
    }
}

// -- Affichage de la barre de balance L/R --
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

// -- ALSA audio capture --
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
    unsigned int rate = 48000;
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

void update_audio_levels(snd_pcm_t *handle, float *lr_values, float *balance_value) {
    short buffer[AUDIO_BUFFER_SIZE * 2];
    int nframes;

    while (1) {
        nframes = snd_pcm_readi(handle, buffer, AUDIO_BUFFER_SIZE);

        if (nframes == -EAGAIN) {
            SDL_Delay(1);
            continue;
        }
        if (nframes == -EPIPE) {
            snd_pcm_prepare(handle);
            continue;
        }
        if (nframes < 0) {
            fprintf(stderr, "Erreur de lecture ALSA: %s\n", snd_strerror(nframes));
            lr_values[0] = lr_values[1] = lr_values[2] = 0.0f;
            *balance_value = 0.0f;
            return;
        }
        break;
    }
    double sum_l = 0, sum_r = 0;
    for (int i = 0; i < nframes; ++i) {
        float sample_l = buffer[2*i]   / 32768.0f;
        float sample_r = buffer[2*i+1] / 32768.0f;
        sum_l += sample_l * sample_l;
        sum_r += sample_r * sample_r;
    }
    float rms_l = sqrt(sum_l / nframes);
    float rms_r = sqrt(sum_r / nframes);

    lr_values[0] = rms_l;
    lr_values[2] = rms_r;
    lr_values[1] = (rms_l + rms_r) * 0.5f;
    float total = rms_l + rms_r;
    if (total > 1e-6)
        *balance_value = (rms_r - rms_l) / total;
    else
        *balance_value = 0.0f;

    pthread_mutex_lock(&fft_mutex);
    memcpy(fft_samples, buffer, sizeof(short) * 2 * FFT_SIZE);
    new_block_ready = 1;
    pthread_cond_signal(&fft_cond);
    pthread_mutex_unlock(&fft_mutex);
}

// -- Gestion de la fft --
void* fft_thread_func(void *arg) {
    fftwf_complex out[FFT_SIZE];
    float in[FFT_SIZE];
    short block_stereo[FFT_SIZE * 2];

    fftwf_plan plan = fftwf_plan_dft_r2c_1d(FFT_SIZE, in, out, FFTW_ESTIMATE); //FFTW_MEASURE

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

        for (int i = 0; i < FFT_SIZE; i++) {
            float L = block_stereo[2*i]   / 32768.0f;
            float R = block_stereo[2*i+1] / 32768.0f;
            in[i] = 0.5f * (L + R);
        }

        for (int i = 0; i < FFT_SIZE; i++)
            in[i] *= 0.5f * (1.0f - cosf(2.0f * M_PI * i / (FFT_SIZE - 1)));

        fftwf_execute(plan);

        float mag[FFT_SIZE / 2];
        for (int i = 0; i < FFT_SIZE / 2; i++)
            mag[i] = hypotf(out[i][0], out[i][1]);

		float band_freqs[NUM_COLS] = {
			25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250,
			315, 400, 500, 630, 800, 1000, 1250, 1600, 2000,
			2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500, 16000, 20000
		};
		float sample_rate = 48000.0f;
		float freq_per_bin = sample_rate / FFT_SIZE;

		for (int c = 0; c < NUM_COLS; c++) {
			float f_low, f_high;
			if (c == 0)
				f_low = 0;
			else
				f_low = (band_freqs[c-1] + band_freqs[c]) / 2.0f;
			if (c == NUM_COLS-1)
				f_high = sample_rate / 2.0f;
			else
				f_high = (band_freqs[c] + band_freqs[c+1]) / 2.0f;

			int bin_start = (int)(f_low / freq_per_bin);
			int bin_end   = (int)(f_high / freq_per_bin);
			if (bin_end > FFT_SIZE / 2) bin_end = FFT_SIZE / 2;
			float sum = 0.0f;
			int n_bins = bin_end - bin_start;
			for (int i = bin_start; i < bin_end; i++) sum += mag[i];
			col_values[c] = (n_bins > 0) ? (sum / n_bins) : 0.0f;
		}

        float maxval = 0.0f;
        for (int c = 0; c < NUM_COLS; c++)
            if (col_values[c] > maxval) maxval = col_values[c];
        if (maxval > 1e-6f)
            for (int c = 0; c < NUM_COLS; c++)col_values[c] /= maxval;

        usleep(1000);
    }
    fftwf_destroy_plan(plan);
    return NULL;
}
void init_gradientLUT() {
    for (int dy = 0; dy < bar_h; ++dy) {
        float val = (float)dy / (bar_h - 1);
        gradientLUT[dy] = getGradientColor(val);
    }
    gradientLUT_initialized = 1;
}
// ====== Boucle principale ======
int main(int argc, char *argv[]) {

    const char *alsa_device = "default";
    for (int i = 1; i < argc - 1; ++i) {
        if (strcmp(argv[i], "-D") == 0) {
            alsa_device = argv[i+1];
            break;
        }
    }
	
	init_gradientLUT();

    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        fprintf(stderr, "SDL_Init error: %s\n", SDL_GetError());
        return 1;
    }

    if (!(IMG_Init(IMG_INIT_PNG) & IMG_INIT_PNG)) {
        fprintf(stderr, "IMG_Init error: %s\n", IMG_GetError());
        SDL_Quit(); return 1;
    }

    SDL_Window *win = SDL_CreateWindow("Audio Analyser",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
        1920, 1080, SDL_WINDOW_SHOWN);
    SDL_SetWindowFullscreen(win, SDL_WINDOW_FULLSCREEN_DESKTOP);
    SDL_Renderer *ren = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED);

    SDL_Texture *bg = loadPngFromMemory(ren);
    if (!bg) {
        fprintf(stderr, "Load memory PNG error: %s\n", SDL_GetError());
        SDL_DestroyRenderer(ren); SDL_DestroyWindow(win); IMG_Quit(); SDL_Quit(); return 1;
    }

    // ----------- INIT carte son -----------
    snd_pcm_t *pcm_handle = init_alsa_capture(alsa_device);
    if (!pcm_handle) {
        SDL_DestroyRenderer(ren);
        SDL_DestroyWindow(win);
        IMG_Quit(); SDL_Quit();
        return 1;
    }

    pthread_t fft_thread;
    pthread_create(&fft_thread, NULL, fft_thread_func, NULL);

    signal(SIGINT, handle_sigint);
    while (running) {

        update_audio_levels(pcm_handle, lr_values, &balance_value);

        SDL_RenderClear(ren);
        SDL_RenderCopy(ren, bg, NULL, NULL);

        // Spectro (31 barres)
        draw_spectrogram_bars(ren, col_values);

        // Barres L S R
        draw_lr_bars(ren, lr_values);

        // Balance L/R
        draw_balance_bar(ren, balance_value);

        SDL_RenderPresent(ren);
        SDL_Delay(1);
    }

    running = 0;
    pthread_cond_signal(&fft_cond); // Réveille le thread FFT si endormi
    pthread_join(fft_thread, NULL);
    snd_pcm_close(pcm_handle);
    SDL_DestroyTexture(bg);
    SDL_DestroyRenderer(ren);
    SDL_DestroyWindow(win);
    IMG_Quit(); SDL_Quit();
    return 0;
}
