#include <stdio.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include "incbin.h"
#include <time.h>
#include <math.h>

// ====== INCBIN : Chargement du PNG de fond ======
INCBIN(background_png, "backgroundVibeScope.png");

// ----------- Outils PNG SDL2 depuis mémoire -----------
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

// ====== Définition des barres du spectrogramme ======
#define NUM_COLS 31
const int bar_x[NUM_COLS] = {
    120,160,200,240,280,320,360,400,440,480,520,560,600,640,680,720,
    760,800,840,880,920,960,1000,1040,1080,1120,1160,1200,1240,1280,1320
};
const int bar_y = 192;
const int bar_w = 22;
const int bar_h = 800;

// ====== Données de spectre ======
float col_values[NUM_COLS] = {0.0f};
int selected_col = 0;

// ====== Barres L, S, R ======
#define LR_BAR_Y   192
#define LR_BAR_H   800
const int lr_bar_x[3] = { 1482, 1527, 1558 };
const int lr_bar_w[3] = { 36, 22, 36 };
float lr_values[3] = {0.6f, 0.4f, 0.8f}; // Exemple
const char *lr_labels[3] = {"L", "S", "R"};

// ====== Gestion des couleurs/degradés (identique SVG) ======
typedef struct { Uint8 r, g, b; } Color;

Color gradientStops[8] = {
    {0, 0, 0},      {0, 0, 128},    {0, 0, 255},    {0, 128, 0},
    {0, 255, 0},    {255, 255, 0},  {255, 128, 0},  {255, 0, 0}
};

// ----------- Dégradé vertical barres spectro -----------
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

// ====== Gestion du clip (cercles d'indicateur au dessus des barres) ======
int circle_cx(int i) { return bar_x[i] + bar_w / 2; }
const int circle_cy = 172;
const int circle_radius = 11;

// ----------- Couleur du cercle d’après niveau -----------
void getCircleColor(float v, Uint8 *r, Uint8 *g, Uint8 *b) {
    // Conversion linéaire : v=0 => -60dB, v=1 => +6dB
    float level_db = -60.0f + v * 66.0f;
    if (level_db < 0.0f) {
        *r = 0; *g = 255; *b = 0;     // Vert
    } else if (level_db < 3.0f) {
        *r = 255; *g = 140; *b = 0;   // Orange
    } else {
        *r = 255; *g = 0; *b = 0;     // Rouge
    }
}

// ----------- Cercle rempli (pour indicateurs clip) -----------
void drawFilledCircle(SDL_Renderer *renderer, int x0, int y0, int radius, Uint8 r, Uint8 g, Uint8 b, Uint8 a) {
    SDL_SetRenderDrawColor(renderer, r, g, b, a);
    for (int y = -radius; y <= radius; y++) {
        for (int x = -radius; x <= radius; x++) {
            if (x*x + y*y <= radius*radius)
                SDL_RenderDrawPoint(renderer, x0 + x, y0 + y);
        }
    }
}

// ====== Dégradés pour la balance L/R (rectangles horizontaux) ======

// ----------- Dégradé gauche (rouge-jaune-vert) -----------
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
// ----------- Dégradé droite (vert-jaune-rouge) -----------
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

// ====== Outils valeurs spectro ======
void set_column_value(int col, float value) {
    if (col < 0 || col >= NUM_COLS) return;
    if (value < 0.0f) value = 0.0f;
    if (value > 1.0f) value = 1.0f;
    col_values[col] = value;
    selected_col = col;
}

// ====== Affichage du spectrogramme ======
void draw_spectrogram_bars(SDL_Renderer *ren, const float *col_values) {
    for (int i = 0; i < NUM_COLS; ++i) {
        // Contour blanc
        SDL_Rect rc = {bar_x[i], bar_y, bar_w, bar_h};
        SDL_SetRenderDrawColor(ren, 255, 255, 255, 255);
        SDL_RenderDrawRect(ren, &rc);

        // Remplissage du dégradé
        int filled_h = (int)(bar_h * col_values[i]);
        for (int dy = 0; dy < filled_h; ++dy) {
            float val = (float)dy / (bar_h - 1);
            Color col = getGradientColor(val);
            SDL_SetRenderDrawColor(ren, col.r, col.g, col.b, 255);
            SDL_Rect seg = {bar_x[i] + 1, bar_y + bar_h - dy - 1, bar_w - 2, 1};
            SDL_RenderFillRect(ren, &seg);
        }

        // Cercles de niveau en haut
        Uint8 r, g, b;
        getCircleColor(col_values[i], &r, &g, &b);
        int cx = circle_cx(i);
        int cy = circle_cy;
        drawFilledCircle(ren, cx, cy, circle_radius, r, g, b, 255);
    }
}

// ====== Affichage des barres L, S, R (même style) ======
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

// ====== Affichage de la barre de balance L/R (-1 à +1) ======
void draw_balance_bar(SDL_Renderer *ren, float value)
{
    // Clamp [-1, 1]
    if (value < -1.0f) value = -1.0f;
    if (value >  1.0f) value =  1.0f;

    // Zones identiques SVG
    int left_x  = 1482, left_w  = 56;
    int right_x = 1538, right_w = 56;
    int bar_y   = 164,  bar_h   = 22;

    // Partie gauche
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
    // Partie droite
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

// ====== Boucle principale ======
int main(int argc, char *argv[]) {
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

    // ----------- INIT valeurs test -----------
    srand((unsigned)time(NULL));
    for (int i = 0; i < NUM_COLS; ++i)
        col_values[i] = (float)rand() / RAND_MAX;

    int running = 1;
    while (running) {
        SDL_RenderClear(ren);
        SDL_RenderCopy(ren, bg, NULL, NULL);

        // Spectro (31 barres)
        draw_spectrogram_bars(ren, col_values);

        // Barres L S R
        draw_lr_bars(ren, lr_values);

        // Balance L/R (exemple: -0.6f)
        float balance_value = -0.6f;
        draw_balance_bar(ren, balance_value);

        SDL_RenderPresent(ren);
        SDL_Delay(16);
    }

    SDL_DestroyTexture(bg);
    SDL_DestroyRenderer(ren);
    SDL_DestroyWindow(win);
    IMG_Quit(); SDL_Quit();
    return 0;
}
