#ifndef SDL_CLASS_H
#define SDL_CLASS_H

#include <SDL2/SDL.h>
#include <iostream>
#include <algorithm>
#include "shapes.h"

/// @brief Classe pour gérer la SDL
class Sdl {
private:
    SDL_Window* window;
    SDL_Renderer* renderer;
    SDL_Texture* texture;
    int width, height;

public:
    /// @brief Init SDL et fenetre
    Sdl(const char* title, int w, int h) : width(w), height(h) {
        if (SDL_Init(SDL_INIT_VIDEO) < 0) exit(1); // Pas d'exception, on quitte direct

        window = SDL_CreateWindow(title, 0, 0, w, h, SDL_WINDOW_SHOWN);
        if (!window) exit(1);

        // Hardware accel
        renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
        if (!renderer) exit(1);
        
        SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_NONE);

        // Texture pour CUDA
        texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA32, 
                                    SDL_TEXTUREACCESS_STREAMING, w, h);
    }

    /// @brief Destructeur
    ~Sdl() {
        SDL_DestroyTexture(texture);
        SDL_DestroyRenderer(renderer);
        SDL_DestroyWindow(window);
        SDL_Quit();
    }

    /// @brief Efface l'ecran (Hack linux inclus)
    void clear() {
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderClear(renderer);

        // Fix pour le ghosting: on dessine un rect noir par dessus
        SDL_Rect r = {0, 0, width, height};
        SDL_RenderFillRect(renderer, &r);
    }

    /// @brief Affiche le rendu
    void present() {
        SDL_RenderPresent(renderer);
    }

    /// @brief Change couleur
    void set_color(int r, int g, int b) {
        SDL_SetRenderDrawColor(renderer, r, g, b, 255);
    }

    // RASTERISATION

    /// @brief Interpolation X
    int interp(int y, Point2d p1, Point2d p2) {
        if (p1.y == p2.y) return p1.x;
        float t = (float)(y - p1.y) / (float)(p2.y - p1.y);
        return p1.x + (int)(t * (p2.x - p1.x));
    }

    /// @brief Remplissage triangle scanline
    void draw_filled_triangle(Point2d p1, Point2d p2, Point2d p3) {
        // Tri
        if (p1.y > p2.y) std::swap(p1, p2);
        if (p1.y > p3.y) std::swap(p1, p3);
        if (p2.y > p3.y) std::swap(p2, p3);

        // Haut
        for (int y = p1.y; y <= p2.y; y++) {
            if (y < 0 || y >= height) continue;
            int xa = interp(y, p1, p3);
            int xb = interp(y, p1, p2);
            if (xa > xb) std::swap(xa, xb);
            SDL_RenderDrawLine(renderer, xa, y, xb, y);
        }
        // Bas
        for (int y = p2.y; y <= p3.y; y++) {
            if (y < 0 || y >= height) continue;
            int xa = interp(y, p1, p3);
            int xb = interp(y, p2, p3);
            if (xa > xb) std::swap(xa, xb);
            SDL_RenderDrawLine(renderer, xa, y, xb, y);
        }
    }

    // CUDA

    /// @brief Affiche buffer pixels
    void draw_buffer(void* pixels, int pitch) {
        SDL_UpdateTexture(texture, NULL, pixels, pitch);
        SDL_RenderCopy(renderer, texture, NULL, NULL);
    }
    
    /// @brief Events
    bool poll_quit() {
        SDL_Event e;
        while(SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) return true;
            if (e.type == SDL_KEYDOWN && e.key.keysym.sym == SDLK_ESCAPE) return true;
        }
        return false;
    }
};

#endif