#include <iostream>
#include <complex>
#include <chrono>
#include <vector>

#include <stdint.h>

#define SDL_MAIN_HANDLED
#include "SDL.h"

struct MandelLog {
    MandelLog() {

    }
};

std::complex<long double> f_c(std::complex<long double> z, std::complex<long double> c)
{
    return z * z + c;
}

const int x_size = 800;
const int y_size = 600;

struct complex_plane {
    long double re_range[2];
    long double im_range[2];
};

unsigned long max_iter = 100;

complex_plane get_plane(long double re_inf, long double re_sup, long double im_inf, long double im_sup)
{
    complex_plane plane;

    plane.re_range[0] = re_inf;
    plane.re_range[1] = re_sup;
    plane.im_range[0] = im_inf;
    plane.im_range[1] = im_sup;

    return plane;
}

int main(int argc, char **argv)
{
    std::cout << "Mandelbrot\n";

    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        std::cout << "SDL Error\n";
        return -1;
    }

    int mouse_pos_x;
    int mouse_pos_y;

    long double mouse_pos_plane_x;
    long double mouse_pos_plane_y;

    SDL_Window *win = SDL_CreateWindow("Mandelbrot",
                                       SDL_WINDOWPOS_UNDEFINED,
                                       SDL_WINDOWPOS_UNDEFINED,
                                       x_size,
                                       y_size,
                                       SDL_WINDOW_OPENGL | SDL_WINDOW_FULLSCREEN);
                                       // SDL_WINDOW_OPENGL);
    SDL_Renderer *renderer = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED);
    SDL_Texture *texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGB888, SDL_TEXTUREACCESS_TARGET, x_size, y_size);

    complex_plane plane;
    if (argc > 1) {
        plane = get_plane(atof(argv[1]), atof(argv[2]), atof(argv[3]), atof(argv[4]));
    } else {
        plane = get_plane(-2.0, 1.0, -1.0, 1.0);
    }

    SDL_Event evt;
    uint32_t *pixels = (uint32_t *)malloc(x_size * y_size * sizeof(uint32_t));
    int **mandel = (int **)malloc(x_size * sizeof(int *));
    if (!mandel)
    {
        std::cout << "Allocation error\n";

        return -1;
    }
    for (int idx = 0; idx < x_size; ++idx)
    {
        mandel[idx] = (int *)malloc(y_size * sizeof(int));
    }

    long double re_step;
    long double im_step;

    while (1)
    {
        SDL_UpdateTexture(texture, NULL, pixels, x_size * sizeof(uint32_t));
        SDL_WaitEvent(&evt);

        switch (evt.type)
        {
        case SDL_MOUSEBUTTONDOWN:
            if (evt.button.button == SDL_BUTTON_LEFT)
            {
                mouse_pos_x = evt.button.x;
                mouse_pos_y = evt.button.y;
                // mouse_pos_plane_x = (long double)(mouse_pos_x) * re_step;
                // mouse_pos_plane_y = (long double)(mouse_pos_y) * im_step;
                mouse_pos_x = mouse_pos_x - (x_size / 2);
                mouse_pos_y = mouse_pos_y - (y_size / 2);
                mouse_pos_plane_x = (long double)(mouse_pos_x) * re_step;
                mouse_pos_plane_y = (long double)(mouse_pos_y) * im_step;
                plane.re_range[0] += mouse_pos_plane_x;
                plane.re_range[1] += mouse_pos_plane_x;
                plane.im_range[0] += mouse_pos_plane_y;
                plane.im_range[1] += mouse_pos_plane_y;
                printf("mouse_pos_x: %d -- mouse_pos_y: %d\n", mouse_pos_x, mouse_pos_y);
            }
            break;
            /*
        case SDL_MOUSEMOTION:
            SDL_GetMouseState(&mouse_pos_x, &mouse_pos_y);
            // mouse_pos_plane_x = (long double)(mouse_pos_x) * re_step;
            // mouse_pos_plane_y = (long double)(mouse_pos_y) * im_step;
            mouse_pos_x = mouse_pos_x - (x_size / 2);
            mouse_pos_y = mouse_pos_y - (y_size / 2);
            mouse_pos_plane_x = (long double)(mouse_pos_x) * re_step;
            mouse_pos_plane_y = (long double)(mouse_pos_y) * im_step;
            break;
            */
        case SDL_KEYDOWN:
            // std::cout << plane.re_range[0] << "\n";
            // std::cout << plane.re_range[1] << "\n";

            if (evt.key.keysym.sym == SDLK_KP_PLUS)
            {
                plane.re_range[0] = plane.re_range[0] + abs(plane.re_range[0] - plane.re_range[1]) / 8.0;
                plane.re_range[1] = plane.re_range[1] - abs(plane.re_range[0] - plane.re_range[1]) / 8.0;
                plane.im_range[0] = plane.im_range[0] + abs(plane.im_range[0] - plane.im_range[1]) / 8.0;
                plane.im_range[1] = plane.im_range[1] - abs(plane.im_range[0] - plane.im_range[1]) / 8.0;
                // max_iter += 128;
            }
            else if (evt.key.keysym.sym == SDLK_KP_MINUS)
            {
                plane.re_range[0] = plane.re_range[0] - abs(plane.re_range[0] - plane.re_range[1]) / 8;
                plane.re_range[1] = plane.re_range[1] + abs(plane.re_range[0] - plane.re_range[1]) / 8;
                plane.im_range[0] = plane.im_range[0] - abs(plane.im_range[0] - plane.im_range[1]) / 8;
                plane.im_range[1] = plane.im_range[1] + abs(plane.im_range[0] - plane.im_range[1]) / 8;
                // max_iter -= 128;
            } else if (evt.key.keysym.sym == SDLK_s) {
                printf("re_step: %.30f\n", re_step);
                printf("im_step: %.30f\n", im_step);
            } else if (evt.key.keysym.sym == SDLK_j) {
                max_iter += 100;
                printf("max_iter: %u\n", max_iter);
            } else if (evt.key.keysym.sym == SDLK_k) {
                max_iter -= 100;
            } else if (evt.key.keysym.sym == SDLK_q) {
                SDL_Quit();
                exit(0);
            }
            break;
        }

        re_step = abs(plane.re_range[0] - plane.re_range[1]) / (long double)x_size;
        im_step = abs(plane.im_range[0] - plane.im_range[1]) / (long double)y_size;

        std::complex<long double> test_z(0, 0);
        long double magnitude;

        long double start_x = plane.re_range[0];
        long double start_y = plane.im_range[0];

        unsigned long int i = 0;
        auto s_mandel_time = std::chrono::steady_clock::now();
#pragma omp parallel for
        for (int x = 0; x < x_size; ++x)
        {
            for (int y = 0; y < y_size; ++y)
            {
                std::complex<long double> tmp(0.0, 0.0);
                std::complex<long double> test_c(start_x + (long double)x * re_step, start_y + (long double)y * im_step);
                for (i = 0; i < max_iter; ++i)
                {
                    tmp = f_c(tmp, test_c);
                    // abs is very slow but safer
                    // magnitude = abs(tmp);
                    magnitude = tmp.real() * tmp.real() + tmp.imag() * tmp.imag();
                    // std::cout << "magnitude: " << magnitude << std::endl;

                    if (magnitude >= 4)
                        break;
                }
                mandel[x][y] = i;
            }
        }

        // auto mandel_time = (std::chrono::steady_clock::now() - s_mandel_time);
        // std::cout << "mandel_time: " << std::chrono::duration_cast<std::chrono::milliseconds>(mandel_time).count() << "\n";

        Uint32 format = SDL_GetWindowPixelFormat(win);
        SDL_PixelFormat *mapping_format = SDL_AllocFormat(format);

        // memset(pixels, 255, x_size * y_size * sizeof(uint32_t));

        int iter;
        for (int j = 0; j < y_size; ++j)
        {
            for (int i = 0; i < x_size; ++i)
            {
                iter = mandel[i][j];
                if (iter >= max_iter)
                {
                    pixels[(j * x_size) + i] = 0;
                }
                else
                {
                    pixels[(j * x_size) + i] = SDL_MapRGB(mapping_format,
                                                          200 * cos(log(iter)),
                                                          250 * sin(log(iter)),
                                                          128 * cos(1.0 - log(iter)));
                    // pixels[(j * x_size) + i] = SDL_MapRGB(mapping_format,
                    //                                       255 * iter / max_iter,
                    //                                       100 * iter / max_iter,
                    //                                       50 * iter / max_iter);
                }
            }
        }

        SDL_RenderClear(renderer);
        SDL_RenderCopy(renderer, texture, NULL, NULL);
        SDL_RenderPresent(renderer);
    }

    // for (int j = 0; j < y_size; ++j)
    // {
    //     for (int i = 0; i < x_size; ++i)
    //     {
    //         if (mandel[i][j] >= max_iter) {
    //             std::cout << "x";
    //         } else {
    //             std::cout << " ";
    //         }
    //     }
    //     std::cout << std::endl;
    // }

    return 0;
}