#ifndef RENDERER3D_H
#define RENDERER3D_H

#include <stdint.h>

#include <armadillo>

#include "world.h"
#include "pixbuf.h"

class Renderer3D {
public:
    Renderer3D(int width, int height, uint32_t base_color);
    Renderer3D(const Renderer3D&) = delete;

    void render(World* world);
    const uint32_t* get_buffer() { return pixBuf.buf; }
    int get_width() { return width; }
    int get_height() { return height; }

private:
    int width, height;
    PixBuf pixBuf;
    float center_x, center_y;
    float aspect;
    float fov;
    float near, far;
    float x_mult, y_mult;
    World *world; // updated on render


    uint32_t calc_color(float x, float y, float z, float tu, float tv, const arma::vec& n, const material &mtl);
    void texture_bottom_tri(render_point v[], const material &mtl, const arma::vec &normal);
    void texture_top_tri(render_point v[], const material &mtl, const arma::vec &normal);
    void draw_textured_triangle(render_point v[], const material &mtl, const arma::vec &normal);

    bool is_inside(int x, int y) const {
        return (x>=0) && (x<width) && (y>=0) && (y<height);
    }
    void bresenham2d(int x1, int y1, int x2, int y2, uint32_t color);
    char cohen_sutherland_clipping(int &x0, int &y0, int &x1, int &y1);

    void draw_ellipse(int x, int y, double r1, double r2, double angle, uint32_t color, int s);
};

#endif // RENDERER3D_H
