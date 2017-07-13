#ifndef PIXBUF_H
#define PIXBUF_H

#include <QWidget>
#include <armadillo>

#include "model3d.h"
#include "camera.h"

#include "structures.h"

struct PixBuf {

    int width, height;
    uint32_t* buf;
    float* z_buf;
    uint32_t baseColor;


    PixBuf(int width, int height, uint32_t color);
    PixBuf(const PixBuf&) = delete;
    ~PixBuf();

    bool resize(int new_width, int new_height);

    void clear();
    void clear_zbuf(float far);

//    bool is_inside(int x, int y) const {
//        return (x>=0) && (x<width) && (y>=0) && (y<height);
//    }
//    void set_pixel(int x, int y, uint32_t color) {
//        if(is_inside(x, y))
//            buf[y*width+x] = color;
//    }

};

#endif // PIXBUF_H
