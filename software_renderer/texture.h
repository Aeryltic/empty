#ifndef TEXTURE_H
#define TEXTURE_H

#include <QWidget>
#include "vertex.h"
#include "structures.h"

struct texture {
    uint32_t* buf;
    int width, height;

    texture(const texture& other);
    texture(); //blank 1x1
    virtual ~texture();

    void resize(int new_width, int new_height);
    bool load_from_file(QString filepath);

    uint32_t get(float u, float v) const { // w przedziale [0, 1]
        int x = (u) * (width - 1),
            y = (1 - v) * (height - 1); // dlaczego są odwrócone?
        return buf[y * width + x];
    }
};

struct material {
    texture tex;
    rgb emission;
    rgb ambient;
    rgb diffuse;
    rgb specular;
    float ns;
    float d;
    unsigned illum;

    material();
    virtual ~material();
};

#endif // TEXTURE_H
