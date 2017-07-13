#include "texture.h"

#include <QImage>

texture::texture() {
    width = height = 3;
    buf = new uint32_t[width * height];
    for(int i = 0; i < width * height; i++) {
        buf[i] = 0xffffffff;
    }
}


texture::texture(const texture &other) {
    width = other.width;
    height = other.height;
    buf = new uint32_t[width * height];
    for(int i = 0; i < width * height; i++) {
        buf[i] = other.buf[i];
    }
}

bool texture::load_from_file(QString filepath) {
    QImage img;
    if(!img.load(filepath)) {
        return false;
    } else {
        qDebug("texture loaded: %s", filepath.toStdString().c_str());

        const uchar* bits = img.constBits();
        resize(img.width(), img.height());

        for(int y = 0; y < height; y++) {
            for(int x = 0; x < width; x++) {
                int i = (y * width + x) * 4;
                uint32_t a = bits[i+3],
                         r = bits[i+2],
                         g = bits[i+1],
                         b = bits[i];
                uint32_t color = a << 24 | r << 16 | g << 8 | b;
                buf[y * width + x] = color;
            }
        }
        qDebug("size: %u, %u", width, height);
    }
    return true;
}

void texture::resize(int new_width, int new_height) {
    if(width != new_width || height != new_height) {
        delete [] buf;
        buf = new uint32_t[new_width * new_height];
        width = new_width;
        height = new_height;
    }
}

texture::~texture() {
    if(buf != nullptr) {
        delete [] buf;
        buf = nullptr;
    }
}

//---------MATERIAL--------------------------

material::material() : ambient(.2, .2, .2) {
    d = 1;
    illum = 0;
}

material::~material() {
}
