#ifndef LIGHT_SOURCE_H
#define LIGHT_SOURCE_H

#include "structures.h"

struct light_source {
    float real_x, real_y, real_z; // rzeczywiste pozycje
    float x, y, z; // przeliczone wzgl kamery

    rgb ambient;
    rgb diffuse;
    rgb specular;

//    float radius;
//    float a, b, c;
//    float intensity;
//    float d_max;

    light_source(float x, float y, float z, rgb a, rgb d, rgb s);//, float radius);
};

#endif // LIGHT_SOURCE_H
