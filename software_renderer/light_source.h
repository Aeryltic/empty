#ifndef LIGHT_SOURCE_H
#define LIGHT_SOURCE_H

#include "structures.h"
#include "model3d.h"

struct light_source {
    float real_x, real_y, real_z; // rzeczywiste pozycje (lub wzgl parenta)
    float x, y, z; // wzgl kamery

    rgb ambient;
    rgb diffuse;
    rgb specular;

//    float radius;
//    float a, b, c;
//    float intensity;
//    float d_max;

    float att_a, att_b, att_c;

    model3d* parent;

    light_source(float x, float y, float z, rgb a, rgb d, rgb s, model3d* parent = nullptr);//, float radius);
};

struct inf_light_source {
    float ray_x, ray_y, ray_z;
    float x, y, z;

    rgb ambient;
    rgb diffuse;
    rgb specular;

    float intensity;
};

#endif // LIGHT_SOURCE_H
