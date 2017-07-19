#include "light_source.h"

point_light::point_light(float x, float y, float z, rgb a, rgb d, rgb s, model3d *parent) {
    this->real_x = x;
    this->real_y = y;
    this->real_z = z;

    this->ambient = a;
    this->diffuse = d;
    this->specular = s;

    att_a = 1;
    att_b = 0.001;
    att_c = 0.0001;

    this->parent = parent;
}
