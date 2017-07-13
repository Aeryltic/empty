#ifndef CAMERA_H
#define CAMERA_H

#include <cmath>

struct camera
{
    float px, py, pz,
          rx, ry, rz,
          v;

    camera();
    void update(int ms);

    bool moving() {return abs(v)>0.001;}
};

#endif // CAMERA_H
