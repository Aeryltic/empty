#include "camera.h"

camera::camera() {

}

void camera::update(int ms) {
    float vt = this->v * ms / 1000.0;
    float vx = vt * cos(ry)*cos(rx);
    float vz = vt * sin(ry)*cos(rx);
    float vy = vt * sin(rx);
    px += vx;
    py += vy;
    pz += vz;
}
