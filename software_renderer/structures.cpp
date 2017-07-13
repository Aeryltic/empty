#include "structures.h"

float dot(const vec3d& a, const vec3d& b){
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
