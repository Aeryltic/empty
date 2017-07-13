#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <cmath>

struct vec2d {
    float x, y;

//    vec2d operator+(const vec2d& other){
//        return {x + other.x, y + other.y};
//    }
//    vec2d operator-(const vec2d& other){
//        return {x - other.x, y - other.y};
//    }
};

struct vec3d {
    float x, y, z;

    vec3d operator+(const vec3d& other) const {
        return {x + other.x, y + other.y, z + other.z};
    }
    vec3d operator-(const vec3d& other) const {
        return {x - other.x, y - other.y, z - other.z};
    }

    float len() const {
        return std::sqrt(x * x + y * y + z * z);
    }
};

float dot(const vec3d& a, const vec3d& b);

typedef vec2d tri2d[3];
typedef vec3d triangle[3];

//struct bar_coor {
//    float u, v, w;

//    bar_coor(vec3d& p, triangle& t) {
//        vec3d v0 = t[1] - t[0], v1 = t[2] - t[0], v2 = p - t[0];
//        float d00 = dot(v0, v0);
//        float d01 = dot(v0, v1);
//        float d11 = dot(v1, v1);
//        float d20 = dot(v2, v0);
//        float d21 = dot(v2, v1);
//        float denom = d00 * d11 - d01 * d01;
//        v = (d11 * d20 - d01 * d21) / denom;
//        w = (d00 * d21 - d01 * d20) / denom;
//        u = 1.0f - v - w;
//    }
//};

struct rgb {
    float r, g, b; // w przedziale [0, 1]

    rgb() : r(0), g(0), b(0) {}
    rgb(float r, float g, float b) : r(r), g(g), b(b) {}

    rgb operator*(const rgb& other) const {
        return rgb(r * other.r, g * other.g, b * other.b);
    }
    rgb operator*(float v) const {
        return rgb(r * v, g * v, b * v);
    }
    rgb operator/(float v) const {
        return rgb(r / v, g / v, b / v);
    }
    rgb operator+(const rgb& other) const {
        return rgb(r + other.r, g + other.g, b + other.b);
    }
    void operator+=(const rgb& other) {
        r += other.r;
        g += other.g;
        b += other.b;
    }
};


#endif // STRUCTURES_H
