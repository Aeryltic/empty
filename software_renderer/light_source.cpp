#include "light_source.h"

light_source::light_source(float x, float y, float z, rgb a, rgb d, rgb s) {//, float radius) {
    this->real_x = x;
    this->real_y = y;
    this->real_z = z;

    this->ambient = a;
    this->diffuse = d;
    this->specular = s;

//    this->radius = radius < 1 ? 1 : radius;
}

/*
vec3 DirectIllumination(vec3 P, vec3 N, vec3 lightCentre, float lightRadius, vec3 lightColour, float cutoff)
{
    // calculate normalized light vector and distance to sphere light surface
    float r = lightRadius;
    vec3 L = lightCentre - P;
    float distance = length(L);
    float d = max(distance - r, 0);
    L /= distance;

    // calculate basic attenuation
    float denom = d/r + 1;
    float attenuation = 1 / (denom*denom);

    // scale and bias attenuation such that:
    //   attenuation == 0 at extent of max influence
    //   attenuation == 1 when d == 0
    attenuation = (attenuation - cutoff) / (1 - cutoff);
    attenuation = max(attenuation, 0);

    float dot = max(dot(L, N), 0);
    return lightColour * dot * attenuation;
}
*/
