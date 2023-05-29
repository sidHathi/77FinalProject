#include "ray.glsl"

#define metal 1
#define lambertian 2
#define dielectric 3

struct Material {
    vec3 albedo;
    int type;
    float matparam;
};