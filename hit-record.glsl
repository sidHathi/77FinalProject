#include "material.glsl"

struct HitRecord {
    vec3 p;
    vec3 normal;
    float t;
    bool front_face;
    float dm;
    Material mat;
};

void hr_set_face_normal(inout HitRecord hr, Ray r, vec3 outward_normal) {
    hr.front_face = dot(r.direction, outward_normal) < 0.;
    hr.normal = hr.front_face ? outward_normal : -outward_normal;
}