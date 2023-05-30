#include "hit-record.glsl"

struct Sphere {
    vec3 center;
    float radius;
};

bool sphere_hit(in Sphere s, Ray r, float t_min, float t_max, inout HitRecord rec) {
    vec3 oc = r.origin - s.center;
    float a = pow(length(r.direction), 2.);
    float hb = dot(oc, r.direction);
    float c = pow(length(oc), 2.) - s.radius*s.radius;

    float discriminant = hb*hb - a*c;
    if (discriminant < 0.0001) {
        return false;
    } 
    float sqd = sqrt(discriminant);

    float root = (-hb - sqd) / a;
    if (root <= t_min || t_max <= root) {
        root = (-hb + sqd) / a;
        if (root <= t_min || t_max <= root) {
            return false;
        }
    }

    rec.t = root;
    rec.p = ray_pointAt(r, rec.t);
    vec3 outward_normal = normalize((rec.p - s.center) / s.radius);
    hr_set_face_normal(rec, r, outward_normal);
    rec.dm = 1.;
    // rec.mat = s.mat;

    return true;
}