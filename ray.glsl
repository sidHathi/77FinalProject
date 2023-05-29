struct Ray {
    vec3 origin;
    vec3 direction;
};

vec3 ray_pointAt(Ray ray, float t) {
    return ray.origin + ray.direction * t;
}