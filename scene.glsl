#include "sphere.glsl"

struct Scene {
    Sphere sphere;
};

Scene sc_init() {
    Sphere sphere = Sphere(vec3(0., 0., -1.), 0.5);
    return Scene(
        sphere
    );
}

vec3 sc_ray_color(in Scene s, in Ray r) {
    HitRecord rec;
    if (sphere_hit(s.sphere, r, 0., 1./0., rec)) {
        return 0.5 * (rec.normal + vec3(1., 1., 1.));
    }
    vec3 unit_dir = normalize(r.direction);
    float t = 0.5*(unit_dir.y + 1.0);
    return (1.0 - t) * vec3(1., 1., 1.) + t*vec3(0.5, 0.7, 1.0);
}

vec3 sc_col(in Scene s) {
    /*
        Initial implementation:
        Create a ray from the origin through the uv coords of the screen
        return the color provided by the ray_color function
    */

    float aspect_ratio = iResolution.x / iResolution.y;
    float viewport_height = 2.;
    float viewport_width = aspect_ratio*viewport_height;
    float focal_length = 1.;

    vec3 origin = vec3(0., 0., 0.);
    vec3 horizontal = vec3(viewport_width, 0., 0.);
    vec3 vertical = vec3(0., viewport_height, 0.);
    vec3 lower_left_corner = origin - (horizontal/2.) - (vertical/2.) - vec3(0., 0., focal_length);

    vec2 uv = gl_FragCoord.xy / iResolution.xy;
    float u = uv.x;
    float v = uv.y;
    Ray ray = Ray(origin, lower_left_corner + u*horizontal + v*vertical - origin);
    return sc_ray_color(s, ray);
}