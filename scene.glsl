#include "sphere.glsl"
#include "fluid.glsl"

struct Scene {
    Sphere sphere;
    FluidCube cube;
};

Scene sc_init() {
    Sphere sphere = Sphere(vec3(0., 0., -1.), 0.5);
    FluidCube cube = fc_create(0.1, vec3(0., 0., -0.5), 64, 0.2, 0., 0.000001);
    fc_add_density(cube, 0, 0, 0, 0.2);
    fc_add_velocity(cube, 0, 0, 0, 10., 10., 10.);
    return Scene(
        sphere,
        cube
    );
}

bool rayFluidInteract(in FluidCube cube, in Ray r, float t_min, float t_max, inout HitRecord rec) {
    // Strategy -> treat the cube like an actual cube
    // find closest hit point and normal
    // reflect based on the density of the cube at the hit point coord
    // requires being able to look up the density at a point
    int iters = 1;
    int max_iter = 1000;
    vec3 loc = r.origin;
    bool hit = false;
    float step_size = 1.;
    while (!hit && iters < max_iter) {
        float half_length = cube.scale * float(cube.size) / 2.;
        vec3 half_bounds = vec3(
            cube.center.x + half_length,
            cube.center.y + half_length,
            cube.center.z + half_length
        );
        step_size = sdBox(loc, half_bounds);
        loc = loc + r.direction * step_size;
        if (step_size < EPSILON) {
            hit = true;
            float density = fc_point_density(cube, loc);
            rec.t = length((loc  - r.origin))/length(r.direction);
            rec.p = loc;
            vec3 outward_normal = normalize((rec.p - cube.center));
            hr_set_face_normal(rec, r, outward_normal);
            return true;
        }
        iters += 1;
    }
    return false;
}

vec3 sc_ray_color(in Scene s, in Ray r) {
    HitRecord rec;
    if (sphere_hit(s.sphere, r, 0., 1./0., rec) || rayFluidInteract(s.cube, r, 0., 1./0., rec)) {
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
    fc_step(s.cube);
    return sc_ray_color(s, ray);
}