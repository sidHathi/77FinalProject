#include "fluid-utils.glsl"

struct FluidCube {
    float scale;
    vec3 center;
    int size;
    float dt;
    float diff;
    float visc;

    float pdensity[MAX_ARR_SIZE];
    float density[MAX_ARR_SIZE];

    float Vx[MAX_ARR_SIZE];
    float Vy[MAX_ARR_SIZE];
    float Vz[MAX_ARR_SIZE];

    float Vx0[MAX_ARR_SIZE];
    float Vy0[MAX_ARR_SIZE];
    float Vz0[MAX_ARR_SIZE];
};

FluidCube fc_create(float scale, vec3 center, int size, float diffusion, float viscosity, float dt) {
    float pdensity[MAX_ARR_SIZE];
    float density[MAX_ARR_SIZE];
    float Vx[MAX_ARR_SIZE];
    float Vy[MAX_ARR_SIZE];
    float Vz[MAX_ARR_SIZE];

    float Vx0[MAX_ARR_SIZE];
    float Vy0[MAX_ARR_SIZE];
    float Vz0[MAX_ARR_SIZE];

    FluidCube c = FluidCube(
        scale,
        center,
        size,
        dt,
        diffusion,
        viscosity,
        pdensity,
        density,
        Vx, 
        Vy,
        Vz,
        Vx0,
        Vy0,
        Vz0
    );
    return c;
}

void fc_add_density(inout FluidCube cube, int x, int y, int z, float amt) {
    int N = cube.size;
    int index = IX(x, y, z, N);
    cube.density[index] += amt;
}

void fc_add_velocity(inout FluidCube cube, int x, int y, int z, float amtX, float amtY, float amtZ) {
    int N = cube.size;
    int index = IX(x, y, z, N);

    cube.Vx[index] += amtX;
    cube.Vy[index] += amtY;
    cube.Vz[index] += amtZ;
}

void fc_step(inout FluidCube cube) {
    int N          = cube.size;
    float visc     = cube.visc;
    float diff     = cube.diff;
    float dt       = cube.dt;

    float Vx[MAX_ARR_SIZE] = cube.Vx;
    float Vy[MAX_ARR_SIZE] = cube.Vy;
    float Vz[MAX_ARR_SIZE] = cube.Vz;

    float Vx0[MAX_ARR_SIZE] = cube.Vx0;
    float Vy0[MAX_ARR_SIZE] = cube.Vy0;
    float Vz0[MAX_ARR_SIZE] = cube.Vz0;

    float s[MAX_ARR_SIZE] = cube.pdensity;;
    float density[MAX_ARR_SIZE] = cube.density;
    
    diffuse(1, Vx0, Vx, visc, dt, 4, N);
    diffuse(2, Vy0, Vy, visc, dt, 4, N);
    diffuse(3, Vz0, Vz, visc, dt, 4, N);
    
    project(Vx0, Vy0, Vz0, Vx, Vy, 4, N);
    
    advect(1, Vx, Vx0, Vx0, Vy0, Vz0, dt, N);
    advect(2, Vy, Vy0, Vx0, Vy0, Vz0, dt, N);
    advect(3, Vz, Vz0, Vx0, Vy0, Vz0, dt, N);
    
    project(Vx, Vy, Vz, Vx0, Vy0, 4, N);
    
    diffuse(0, s, density, diff, dt, 4, N);
    advect(0, density, s, Vx, Vy, Vz, dt, N);
}

float fc_point_density(in FluidCube cube, in vec3 point) {
    int N = cube.size;
    // naive -> improve to binary search later
    vec3 bottomLeft = cube.center - vec3(float(N)*cube.scale / 2.);
    for (int i = 0; i < N; i ++) {
        for (int j = 0; j < N; j ++) {
            for (int k = 0; k < N; k ++) {
                float xLB = bottomLeft.x + (cube.scale * float(i)) + cube.scale * 0.;
                float xGB = bottomLeft.x + (cube.scale * float(i)) + cube.scale * 1.0;
                float yLB = bottomLeft.y + (cube.scale * float(j)) + cube.scale * 0.;
                float yGB = bottomLeft.y + (cube.scale * float(j)) + cube.scale * 1.0;
                float zLB = bottomLeft.z + (cube.scale * float(k)) + cube.scale * 0.;
                float zGB = bottomLeft.z + (cube.scale * float(k)) + cube.scale * 1.0;
                if (point.x > xLB && point.x < xGB && point.y > yLB && point.y < yGB && point.z > zLB && point.z < zGB) {
                    return cube.density[IX(i, j, k, N)];
                }
            }
        }
    }
    return 1.;
}

float sdBox(vec3 p, vec3 half_bounds) {
    float d = 0.;
    if (abs(p.x) < half_bounds.x && abs(p.y) < half_bounds.y && abs(p.z) < half_bounds.z) {
        return max(max(abs(p.x) - half_bounds.x, abs(p.y) - half_bounds.y), abs(p.z) - half_bounds.z);
    } else if (abs(p.x) < half_bounds.x) {
        if (abs(p.y) < half_bounds.y) {
            d = abs(p.z) - half_bounds.z;
        } else {
            if (abs(p.z) < half_bounds.z) {
                d = abs(p.y) - half_bounds.y;
            } else {
                d = length(vec2(abs(p.y) - half_bounds.y, abs(p.z) - half_bounds.z));
            }
        }
    } else if (abs(p.y) <= half_bounds.y) {
        if (abs(p.z) <= half_bounds.z) {
            d = abs(p.x) - half_bounds.x;
        } else {
            d = length(vec2(abs(p.x) - half_bounds.x, abs(p.z) - half_bounds.z));
        }
    } else {
        if (abs(p.z) <= half_bounds.z) {
            d = length(vec2(abs(p.x) - half_bounds.x, abs(p.y) - half_bounds.y));
        } else {
            d = length(vec3(abs(p.x) - half_bounds.x, abs(p.y) - half_bounds.y, abs(p.z) - half_bounds.z));
        }
    }
    return d;
}