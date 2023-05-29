#define EPSILON 1e-4
#define PI 3.1415926535897932384626433832795

// This is a constraint of GLSL and means N can be no larger than 64
const int MAX_ARR_SIZE = 262144;

int IX(int x, int y, int z, int N) {
    return x + y*N + z*N*N;
}

void set_bnd(in int b, inout float x[MAX_ARR_SIZE], int N) {
    for(int j = 1; j < N - 1; j++) {
        for(int i = 1; i < N - 1; i++) {
            x[IX(i, j, 0, N)] = b == 3 ? -x[IX(i, j, 1, N)] : x[IX(i, j, 1, N)];
            x[IX(i, j, N-1, N)] = b == 3 ? -x[IX(i, j, N-2, N)] : x[IX(i, j, N-2, N)];
        }
    }
    for(int k = 1; k < N - 1; k++) {
        for(int i = 1; i < N - 1; i++) {
            x[IX(i, 0  , k, N)] = b == 2 ? -x[IX(i, 1  , k, N)] : x[IX(i, 1  , k, N)];
            x[IX(i, N-1, k, N)] = b == 2 ? -x[IX(i, N-2, k, N)] : x[IX(i, N-2, k, N)];
        }
    }
    for(int k = 1; k < N - 1; k++) {
        for(int j = 1; j < N - 1; j++) {
            x[IX(0  , j, k, N)] = b == 1 ? -x[IX(1  , j, k, N)] : x[IX(1  , j, k, N)];
            x[IX(N-1, j, k, N)] = b == 1 ? -x[IX(N-2, j, k, N)] : x[IX(N-2, j, k, N)];
        }
    }
    
    x[IX(0, 0, 0, N)]       = 0.33 * (x[IX(1, 0, 0, N)]
                                  + x[IX(0, 1, 0, N)]
                                  + x[IX(0, 0, 1, N)]);
    x[IX(0, N-1, 0, N)]     = 0.33 * (x[IX(1, N-1, 0, N)]
                                  + x[IX(0, N-2, 0, N)]
                                  + x[IX(0, N-1, 1, N)]);
    x[IX(0, 0, N-1, N)]     = 0.33 * (x[IX(1, 0, N-1, N)]
                                  + x[IX(0, 1, N-1, N)]
                                  + x[IX(0, 0, N, N)]);
    x[IX(0, N-1, N-1, N)]   = 0.33 * (x[IX(1, N-1, N-1, N)]
                                  + x[IX(0, N-2, N-1, N)]
                                  + x[IX(0, N-1, N-2, N)]);
    x[IX(N-1, 0, 0, N)]     = 0.33 * (x[IX(N-2, 0, 0, N)]
                                  + x[IX(N-1, 1, 0, N)]
                                  + x[IX(N-1, 0, 1, N)]);
    x[IX(N-1, N-1, 0, N)]   = 0.33 * (x[IX(N-2, N-1, 0, N)]
                                  + x[IX(N-1, N-2, 0, N)]
                                  + x[IX(N-1, N-1, 1, N)]);
    x[IX(N-1, 0, N-1, N)]   = 0.33 * (x[IX(N-2, 0, N-1, N)]
                                  + x[IX(N-1, 1, N-1, N)]
                                  + x[IX(N-1, 0, N-2, N)]);
    x[IX(N-1, N-1, N-1, N)] = 0.33 * (x[IX(N-2, N-1, N-1, N)]
                                  + x[IX(N-1, N-2, N-1, N)]
                                  + x[IX(N-1, N-1, N-2, N)]);
}

void lin_solve(int b, inout float x[MAX_ARR_SIZE], inout float x0[MAX_ARR_SIZE], float a, float c, int iter, int N) {
    float cRecip = 1.0 / c;
    for (int k = 0; k < iter; k++) {
        for (int m = 1; m < N - 1; m++) {
            for (int j = 1; j < N - 1; j++) {
                for (int i = 1; i < N - 1; i++) {
                    x[IX(i, j, m, N)] =
                        (x0[IX(i, j, m, N)]
                            + a*(    x[IX(i+1, j  , m, N)]
                                    +x[IX(i-1, j  , m, N)]
                                    +x[IX(i  , j+1, m, N)]
                                    +x[IX(i  , j-1, m, N)]
                                    +x[IX(i  , j  , m+1, N)]
                                    +x[IX(i  , j  , m-1, N)]
                           )) * cRecip;
                }
            }
        }
        set_bnd(b, x, N);
    }
}

void diffuse(int b, inout float x[MAX_ARR_SIZE], inout float x0[MAX_ARR_SIZE], float diff, float dt, int iter, int N) {
    float a = dt * diff * float(N - 2) * float(N - 2);
    lin_solve(b, x, x0, a, 1. + 6. * a, iter, N);
}

void project(inout float velocX[MAX_ARR_SIZE], inout float velocY[MAX_ARR_SIZE], inout float velocZ[MAX_ARR_SIZE], inout float p[MAX_ARR_SIZE], inout float div[MAX_ARR_SIZE], int iter, int N) {
    for (int k = 1; k < N - 1; k++) {
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                div[IX(i, j, k, N)] = -0.5f*(
                         velocX[IX(i+1, j  , k  , N)]
                        -velocX[IX(i-1, j  , k  , N)]
                        +velocY[IX(i  , j+1, k  , N)]
                        -velocY[IX(i  , j-1, k  , N)]
                        +velocZ[IX(i  , j  , k+1, N)]
                        -velocZ[IX(i  , j  , k-1, N)]
                    )/float(N);
                p[IX(i, j, k, N)] = 0.;
            }
        }
    }
    set_bnd(0, div, N); 
    set_bnd(0, p, N);
    lin_solve(0, p, div, 1., 6., iter, N);
    
    for (int k = 1; k < N - 1; k++) {
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                velocX[IX(i, j, k, N)] -= 0.5* (  p[IX(i+1, j, k, N)]
                                                -p[IX(i-1, j, k, N)]) * float(N);
                velocY[IX(i, j, k, N)] -= 0.5 * (  p[IX(i, j+1, k, N)]
                                                -p[IX(i, j-1, k, N)]) * float(N);
                velocZ[IX(i, j, k, N)] -= 0.5 * (  p[IX(i, j, k+1, N)]
                                                -p[IX(i, j, k-1, N)]) * float(N);
            }
        }
    }
    set_bnd(1, velocX, N);
    set_bnd(2, velocY, N);
    set_bnd(3, velocZ, N);
}

void advect(int b, inout float d[MAX_ARR_SIZE], inout float d0[MAX_ARR_SIZE],  inout float velocX[MAX_ARR_SIZE], inout float velocY[MAX_ARR_SIZE], inout float velocZ[MAX_ARR_SIZE], float dt, int N)
{
    float i0, i1, j0, j1, k0, k1;
    
    float dtx = dt * float(N - 2);
    float dty = dt * float(N - 2);
    float dtz = dt * float(N - 2);
    
    float s0, s1, t0, t1, u0, u1;
    float tmp1, tmp2, tmp3, x, y, z;
    
    float Nfloat = float(N);
    float ifloat, jfloat, kfloat;
    int i, j, k;
    
    for(k = 1, kfloat = 1.; k < N - 1; k++, kfloat++) {
        for(j = 1, jfloat = 1.; j < N - 1; j++, jfloat++) { 
            for(i = 1, ifloat = 1.; i < N - 1; i++, ifloat++) {
                tmp1 = dtx * velocX[IX(i, j, k, N)];
                tmp2 = dty * velocY[IX(i, j, k, N)];
                tmp3 = dtz * velocZ[IX(i, j, k, N)];
                x    = ifloat - tmp1; 
                y    = jfloat - tmp2;
                z    = kfloat - tmp3;
                
                if(x < 0.5) x = 0.5; 
                if(x > Nfloat + 0.5) x = Nfloat + 0.5; 
                i0 = floor(x); 
                i1 = i0 + 1.0;
                if(y < 0.5) y = 0.5; 
                if(y > Nfloat + 0.5) y = Nfloat + 0.5; 
                j0 = floor(y);
                j1 = j0 + 1.0f; 
                if(z < 0.5) z = 0.5;
                if(z > Nfloat + 0.5) z = Nfloat + 0.5;
                k0 = floor(z);
                k1 = k0 + 1.0;
                
                s1 = x - i0; 
                s0 = 1.0 - s1; 
                t1 = y - j0; 
                t0 = 1.0 - t1;
                u1 = z - k0;
                u0 = 1.0 - u1;
                
                int i0i = int(i0);
                int i1i = int(i1);
                int j0i = int(j0);
                int j1i = int(j1);
                int k0i = int(k0);
                int k1i = int(k1);
                
                d[IX(i, j, k, N)] = 
                
                    s0 * ( t0 * (u0 * d0[IX(i0i, j0i, k0i, N)]
                                +u1 * d0[IX(i0i, j0i, k1i, N)])
                        +( t1 * (u0 * d0[IX(i0i, j1i, k0i, N)]
                                +u1 * d0[IX(i0i, j1i, k1i, N)])))
                   +s1 * ( t0 * (u0 * d0[IX(i1i, j0i, k0i, N)]
                                +u1 * d0[IX(i1i, j0i, k1i, N)])
                        +( t1 * (u0 * d0[IX(i1i, j1i, k0i, N)]
                                +u1 * d0[IX(i1i, j1i, k1i, N)])));
            }
        }
    }
    set_bnd(b, d, N);
}

