#include "scene.glsl"

void main() {
    Scene scene = sc_init();
    vec3 col = sc_col(scene);
    gl_FragColor = vec4(col, 1.);
}