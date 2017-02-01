/**
 * @file test_particle.c
 * @brief Test program for particle structs and functions
 */
#include <stdio.h>
#include "B_GS.h"
#include "particle.h"
#include "math.h"

int main(int argc, char** argv) {
    B_GS_offload_data offload_data;
    real* offload_array;
    B_GS_init_offload(&offload_data, &offload_array);
    B_GS_data Bdata;
    B_GS_init(&Bdata, &offload_data, offload_array);

    particle p;
    p.r = 6;
    p.phi = 0.0;
    p.z = -0.2;
    p.rdot = 1.6e6;
    p.phidot = -7.5e6;
    p.zdot = -1.0e7;
    p.mass = 4*1.66e-27;
    p.charge = 2*1.602e-19;
    p.weight = 1.0;
    p.time = 0.0;
    p.id = 0;
    p.running = 1;
    p.endcond = 0;
    p.walltile = -1;
    fprintf(stderr, "prpz: %lf %lf %lf\n", p.r, p.phi, p.z);
    particle_simd_gc p_gc;
    particle_to_gc(&p, &p_gc, 0, &Bdata);
    fprintf(stderr, "Brpz: %lf %lf %lf\n", p_gc.B_r[0], p_gc.B_phi[0],
            p_gc.B_z[0]);
    fprintf(stderr, "pgcrpz: %lf %lf %lf\n", p_gc.r[0], p_gc.phi[0], p_gc.z[0]);
    gc_to_particle(&p_gc, 0, &p, &Bdata);
    fprintf(stderr, "prpz:    %lf %lf %lf\n", p.r, p.phi, p.z);
    real theta;
    for(theta = 0; theta < 2*math_pi; theta+=0.1) {
        p_gc.theta[0] = theta;
        gc_to_particle(&p_gc, 0, &p, &Bdata);
        printf("%lf %lf %lf\n", p.r, p.phi, p.z);
    }
    return 0;
}
