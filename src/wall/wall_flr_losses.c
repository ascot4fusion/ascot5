#include "wall_flr_losses.h"
#include <stdlib.h>
#include <math.h>
#include "../consts.h"
#include "../gctransform.h"
#include "../B_field.h"
#include "../wall.h"
#include "../endcond.h"
#include "../error.h"

/**
 * @brief Evaluate FLR wall loss for a single guiding center marker.
 *
 * Samples a random gyrophase, performs first-order guiding-center to particle
 * transform, and checks if the line segment from GC position to particle
 * position intersects the wall. Does not mutate global particle SIMD arrays;
 * instead returns intersection info through output pointers.
 *
 * @param r      Guiding center R [m]
 * @param phi    Guiding center phi [rad]
 * @param z      Guiding center Z [m]
 * @param ppar   Parallel momentum [kg m/s]
 * @param mu     Magnetic moment [J/T]
 * @param mass   Particle mass [kg]
 * @param charge Particle charge [C]
 * @param time   Simulation time [s] (for time-dependent B)
 * @param B      Magnetic field data
 * @param wall   Wall data
 * @param rnd    Random generator state.
 * @param walltile_out  (output) Wall element ID if lost (0 otherwise)
 * @param err_out       (output) Error flag (set if B-field eval fails)
 * @return 1 if FLR loss occurred, 0 otherwise
 */
int flr_losses_eval(real r, real phi, real z, real ppar, real mu,
                    real mass, real charge, real time,
                    B_field_data* B, wall_data* wall, random_data* rnd,
                    int* walltile_out, int* err_out) {

    if(walltile_out) *walltile_out = 0;
    if(err_out) *err_out = 0;
 
    /* Sample random gyrophase */
    real zeta_rand;
    zeta_rand = random_uniform(rnd) * CONST_2PI;

    /* Evaluate magnetic field and gradients at GC position */
    real B_dB[15];
    a5err err = B_field_eval_B_dB(B_dB, r, phi, z, time, B);
    if(err) {
        if(err_out) *err_out = err;
        return 0; /* Cannot proceed */
    }

    /* Transform GC -> particle (first order) for sampled gyrophase */
    real rprt, phiprt, zprt;
    real pparprt, muprt, zetaprt; /* Unused outputs */
    gctransform_guidingcenter2particle(mass, charge, B_dB,
                                       r, phi, z, ppar, mu, zeta_rand,
                                       &rprt, &phiprt, &zprt,
                                       &pparprt, &muprt, &zetaprt);

    /* Check wall intersection */
    real w_coll = 0.0;
    int tile = wall_hit_wall(r, phi, z, rprt, phiprt, zprt, wall, &w_coll);
    if(tile > 0) {
        if(walltile_out) *walltile_out = tile;
        return 1;
    }
    return 0;
}