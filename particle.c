/**
 * @file particle.c
 * @brief Particle representations and helper functions
 */

#include <math.h>
#include "ascot5.h"
#include "math.h"
#include "particle.h"
#include "B_field.h"

void particle_to_gc(particle* p, int i, particle_simd_gc* p_gc, int j,
                    B_field_data* Bdata) {
    real B[3];
    B_field_eval_B(B, p->r, p->phi, p->z, Bdata);
    real normB = sqrt(math_dot(B, B));

    real v[3];
    v[0] = p->rdot;
    v[1] = p->phidot;
    v[2] = p->zdot;
    real normv = sqrt(math_dot(v, v));

    real vcrossB[3];
    math_cross(v, B, vcrossB);

    p_gc->r[j] = p->r + p->mass * vcrossB[0]
                               / (p->charge * normB * normB);
    p_gc->phi[j] = p->phi + p->mass * vcrossB[1]
                                   / (p->charge * normB * normB);
    p_gc->z[j] = p->z + p->mass * vcrossB[2]
                               / (p->charge * normB * normB);
    p_gc->vpar[j] = math_dot(v, B) / normB;
    p_gc->mu[j] = p->mass / (2 * normB)
                          * (normv*normv - p_gc->vpar[j]*p_gc->vpar[j]);

    real ez[3];
    ez[0] = -B[2]/normB * B[0]/normB;
    ez[1] = -B[2]/normB * B[1]/normB;
    ez[2] = 1 - B[2]/normB * B[2]/normB;
    real rho[3];
    rho[0] = p->r - p_gc->r[j];
    rho[1] = p->phi - p_gc->phi[j];
    rho[2] = p->z - p_gc->z[j];
    real ezcrossrho[3];
    math_cross(ez, rho, ezcrossrho);
    p_gc->theta[j] = atan2(math_norm(ezcrossrho), math_dot(ez, rho));
/*    p_gc->theta[j] = acos(math_dot(ez, rho)
                          / (math_norm(ez)*math_norm(rho))); */

    p_gc->B_r[j] = B[0];
    p_gc->B_phi[j] = B[1];
    p_gc->B_z[j] = B[2];
    p_gc->mass[j] = p->mass;
    p_gc->charge[j] = p->charge;
    p_gc->weight[j] = p->weight;
    p_gc->time[j] = 0;
    p_gc->id[j] = p->id;
    p_gc->running[j] = p->running;
    p_gc->B_r[j] = B[0];
    p_gc->B_phi[j] = B[1];
    p_gc->B_z[j] = B[2];
    p_gc->index[j] = i;
}

void particle_to_gc_dummy(particle_simd_gc* p_gc, int j) {
    p_gc->r[j] = 1;
    p_gc->phi[j] = 1;
    p_gc->z[j] = 1;
    p_gc->vpar[j] = 1;
    p_gc->mu[j] = 1;
    p_gc->theta[j] = 1;
    p_gc->B_r[j] = 1;
    p_gc->B_phi[j] = 1;
    p_gc->B_z[j] = 1;
    p_gc->mass[j] = 1;
    p_gc->charge[j] = 1;
    p_gc->time[j] = 0;
    p_gc->weight[j] = 0;
    p_gc->id[j] = -1;
    p_gc->running[j] = 0;
    p_gc->index[j] = -1;
}

void gc_to_particle(particle_simd_gc* p_gc, int j, particle* p) {
    real B[3];
    B[0] = p_gc->B_r[j];
    B[1] = p_gc->B_phi[j];
    B[2] = p_gc->B_z[j];
    real normB = sqrt(math_dot(B, B));
    real Bxyz[3];
    math_vec_rpz2xyz(B, Bxyz, p_gc->phi[j]);

    real ez[3];
    ez[0] = -Bxyz[2]/normB * Bxyz[0]/normB;
    ez[1] = -Bxyz[2]/normB * Bxyz[1]/normB;
    ez[2] = 1 - Bxyz[2]/normB * Bxyz[2]/normB;

    real a[3];
    a[0] = ez[0]/math_norm(ez);
    a[1] = ez[1]/math_norm(ez);
    a[2] = ez[2]/math_norm(ez);
    
    real b[3];
    math_cross(B,ez,b);
    real tmp = math_norm(b);
    b[0] = b[0]/tmp;
    b[1] = b[1]/tmp;
    b[2] = b[2]/tmp;

    real rho = fabs(sqrt(2*p_gc->mu[j]*p_gc->mass[j]/normB)/p_gc->charge[j]);
    real rhovec[3];
    rhovec[0] = rho*(sin(p_gc->theta[j]) * a[0] + cos(p_gc->theta[j]) * b[0]);
    rhovec[1] = rho*(sin(p_gc->theta[j]) * a[1] + cos(p_gc->theta[j]) * b[1]);
    rhovec[2] = rho*(sin(p_gc->theta[j]) * a[2] + cos(p_gc->theta[j]) * b[2]);

    real rhovecrpz[3];
    math_vec_xyz2rpz(rhovec,rhovecrpz,p_gc->phi[j]);

/*    p->r = p_gc->r[j] + rhovecrpz[0];
    p->phi = p_gc->phi[j] + rhovecrpz[1];
    p->z = p_gc->z[j] + rhovecrpz[2];*/

    p->r = p_gc->r[j];
    p->phi = p_gc->phi[j];
    p->z = p_gc->z[j];
    p->rdot = p_gc->mu[j];
    p->phidot = p_gc->vpar[j];
    p->zdot = 0;
    p->time = p_gc->time[j];
    p->running = p_gc->running[j];
    p->endcond = p_gc->endcond[j];
    p->walltile = p_gc->walltile[j];
}
