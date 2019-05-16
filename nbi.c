/**
 * @file nbi.c
 * @brief Functions for NBI simulation and particle generation
 */
#define _XOPEN_SOURCE 500 /**< drand48 requires POSIX 1995 standard */
#include <math.h>
#include "B_field.h"
#include "consts.h"
#include "math.h"
#include "nbi.h"
#include "particle.h"
#include "random.h"
#include "plasma.h"
#include "suzuki.h"

void nbi_inject(nbi_injector* n, real* x, real* y, real* z, real* vx, real* vy,
                real* vz, real* anum, real* znum, random_data* rng) {
    int i_beamlet = floor(random_uniform(rng) * n->n_beamlet);

    *x = n->beamlet_x[i_beamlet];
    *y = n->beamlet_y[i_beamlet];
    *z = n->beamlet_z[i_beamlet];

    real energy;
    real r = random_uniform(rng);
    if(r < n->efrac[0]) {
        energy = n->energy;
    } else if(r < n->efrac[0] + n->efrac[1]) {
        energy = n->energy / 2;
    } else {
        energy = n->energy / 3;
    }

    real absv = sqrt(2 * energy / (n->anum * CONST_U));
    *vx = absv * n->beamlet_dx[i_beamlet];
    *vy = absv * n->beamlet_dy[i_beamlet];
    *vz = absv * n->beamlet_dz[i_beamlet];

    *anum = n->anum;
    *znum = n->znum;
}

void nbi_ionize(real* xyz, real* vxyz, int* shinethrough, int anum, int znum,
                B_field_data* Bdata, plasma_data* plsdata, random_data* rng) {
    a5err err;

    real absv = math_norm(vxyz);
    real energy = 0.5 * anum * CONST_U * absv*absv / CONST_E / 1000;
    real vhat[3];
    vhat[0] = vxyz[0] / absv;
    vhat[1] = vxyz[1] / absv;
    vhat[2] = vxyz[2] / absv;

    real remaining = 1;
    real threshold = random_uniform(rng);
    real ds = 1e-3;

    int n_species = plasma_get_n_species(plsdata);
    const real* pls_mass = plasma_get_species_mass(plsdata);
    const real* pls_charge = plasma_get_species_charge(plsdata);

    real pls_temp[MAX_SPECIES], pls_dens[MAX_SPECIES];
    int pls_anum[MAX_SPECIES], pls_znum[MAX_SPECIES];

    for(int i = 0; i < n_species; i++) {
        pls_anum[i] = round(pls_mass[i] / CONST_U);
        pls_znum[i] = round(pls_charge[i] / CONST_E);
    }

    real s = 0.0;
    while(remaining > threshold && s < NBI_MAX_DISTANCE) {
        real rpz[3];
        math_xyz2rpz(xyz, rpz);

        real psi, rho;
        B_field_eval_psi(&psi, rpz[0], rpz[1], rpz[2], 0.0, Bdata);
        B_field_eval_rho(&rho, psi, Bdata);

        err = plasma_eval_densandtemp(pls_dens, pls_temp, rho, rpz[0], rpz[1],
                                      rpz[2], 0.0, plsdata);

        real rate;
        if(!err) {
            rate = pls_dens[0] * 1e-4*suzuki_sigmav(energy / anum, pls_dens[0],
                                                    pls_temp[0] / CONST_E,
                                                    n_species-1,
                                                    pls_dens+1, pls_anum+1,
                                                    pls_znum+1);
        }
        else {
            rate = 0.0; /* probably outside the plasma */
        }

        s += ds;
        xyz[0] += ds * vhat[0];
        xyz[1] += ds * vhat[1];
        xyz[2] += ds * vhat[2];
        remaining *= exp(-rate * ds);
    }

    if(s < NBI_MAX_DISTANCE) {
        *shinethrough = 0;
    }
    else {
        *shinethrough = 1;
    }
}

void nbi_generate(int nprt, particle_state* p, int* shinethrough,
                  nbi_injector* n, B_field_data* Bdata, plasma_data* plsdata,
                  random_data* rng) {
    for(int i = 0; i < nprt; i++) {
        real xyz[3], vxyz[3];
        real anum, znum;

        nbi_inject(n, &xyz[0], &xyz[1], &xyz[2], &vxyz[0], &vxyz[1], &vxyz[2],
                   &anum, &znum, rng);
        nbi_ionize(xyz, vxyz, shinethrough, anum, znum, Bdata, plsdata, rng);

        real rpz[3], vrpz[3];
        math_xyz2rpz(xyz, rpz);
        math_vec_xyz2rpz(vxyz, vrpz, rpz[1]);

        p[i].rprt = rpz[0];
        p[i].phiprt = rpz[1];
        p[i].zprt = rpz[2];
        p[i].rdot = vrpz[0];
        p[i].phidot = vrpz[1];
        p[i].zdot = vrpz[2];
    }
}
