/**
 * @file nbi.c
 * @brief Functions for NBI simulation and particle generation
 */
#include <math.h>
#include "B_field.h"
#include "consts.h"
#include "math.h"
#include "nbi.h"
#include "random.h"
#include "plasma.h"
#include "suzuki.h"

void nbi_inject(nbi_injector* n, real* x, real* y, real* z, real* vx, real* vy,
                real* vz, real* anum, real* znum, random_data* rng) {
    int i_beamlet = random_uniform(rng) * n->n_beamlet;

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
    *vx = absv * n->beamlet_dirx[i_beamlet];
    *vy = absv * n->beamlet_diry[i_beamlet];
    *vz = absv * n->beamlet_dirz[i_beamlet];

    *anum = n->anum;
    *znum = n->znum;
}

void nbi_ionize(real* xyz, real* vxyz, int anum, int znum, B_field_data* Bdata, plasma_data* plsdata, random_data* rng) {
    real absv = math_norm(xyz);
    real energy = 0.5 * anum * CONST_U * absv*absv;
    real vhat[3];
    vhat[0] = vxyz[0] / absv;
    vhat[1] = vxyz[1] / absv;
    vhat[2] = vxyz[2] / absv;

    real remaining = 1;
    real threshold = random_uniform(rng);
    real ds = 1e-3;

    int n_species = plasma_get_n_species(plsdata);
    real* temp = (real*) malloc(n_species * sizeof(real));
    real* dens = (real*) malloc(n_species * sizeof(real));

    int* pls_anum = (int*) malloc(n_species * sizeof(int));;
    int* pls_znum = (int*) malloc(n_species * sizeof(int));;
    real* pls_mass = plasma_get_species_mass(plsdata);
    real* pls_charge = plasma_get_species_charge(plsdata);
    for(int i = 0; i < n_species; i++) {
        pls_anum[i] = round(pls_mass[i] / CONST_U);
        pls_znum[i] = round(pls_charge[i] / CONST_E);
    }

    while(remaining > threshold) {
        real rpz[3];
        math_xyz2rpz(xyz, rpz);

        real psi, rho;
        B_field_eval_psi(&psi, rpz[0], rpz[1], rpz[2], 0.0, Bdata);
        B_field_eval_rho(&rho, psi, Bdata);


        plasma_eval_densandtemp(dens, temp, rho, rpz[0], rpz[1], rpz[2], 0.0,
                                plsdata);
        real rate = dens[0] * 1e-4 * suzuki_sigmav(energy, dens[0], temp[0],
                                             n_species-1, &dens[1], anum, znum);
        xyz[0] += ds * vhat[0];
        xyz[1] += ds * vhat[1];
        xyz[2] += ds * vhat[2];
        remaining *= exp(-rate * ds);
    }
}
