/**
 * @file nbi.c
 * @brief Functions for NBI simulation and particle generation
 */
#define _XOPEN_SOURCE 500 /**< drand48 requires POSIX 1995 standard */
#include <math.h>
#include "B_field.h"
#include "consts.h"
#include "physlib.h"
#include "math.h"
#include "nbi.h"
#include "particle.h"
#include "random.h"
#include "plasma.h"
#include "suzuki.h"
#include "wall.h"

void nbi_inject(nbi_injector* n, real* x, real* y, real* z, real* vx, real* vy,
                real* vz, int* anum, int* znum, real* mass, random_data* rng) {
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

    /* calculate vertical and horizontal normals for divergence */
    real dir[3], normalv[3], normalh[3], tmp[3];

    dir[0] = n->beamlet_dx[i_beamlet];
    dir[1] = n->beamlet_dy[i_beamlet];
    dir[2] = n->beamlet_dz[i_beamlet];

    real phi = atan2(dir[1], dir[0]);
    real theta = acos(dir[2]);

    normalv[0] = sin(theta+CONST_PI/2) * cos(phi);
    normalv[1] = sin(theta+CONST_PI/2) * sin(phi);
    normalv[2] = cos(theta+CONST_PI/2);

    math_cross(dir, normalv, tmp);
    math_unit(tmp, normalh);

    /* assuming isotropic divergence using horizontal value,
       ignoring halo divergence for now */
    r = random_uniform(rng);
    real div = n->div_h * sqrt(log(1/(1-r)));

    real angle = random_uniform(rng) * 2 * CONST_PI;
    tmp[0] = dir[0] + div * (cos(angle)*normalh[0] + sin(angle)*normalv[0]);
    tmp[1] = dir[1] + div * (cos(angle)*normalh[1] + sin(angle)*normalv[1]);
    tmp[2] = dir[2] + div * (cos(angle)*normalh[2] + sin(angle)*normalv[2]);

    math_unit(tmp, dir);

    real absv = sqrt(2 * energy / (n->anum * CONST_U));
    *vx = absv * dir[0];
    *vy = absv * dir[1];
    *vz = absv * dir[2];

    *anum = n->anum;
    *znum = n->znum;
    *mass = n->mass;
}

void nbi_ionize(real* xyz, real* vxyz, int* shinethrough, int anum, int znum,
                B_field_data* Bdata, plasma_data* plsdata, wall_data* walldata,
                random_data* rng) {
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
    int entered_plasma = 0;
    int exited_plasma = 0;

    while(remaining > threshold && s < NBI_MAX_DISTANCE) {
        err = 0;
        real rate = 0.0;

        real rpz[3];
        math_xyz2rpz(xyz, rpz);

        real psi, rho;
        err = B_field_eval_psi(&psi, rpz[0], rpz[1], rpz[2], 0.0, Bdata);

        if(!err) {
            err = B_field_eval_rho(&rho, psi, Bdata);

            /* check for wall collisions after passing through separatrix
             * twice */
            if(!entered_plasma && rho < 1.0) {
                entered_plasma = 1;
            }
            if(entered_plasma && !exited_plasma && rho >= 1.0) {
                exited_plasma = 1;
            }

            err = plasma_eval_densandtemp(pls_dens, pls_temp, rho, rpz[0],
                                          rpz[1], rpz[2], 0.0, plsdata);

            if(!err) {
                rate = pls_dens[0] * 1e-4*suzuki_sigmav(energy / anum,
                                                        pls_dens[0],
                                                        pls_temp[0] / CONST_E,
                                                        n_species-1,
                                                        pls_dens+1,
                                                        pls_anum+1,
                                                        pls_znum+1);
            }
            else {
                rate = 0.0; /* outside the plasma */
            }
        }
        else {
            rate = 0.0; /* outside the magnetic field */
        }

        s += ds;
        xyz[0] += ds * vhat[0];
        xyz[1] += ds * vhat[1];
        xyz[2] += ds * vhat[2];
        remaining *= exp(-rate * ds);

        if(entered_plasma) {
            real rpz2[3]; /* new position, old position already in rpz */
            math_xyz2rpz(xyz, rpz2);
	    real w_coll = 0;
            int tile = wall_hit_wall(rpz[0], rpz[1], rpz[2], rpz2[0], rpz2[1],
                                     rpz2[2], walldata, &w_coll);

            if(tile > 0) {
                /* hit wall */
                *shinethrough = 1;
                return;
            }
        }
    }

    /* check if we've missed the plasma entirely */
    if(s < NBI_MAX_DISTANCE) {
        *shinethrough = 0;
    }
    else {
        *shinethrough = 1;
    }
}

void nbi_generate(int nprt, particle* p, nbi_injector* n,
                  B_field_data* Bdata, plasma_data* plsdata,
                  wall_data* walldata, random_data* rng) {

    real totalShine = 0.0;
    real totalPower = 0.0;

    #pragma omp parallel for
    for(int i = 0; i < nprt; i++) {
        real xyz[3], vxyz[3];
        int anum, znum;
        real mass;

        int shinethrough = 1;
        do {
            nbi_inject(n, &xyz[0], &xyz[1], &xyz[2], &vxyz[0], &vxyz[1],
                       &vxyz[2], &anum, &znum, &mass, rng);
            nbi_ionize(xyz, vxyz, &shinethrough, anum, znum, Bdata, plsdata,
                   walldata, rng);
            totalPower += 0.5 * mass * pow(math_norm(vxyz), 2);
            if(shinethrough == 1) {
                #pragma omp atomic
                totalShine += 0.5 * mass * pow(math_norm(vxyz), 2);

            }
        } while(shinethrough == 1);

        real rpz[3], vrpz[3];
        math_xyz2rpz(xyz, rpz);
        math_vec_xyz2rpz(vxyz, vrpz, rpz[1]);

        real gamma = physlib_gamma_vnorm(math_norm(vrpz));
        p[i].r      = rpz[0];
        p[i].phi    = rpz[1];
        p[i].z      = rpz[2];
        p[i].p_r    = vrpz[0] * gamma * mass;
        p[i].p_phi  = vrpz[1] * gamma * mass;
        p[i].p_z    = vrpz[2] * gamma * mass;
        p[i].anum   = anum;
        p[i].znum   = znum;
        p[i].charge = 1 * CONST_E;
        p[i].mass   = mass;
        p[i].id     = i+1;
        p[i].time   = 0;
    }

    for(int i = 0; i < nprt; i++) {
        p[i].weight = n->power * (1 - totalShine/totalPower) / nprt;
    }
}
