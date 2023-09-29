/**
 * @file nbi.c
 * @brief Functions for NBI simulation and particle generation
 */
#define _XOPEN_SOURCE 500 /**< drand48 requires POSIX 1995 standard */
#include <math.h>
#include "print.h"
#include "ascot5.h"
#include "consts.h"
#include "math.h"
#include "physlib.h"
#include "particle.h"
#include "random.h"
#include "suzuki.h"
#include "B_field.h"
#include "plasma.h"
#include "wall.h"
#include "nbi.h"

#include "diag.h"
#include "diag/dist_5D.h"

void nbi_inject(real* xyz, real* vxyz, real* efrac, nbi_injector* inj,
                random_data* rng);

void nbi_ionize(real* xyz, real* vxyz, real time, int* shinethrough, int anum,
                int znum, real mass, B_field_data* Bdata, plasma_data* plsdata,
                wall_data* walldata, random_data* rng);

/**
 * @brief Load NBI data and prepare parameters for offload.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 *
 * @return zero if initialization succeeded.
 */
int nbi_init_offload(nbi_offload_data* offload_data, real** offload_array) {
    int err = 0;
    print_out(VERBOSE_IO, "\nNBI input\n");
    print_out(VERBOSE_IO, "Number of injectors %d:\n", offload_data->ninj);
    for(int i=0; i<offload_data->ninj; i++) {
        print_out(VERBOSE_IO, "\n  Injector ID %d (%d beamlets) Power: %1.1e\n",
                  offload_data->id[i], offload_data->n_beamlet[i],
                  offload_data->power[i]);
        print_out(VERBOSE_IO,
                  "    Anum %d Znum %d mass %1.1e amu energy %1.1e eV\n",
                  offload_data->anum[i], offload_data->znum[i],
                  offload_data->mass[i] / CONST_U,
                  offload_data->energy[i] / CONST_E);
        print_out(VERBOSE_IO,
                  "    Energy fractions: %1.1e (Full) %1.1e (1/2) %1.1e (1/3)\n",
                  offload_data->efrac[0], offload_data->efrac[1],
                  offload_data->efrac[2]);

        /* Even if halo fraction is zero, the divergences should be nonzero
           to avoid division by zero during evaluation. Do this after the
           input has been printed as to not confuse the user */
        if(offload_data->div_halo_frac[i] == 0) {
            offload_data->div_halo_h[i] = 1e-10;
            offload_data->div_halo_v[i] = 1e-10;
        }
    }
    return err;
}

/**
 * @brief Initialize NBI data struct on target
 *
 * @param nbi pointer to data struct on target
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array
 */
void nbi_init(nbi_data* nbi, nbi_offload_data* offload_data,
              real* offload_array) {
    int idx = 0;
    nbi->ninj = offload_data->ninj;
    for(int i=0; i<nbi->ninj; i++) {
        nbi->inj[i].anum          = offload_data->anum[i];
        nbi->inj[i].znum          = offload_data->znum[i];
        nbi->inj[i].mass          = offload_data->mass[i];
        nbi->inj[i].power         = offload_data->power[i];
        nbi->inj[i].energy        = offload_data->energy[i];
        nbi->inj[i].efrac[0]      = offload_data->efrac[3*i+0];
        nbi->inj[i].efrac[1]      = offload_data->efrac[3*i+1];
        nbi->inj[i].efrac[2]      = offload_data->efrac[3*i+2];
        nbi->inj[i].div_h         = offload_data->div_h[i];
        nbi->inj[i].div_v         = offload_data->div_v[i];
        nbi->inj[i].div_halo_frac = offload_data->div_halo_frac[i];
        nbi->inj[i].div_halo_h    = offload_data->div_halo_h[i];
        nbi->inj[i].div_halo_v    = offload_data->div_halo_v[i];
        nbi->inj[i].id            = offload_data->id[i];
        nbi->inj[i].n_beamlet     = offload_data->n_beamlet[i];

        int n_beamlet = nbi->inj[i].n_beamlet;
        nbi->inj[i].beamlet_x  = &(offload_array[idx + 0*n_beamlet]);
        nbi->inj[i].beamlet_y  = &(offload_array[idx + 1*n_beamlet]);
        nbi->inj[i].beamlet_z  = &(offload_array[idx + 2*n_beamlet]);
        nbi->inj[i].beamlet_dx = &(offload_array[idx + 3*n_beamlet]);
        nbi->inj[i].beamlet_dy = &(offload_array[idx + 4*n_beamlet]);
        nbi->inj[i].beamlet_dz = &(offload_array[idx + 5*n_beamlet]);
        idx += 6 * n_beamlet;
    }
}

/**
 * @brief Free offload array
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
void nbi_free_offload(nbi_offload_data* offload_data,
                      real** offload_array) {
    for(int i=0; i<offload_data->ninj; i++) {
        offload_data->n_beamlet[i] = 0;
    }
    offload_data->ninj = 0;
    free(*offload_array);
}

/**
 * @brief Generate NBI ions by injecting and tracing neutrals from an injector
 *
 * @param pout pointer where generated markers are stored or NULL
 * @param nprt number of markers to be injected or generated
 * @param t0 time when the injector is turned on
 * @param t1 time when the injector is turned off
 * @param inj pointer to injector data
 * @param Bdata pointer to magnetic field data
 * @param plasma_data pointer to plasma data
 * @param wall_data pointer to wall data
 * @param random_data pointer to random number generator data
 * @param diag pointer to diagnostics data
 */
void nbi_generate(particle* pout, int nprt, real t0, real t1, nbi_injector* inj,
                  B_field_data* Bdata, plasma_data* plasma_data,
                  wall_data* wall_data, random_data* random_data,
                  diag_data* diag) {

    int nelements = wall_get_n_elements(wall_data);
    real* eload   = (real*)malloc( nelements * sizeof(real) );

    /* One marker contributes to a distribution as weight*dt. If beams are on
     * for a set time, use the time interval as dt. If beams are on indefinitely
     * (t1=t0) then use dt=1.0 s since then we don't need to normalize the
     * distribution afterwards */
    real dt = t1 - t0;
    if(dt == 0.0) {
        dt = 1.0;
    }
    /* Set marker weights assuming a large number is created so that the energy
     * fractions of generated markers are close to the injector values */
    real f  =     1.0 * inj->efrac[0] + (1.0/2) * inj->efrac[1]
            + (1.0/3) * inj->efrac[2];
    real weight = (inj->power / inj->energy ) / ( f * nprt );

    int ngenerated = 0;
    #pragma omp parallel for
    for(int i = 0; i < nprt; i++) {
        particle ptemp, *p;
        real xyz[3], vxyz[3];
        real charge, gamma, efrac;
        int anum  = inj->anum;
        int znum  = inj->znum;
        real mass = inj->mass;
        real time = t0 + random_uniform(rng) * (t1-t0);

        int walltile = -1;
        nbi_inject(xyz, vxyz, &efrac, inj, random_data);
        nbi_ionize(xyz, vxyz, time, &walltile, anum, znum, mass, Bdata,
                   plasma_data, wall_data, random_data);

        if(walltile > 0) {
            eload[walltile] += weight;
        }
        else {
            int idx;
            #pragma omp critical
            {
                idx = ngenerated;
                ngenerated++;
            }

            real rpz[3], vrpz[3];
            math_xyz2rpz(xyz, rpz);
            math_vec_xyz2rpz(vxyz, vrpz, rpz[1]);

            /* Convert state to input and store in output array if requested */
            gamma = physlib_gamma_vnorm(math_norm(vrpz));
            charge = 1 * CONST_E; // Singly ionized always
            if(pout != NULL) {
                pout[idx].r      = rpz[0];
                pout[idx].phi    = rpz[1];
                pout[idx].z      = rpz[2];
                pout[idx].p_r    = vrpz[0] * gamma * mass;
                pout[idx].p_phi  = vrpz[1] * gamma * mass;
                pout[idx].p_z    = vrpz[2] * gamma * mass;
                pout[idx].anum   = anum;
                pout[idx].znum   = znum;
                pout[idx].charge = charge;
                pout[idx].mass   = mass;
                pout[idx].id     = i+1;
                pout[idx].time   = time;
                pout[idx].weight = weight;
                p = &(pout[idx]);
            }
            else {
                ptemp.r      = rpz[0];
                ptemp.phi    = rpz[1];
                ptemp.z      = rpz[2];
                ptemp.p_r    = vrpz[0] * gamma * mass;
                ptemp.p_phi  = vrpz[1] * gamma * mass;
                ptemp.p_z    = vrpz[2] * gamma * mass;
                ptemp.anum   = anum;
                ptemp.znum   = znum;
                ptemp.charge = charge;
                ptemp.mass   = mass;
                ptemp.id     = i+1;
                ptemp.time   = time;
                ptemp.weight = weight;
                p = &ptemp;
            }

            /* Convert input to state and state to sim. struct for diagnostics*/
            particle_state ps;
            particle_simd_fo pi, pf;
            particle_input_p_to_state(p, &ps, Bdata);
            for(int j=0; j<NSIMD; j++) {
                pi.id[j] = -1;
                pf.id[j] = -1;
            }
            particle_state_to_fo(&ps, 0, &pi, 0, Bdata);
            particle_copy_fo(&pi, 0, &pf, 0);
            pf.time[0] = 1.0;

            diag_update_fo(diag, Bdata, &pf, &pi);
        }
    }

    /* Second loop but this time we don't update diagnostics. Instead we
     * generate markers until the particle struct is filled (skip this loop
     * if markers are not requested). We have do it like this because we can
     * either fix the number of markers to be generated or the marker weight
     * before the simulation loop but not both. Last loop fixed weight since
     * it is required by the diagnostics (markers generated in last loop need
     * to be reweighted). */
    int nionized = 0;
    int nshined  = 0;
    if(pout == NULL) { ngenerated = nprt; }
    #pragma omp parallel for
    for(int i = ngenerated; i < nprt; i++) {
        real xyz[3], vxyz[3];
        real charge, gamma, efrac;
        int anum  = inj->anum;
        int znum  = inj->znum;
        real mass = inj->mass;
        real time = t0 + random_uniform(rng) * (t1-t0);

        int shinethrough = -1;
        do {
            nbi_inject(xyz, vxyz, &efrac, inj, random_data);
            nbi_ionize(xyz, vxyz, time, &shinethrough, anum, znum, mass, Bdata,
                       plasma_data, wall_data, random_data);
            if(shinethrough > 0) {
                #pragma omp atomic
                nshined++;
            }
        } while(shinethrough != 0);
        real rpz[3], vrpz[3];
        math_xyz2rpz(xyz, rpz);
        math_vec_xyz2rpz(vxyz, vrpz, rpz[1]);

        gamma  = physlib_gamma_vnorm(math_norm(vrpz));
        charge = 1 * CONST_E; // Singly ionized always
        pout[i].r      = rpz[0];
        pout[i].phi    = rpz[1];
        pout[i].z      = rpz[2];
        pout[i].p_r    = vrpz[0] * gamma * mass;
        pout[i].p_phi  = vrpz[1] * gamma * mass;
        pout[i].p_z    = vrpz[2] * gamma * mass;
        pout[i].anum   = anum;
        pout[i].znum   = znum;
        pout[i].charge = charge;
        pout[i].mass   = mass;
        pout[i].id     = i+1;
        pout[i].time   = time;
        pout[i].weight = weight;

        #pragma omp atomic
        nionized++;
    }
    weight *= ((real) nprt) / (nprt + nshined + nionized);
    for(int i = 0; i < nprt; i++) {
        pout[i].weight = weight;
    }
}

/**
 * @brief Initialize a neutral marker from an injector.
 *
 * @param inj pointer to injector data
 * @param xyz initialized marker's position in cartesian coordinates [m]
 * @param vxyz initialized marker's velocity in cartesian coordinates [m/s]
 * @param efrac what fraction of injector energy was assigned to this marker
 * @param rng pointer to random number generator data
 */
void nbi_inject(real* xyz, real* vxyz, real* efrac, nbi_injector* inj,
                random_data* rng) {
    /* Pick a random beamlet and initialize marker there */
    int i_beamlet = floor(random_uniform(rng) * inj->n_beamlet);
    xyz[0] = inj->beamlet_x[i_beamlet];
    xyz[1] = inj->beamlet_y[i_beamlet];
    xyz[2] = inj->beamlet_z[i_beamlet];

    /* Pick marker energy based on the energy fractions */
    real energy;
    real r = random_uniform(rng);
    if(r < inj->efrac[0]) {
        energy = inj->energy;
        *efrac = 1.0;
    } else if(r < inj->efrac[0] + inj->efrac[1]) {
        energy = inj->energy / 2;
        *efrac = 1.0/2;
    } else {
        energy = inj->energy / 3;
        *efrac = 1.0/3;
    }

    /* Calculate vertical and horizontal normals for beam divergence */
    real dir[3], normalv[3], normalh[3], tmp[3];
    dir[0] = inj->beamlet_dx[i_beamlet];
    dir[1] = inj->beamlet_dy[i_beamlet];
    dir[2] = inj->beamlet_dz[i_beamlet];

    real phi   = atan2(dir[1], dir[0]);
    real theta = acos(dir[2]);

    normalv[0] = sin(theta+CONST_PI/2) * cos(phi);
    normalv[1] = sin(theta+CONST_PI/2) * sin(phi);
    normalv[2] = cos(theta+CONST_PI/2);

    math_cross(dir, normalv, tmp);
    math_unit(tmp, normalh);

    /* Assuming isotropic divergence using horizontal value,
       ignoring halo divergence for now */
    r = random_uniform(rng);
    real div = inj->div_h * sqrt(log(1/(1-r)));

    real angle = random_uniform(rng) * 2 * CONST_PI;
    tmp[0] = dir[0] + div * ( cos(angle)*normalh[0] + sin(angle)*normalv[0] );
    tmp[1] = dir[1] + div * ( cos(angle)*normalh[1] + sin(angle)*normalv[1] );
    tmp[2] = dir[2] + div * ( cos(angle)*normalh[2] + sin(angle)*normalv[2] );

    math_unit(tmp, dir);

    real gamma = physlib_gamma_Ekin(inj->mass, energy);
    real absv  = sqrt( 1.0 - 1.0 / (gamma * gamma) ) * CONST_C;
    vxyz[0] = absv * dir[0];
    vxyz[1] = absv * dir[1];
    vxyz[2] = absv * dir[2];
}

/**
 * @brief Trace a neutral marker until it has ionized or hit wall
 *
 * This function updates marker's position and velocity coordinates that are
 * given as an input.
 *
 * @param xyz marker initial position in cartesian coordinates [m]
 * @param vxyz marker initial velocity vector in cartesian coordinates [m/s]
 * @param time time instance when marker was born [s]
 * @param shinethrough flag indicating if the marker hit the wall and which tile
 * @param anum marker atomic mass number
 * @param znum marker charge number
 * @param mass marker mass
 * @param Bdata pointer to magnetic field data
 * @param plsdata pointer to plasma data
 * @param walldata pointer to wall data
 * @param rng pointer to random number generator data
 */
void nbi_ionize(real* xyz, real* vxyz, real time, int* shinethrough, int anum,
                int znum, real mass, B_field_data* Bdata, plasma_data* plsdata,
                wall_data* walldata, random_data* rng) {
    a5err err;

    real vnorm = math_norm(vxyz);
    real gamma = physlib_gamma_vnorm(vnorm);
    real ekin  = physlib_Ekin_gamma(mass, gamma);
    real vhat[3];
    vhat[0] = vxyz[0] / vnorm;
    vhat[1] = vxyz[1] / vnorm;
    vhat[2] = vxyz[2] / vnorm;

    int n_species       = plasma_get_n_species(plsdata);
    const int* pls_anum = plasma_get_species_anum(plsdata);
    const int* pls_znum = plasma_get_species_znum(plsdata);
    real pls_temp[MAX_SPECIES], pls_dens[MAX_SPECIES];

    real remaining = 1;
    real threshold = random_uniform(rng);
    real ds = 1e-3;

    real s = 0.0;
    int entered_plasma = 0;
    int exited_plasma  = 0;

    while(remaining > threshold && s < NBI_MAX_DISTANCE) {
        err = 0;
        real rate = 0.0;

        real rpz[3];
        math_xyz2rpz(xyz, rpz);

        real psi, rho[2];
        err = B_field_eval_psi(&psi, rpz[0], rpz[1], rpz[2], 0.0, Bdata);

        if(!err) {
            err = B_field_eval_rho(rho, psi, Bdata);

            /* Check for wall collisions after passing through separatrix
             * twice */
            if(!entered_plasma && rho[0] < 1.0) {
                entered_plasma = 1;
            }
            if(entered_plasma && !exited_plasma && rho[0] >= 1.0) {
                exited_plasma = 1;
            }

            err = plasma_eval_densandtemp(pls_dens, pls_temp, rho[0], rpz[0],
                                          rpz[1], rpz[2], time, plsdata);

            if(!err) {
                real sigmav;
                if( suzuki_sigmav(
                        &sigmav, ekin / anum, pls_dens[0], pls_temp[0],
                        n_species-1, &(pls_dens[1]), pls_anum, pls_znum) ) {
                    err = 1;
                }
                rate = pls_dens[0] * sigmav;
            }
            else {
                rate = 0.0; /* Outside the plasma */
            }
        }
        else {
            rate = 0.0; /* Outside the magnetic field */
        }

        s += ds;
        xyz[0] += ds * vhat[0];
        xyz[1] += ds * vhat[1];
        xyz[2] += ds * vhat[2];
        remaining *= exp(-rate * ds);

        if(exited_plasma) {
            real rpz2[3]; /* New position, old position already in rpz */
            math_xyz2rpz(xyz, rpz2);
            real w_coll;
            int tile = wall_hit_wall(rpz[0], rpz[1], rpz[2], rpz2[0], rpz2[1],
                                     rpz2[2], walldata, &w_coll);

            if(tile > 0) {
                /* Hit wall */
                *shinethrough = tile;
                return;
            }
        }
    }

    /* Check if we've missed the plasma entirely */
    if(s < NBI_MAX_DISTANCE) {
        *shinethrough = 0;
    }
    else {
        *shinethrough = -1;
    }
}
