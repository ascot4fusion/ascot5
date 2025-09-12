/**
 * @file bbnbi5.c
 * @brief Interface to BBNBI5
 *
 * BBNBI5 models neutral beam injectors and is used to evaluate shine-through
 * and beam birth-profile. Neutral markers are generated from injector geometry
 * and traced until they ionize or hit the wall. Several injectors can be
 * modelled simultaneously keeping in mind that in this case the output
 * the injector from which a particle originated is lost.
 */
#include <getopt.h>
#include <math.h>
#ifdef MPI
  #include <mpi.h>
#endif
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "ascot5.h"
#include "consts.h"
#include "gitver.h"
#include "math.h"
#include "physlib.h"
#include "simulate.h"
#include "endcond.h"
#include "random.h"
#include "particle.h"
#include "suzuki.h"
#include "B_field.h"
#include "plasma.h"
#include "neutral.h"
#include "wall.h"
#include "asigma.h"
#include "nbi.h"
#include "diag.h"
#include "bbnbi5.h"

void bbnbi_trace_markers(particle_queue *pq, sim_data* sim);
void bbnbi_inject_markers(particle_state* p, int nprt, int ngenerated, real t0,
                          real t1, nbi_injector* inj, sim_data* sim);

/**
 * @brief Simulate NBI injection
 *
 * This function initializes neutrals and traces them until they have ionized or
 * hit the wall.
 *
 * @param sim pointer to the simulation data structure
 * @param nprt number of markers to be injected
 * @param t1 time instant when the injector is turned on
 * @param t2 time instant when the injector is turned off
 * @param p pointer to the marker array which is allocated here
 */
void bbnbi_simulate(
    sim_data* sim, int nprt, real t1, real t2, particle_state** p) {

    /* Initialize input data */
    //simulate_init(sim);
    random_init(&sim->random_data, time(NULL));

    /* Calculate total NBI power so that we can distribute markers along
     * the injectors according to their power */
    real total_power = 0;
    for(int i=0; i < sim->nbi_data.ninj; i++) {
        total_power += sim->nbi_data.inj[i].power;
    }

    /* Initialize particle struct */
    *p = (particle_state*) malloc(nprt * sizeof(particle_state));

    /* Generate markers at the injectors */
    int nprt_generated = 0;
    for(int i = 0; i < sim->nbi_data.ninj; i++) {

        /* Number of markers generated is proportional to NBI power */
        int nprt_inj = ( sim->nbi_data.inj[i].power / total_power ) * nprt;
        if(i == sim->nbi_data.ninj-1) {
            /* All "remaining" markers goes to the last injector to avoid any
             * rounding issues */
            nprt_inj = nprt - nprt_generated;
        }

        /* Generates markers at the injector location and traces them until
         * they enter the region with magnetic field data */
        bbnbi_inject_markers(&((*p)[nprt_generated]), nprt_inj, nprt_generated,
                             t1, t2, &(sim->nbi_data.inj[i]), sim);

        nprt_generated += nprt_inj;
    }

    /* Place markers in a queue */
    particle_queue pq;
    pq.n = 0;
    for(int i = 0; i < nprt; i++) {
        pq.n++;
    }
    pq.p = (particle_state**) malloc(pq.n * sizeof(particle_state*));
    pq.finished = 0;

    pq.next = 0;
    for(int i = 0; i < nprt; i++) {
        pq.p[pq.next++] = &((*p)[i]);

    }
    pq.next = 0;

    /* Trace neutrals until they are ionized or lost to the wall */
    #pragma omp parallel
    bbnbi_trace_markers(&pq, sim);
}

/**
 * @brief Inject neutrals from an injector
 *
 * This function initializes neutral markers at the beamlet positions and
 * launches them in a (random) direction based on the injector specs.
 * The marker is traced until it enters the magnetic field, at which point
 * the particle struct is filled with only the particle data, and the struct
 * is returned.
 *
 * @param p pointer where generated markers are stored
 * @param nprt number of markers to be injected or generated
 * @param ngenerated number of markers that have already been generated
 * @param t0 time when the injector is turned on
 * @param t1 time when the injector is turned off
 * @param inj pointer to injector data
 * @param sim pointer to the sim struct with initialized data
 */
void bbnbi_inject_markers(particle_state* p, int nprt, int ngenerated, real t0,
                          real t1, nbi_injector* inj, sim_data* sim) {

    /* Set marker weights assuming a large number is created so that the energy
     * fractions of generated markers are close to the injector values */
    real f  =     1.0 * inj->efrac[0] + (1.0/2) * inj->efrac[1]
            + (1.0/3) * inj->efrac[2];
    real weight = (inj->power / inj->energy ) / ( f * nprt );

    /* Inject markers and trace their ballistic trajectories (without any
     * other physics) until they enter the plasma for the first time.     */
    #pragma omp parallel for
    for(int i = 0; i < nprt; i++) {
        real time = t0 + random_uniform(&sim->random_data) * (t1-t0);

        /* Assign initial phase-space coordinates for this marker */
        real xyz[3], vxyz[3], rpz[3], vhat[3];
        nbi_inject(xyz, vxyz, inj, &sim->random_data);
        math_xyz2rpz(xyz, rpz);
        math_unit(vxyz, vhat);

        /* Advance until the marker enters the magnetic field */
        real psi;
        real ds = 1e-3;
        a5err err = B_field_eval_psi(&psi, rpz[0], rpz[1], rpz[2], time,
                                     &sim->B_data);
        while(err) {
            xyz[0] += ds * vhat[0];
            xyz[1] += ds * vhat[1];
            xyz[2] += ds * vhat[2];
            math_xyz2rpz(xyz, rpz);
            err = B_field_eval_psi(&psi, rpz[0], rpz[1], rpz[2], time,
                                   &sim->B_data);
        }

        real vrpz[3];
        math_vec_xyz2rpz(vxyz, vrpz, rpz[1]);
        real gamma = physlib_gamma_vnorm(math_norm(vrpz));

        /* Fill the particle state with particle coordinates */
        p[i].rprt     = rpz[0];
        p[i].phiprt   = rpz[1];
        p[i].zprt     = rpz[2];
        p[i].p_r      = vrpz[0] * gamma * inj->mass;
        p[i].p_phi    = vrpz[1] * gamma * inj->mass;
        p[i].p_z      = vrpz[2] * gamma * inj->mass;
        p[i].mass     = inj->mass;
        p[i].charge   = 0.0;
        p[i].anum     = inj->anum;
        p[i].znum     = inj->znum;
        p[i].weight   = weight;
        p[i].time     = time;
        p[i].mileage  = 0.0;
        p[i].cputime  = 0.0;
        p[i].id       = ngenerated + i + 1;
        p[i].endcond  = 0;
        p[i].walltile = 0;
        p[i].err      = 0;
    }
}

/**
 * @brief Trace a neutral marker until it has ionized or hit wall
 *
 * This function is for the most part identical to simulate_fo with few
 * exceptions relevant for BBNBI.
 *
 * @param pq pointer to the marker queue containing the initial neutrals
 * @param sim pointer to the simu struct with initialized data
 */
void bbnbi_trace_markers(particle_queue *pq, sim_data* sim) {
    int cycle[NSIMD]  __memalign__;
    real hin[NSIMD]  __memalign__;
    int shinethrough[NSIMD] __memalign__;
    real remaining[NSIMD]  __memalign__;
    real threshold[NSIMD]  __memalign__;
    particle_simd_fo p, p0, pdiag;

    int n_species       = plasma_get_n_species(&sim->plasma_data);
    const int* pls_anum = plasma_get_species_anum(&sim->plasma_data);
    const int* pls_znum = plasma_get_species_znum(&sim->plasma_data);

    /* Init dummy markers */
    particle_allocate_fo(&p, NSIMD);
    particle_allocate_fo(&p0, NSIMD);
    particle_allocate_fo(&pdiag, NSIMD);
    for(int i=0; i< NSIMD; i++) {
        p.id[i] = -1;
        p.running[i] = 0;
        hin[i] = 1e-10;
        threshold[i] = random_uniform(&sim->random_data);
        remaining[i] = 1.0;
        shinethrough[i] = 0;
    }

    /* Initialize running particles */
    int n_running = particle_cycle_fo(pq, &p, &sim->B_data, cycle);
    while(n_running > 0) {

        #pragma omp simd
        for(int i=0; i< NSIMD; i++) {
            /* Store marker states */
            particle_copy_fo(&p, i, &p0, i);

            if(p.running[i]) {
                a5err err = 0;

                /* These are needed later */
                real pnorm = math_normc(p.p_r[i], p.p_phi[i], p.p_z[i]);
                real gamma = physlib_gamma_pnorm(p.mass[i], pnorm);
                real ekin  = physlib_Ekin_pnorm(p.mass[i], pnorm);

                /* Advance ballistic trajectory by converting momentum to
                 * cartesian coordinates */
                real prpz[3] = {p.p_r[i], p.p_phi[i], p.p_z[i]};
                real pxyz[3];
                math_vec_rpz2xyz(prpz, pxyz, p.phi[i]);

                real posrpz[3] = {p.r[i], p.phi[i], p.z[i]};
                real posxyz[3], fposxyz[3];
                math_rpz2xyz(posrpz, posxyz);
                fposxyz[0] = posxyz[0] + pxyz[0] * hin[i] / (gamma * p.mass[i]);
                fposxyz[1] = posxyz[1] + pxyz[1] * hin[i] / (gamma * p.mass[i]);
                fposxyz[2] = posxyz[2] + pxyz[2] * hin[i] / (gamma * p.mass[i]);

                /* Back to cylindrical coordinates (note phi is cumulative) */
                p.r[i] = sqrt(fposxyz[0]*fposxyz[0] + fposxyz[1]*fposxyz[1]);
                p.phi[i] += atan2(
                    posxyz[0] * fposxyz[1] - posxyz[1] * fposxyz[0],
                    posxyz[0] * fposxyz[0] + posxyz[1] * fposxyz[1] );
                p.z[i] = fposxyz[2];

                real cosp = cos(p.phi[i]);
                real sinp = sin(p.phi[i]);
                p.p_r[i]   =  pxyz[0] * cosp + pxyz[1] * sinp;
                p.p_phi[i] = -pxyz[0] * sinp + pxyz[1] * cosp;
                p.p_z[i]   =  pxyz[2];

                real ds = hin[i];
                p.mileage[i] += hin[i];

                /* Update background values at the new position */
                real psi, rho[2], pls_dens[MAX_SPECIES], pls_temp[MAX_SPECIES];
                err = B_field_eval_psi(
                    &psi, p.r[i], p.phi[i], p.z[i], p.time[i], &sim->B_data);
                if(!err) {
                    err = B_field_eval_rho(rho, psi, &sim->B_data);
                }
                if(!err && p.rho[i] <= 1.0 && rho[0] > 1.0) {
                    shinethrough[i] = 1;
                }
                p.rho[i] = rho[0];
                /* Update theta value */
                real axisrz[2];
                if(!err) {
                    B_field_get_axis_rz(axisrz, &sim->B_data, p.phi[i]);
                }
                p.theta[i] = atan2(p.z[i]-axisrz[1], p.r[i]-axisrz[0]);

                if(!err) {
                    err = plasma_eval_densandtemp(
                        pls_dens, pls_temp, rho[0], p.r[i], p.phi[i], p.z[i],
                        p.time[i], &sim->plasma_data);
                }

                /* Calculate ionization rate */
                real rate = 0.0;
                if(!err) {
                    real sigmav;
                    if( asigma_eval_bms(
                            &sigmav, p.znum[i], p.anum[i], ekin, p.mass[i],
                            n_species-1, pls_znum, pls_anum, pls_temp[0],
                            &(pls_dens[1]), &sim->asigma_data) ) {
                        err = 1;
                    }
                    rate = pls_dens[0] * sigmav;
                }
                remaining[i] *= exp(-rate * ds);

                /* Check for end conditions */
                if(!err) {
                    real w_coll = 0;
                    int tile = 0;
                    if(shinethrough[i]) {
                        tile = wall_hit_wall(
                            p0.r[i], p0.phi[i], p0.z[i],
                            p.r[i], p.phi[i], p.z[i], &sim->wall_data, &w_coll);
                    }
                    if(tile > 0) {
                        real w = w_coll;
                        p.time[i] = p0.time[i] + w*(p.time[i] - p0.time[i]);
                        p.r[i]    = p0.r[i]    + w*(p.r[i]    - p0.r[i]);
                        p.phi[i]  = p0.phi[i]  + w*(p.phi[i]  - p0.phi[i]);
                        p.z[i]    = p0.z[i]    + w*(p.z[i]    - p0.z[i]);

                        p.walltile[i] = tile;
                        p.endcond[i] |= endcond_wall;
                        p.running[i] = 0;
                    }
                    if(p.mileage[i] > NBI_MAX_DISTANCE) {
                        p.endcond[i] |= endcond_tlim;
                        p.running[i] = 0;
                    }
                    if(remaining[i] < threshold[i]) {
                        p.charge[i] = 1*CONST_E;
                        p.endcond[i] |= endcond_ioniz;
                        p.running[i] = 0;
                    }
                } else {
                    p.err[i] = err;
                    p.running[i] = 0;
                }
            }
        }

        /* Update markers that just finished */
        #pragma omp simd
        for(int i=0; i< NSIMD; i++) {
            /* Use this as a flag for which markers to update in diagnostics */
            pdiag.running[i] = 0;
            if(!p.running[i] && p.id[i] >= 0) {
                p.time[i] += p.mileage[i];

                /* Reset these for the next marker */
                threshold[i] = random_uniform(&sim->random_data);
                remaining[i] = 1.0;
                shinethrough[i] = 0;

                /* Update the magnetic field at the marker position */
                if(!p.err[i]) {
                    real B_dB[15];
                    B_field_eval_B_dB(B_dB, p.r[i], p.phi[i], p.z[i], p.time[i],
                                      &sim->B_data);
                    p.B_r[i]        = B_dB[0];
                    p.B_r_dr[i]     = B_dB[1];
                    p.B_r_dphi[i]   = B_dB[2];
                    p.B_r_dz[i]     = B_dB[3];

                    p.B_phi[i]      = B_dB[4];
                    p.B_phi_dr[i]   = B_dB[5];
                    p.B_phi_dphi[i] = B_dB[6];
                    p.B_phi_dz[i]   = B_dB[7];

                    p.B_z[i]        = B_dB[8];
                    p.B_z_dr[i]     = B_dB[9];
                    p.B_z_dphi[i]   = B_dB[10];
                    p.B_z_dz[i]     = B_dB[11];
                }
            }
            particle_copy_fo(&p, i, &pdiag, i);
            if(!p.running[i] && p.id[i] >= 0) {
                pdiag.running[i] = 1;
            }

            /* Normalize weight with time and add hin so that we don't divide
             * with zero when updating distributions */
            pdiag.time[i]   += hin[i];
            pdiag.weight[i] /= hin[i];
        }

        /* Update distributions for markers that finished */
        //diag_update_fo(&sim->diag_data, &sim->B_data, &pdiag, &p);

        /* Update running particles */
        n_running = particle_cycle_fo(pq, &p, &sim->B_data, cycle);
    }
}
