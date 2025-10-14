/**
 * Simulate guiding centers using adaptive time-step (see simulate.h).
 */
#include "consts.h"
#include "coulomb_collisions.h"
#include "data/bfield.h"
#include "data/boozer.h"
#include "data/diag.h"
#include "data/efield.h"
#include "data/marker.h"
#include "data/mhd.h"
#include "data/plasma.h"
#include "data/rfof.h"
#include "data/wall.h"
#include "datatypes.h"
#include "defines.h"
#include "endcond.h"
#include "orbit_following.h"
#include "simulate.h"
#include "utils/mathlib.h"
#include "utils/physlib.h"
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

DECLARE_TARGET_SIMD_UNIFORM(sim)
real simulate_gc_adaptive_inidt(Simulation *sim, MarkerGuidingCenter *p, int i);

#define DUMMY_TIMESTEP_VAL 1.0 /**< Dummy time step value */

void simulate_gc_adaptive(Simulation *sim, MarkerQueue *pq, size_t vector_size)
{

    /* Wiener arrays needed for the adaptive time step */
    mccc_wienarr wienarr[NSIMD];

    /* Current time step, suggestions for the next time step and next time
     * step                                                                */
    real *hin = (real *)malloc(vector_size * sizeof(real));
    real *hout_orb = (real *)malloc(vector_size * sizeof(real));
    real *hout_col = (real *)malloc(vector_size * sizeof(real));
    real *hout_rfof = (real *)malloc(vector_size * sizeof(real));
    real *hnext = (real *)malloc(vector_size * sizeof(real));

    /* Flag indicateing whether a new marker was initialized */
    size_t *cycle = (size_t *)malloc(vector_size * sizeof(size_t));

    real tol_col = sim->options->adaptive_tolerance_collisions;
    real tol_orb = sim->options->adaptive_tolerance_orbit;

    real cputime, cputime_last; // Global cpu time: recent and previous record

    MarkerGuidingCenter p;  // This array holds current states
    MarkerGuidingCenter p0; // This array stores previous states
    marker_allocate_gc(&p, vector_size);
    marker_allocate_gc(&p0, vector_size);

    rfof_marker rfof_mrk; // RFOF specific data

    for (size_t i = 0; i < vector_size; i++)
    {
        p.id[i] = 0;
        p.running[i] = 0;
    }

    /* Initialize running particles */
    size_t n_running = marker_cycle_gc(pq, &p, &sim->bfield, cycle);

    if (sim->options->enable_icrh)
    {
        rfof_set_up(&rfof_mrk, &sim->rfof);
    }

#pragma omp simd
    for (size_t i = 0; i < vector_size; i++)
    {
        if (cycle[i] > 0)
        {
            /* Determine initial time-step */
            hin[i] = simulate_gc_adaptive_inidt(sim, &p, i);
            if (sim->options->enable_coulomb_collisions)
            {
                /* Allocate array storing the Wiener processes */
                mccc_wiener_initialize(&(wienarr[i]), p.time[i]);
            }
        }
    }

    cputime_last = A5_WTIME;

    /* MAIN SIMULATION LOOP
     * - Store current state
     * - Integrate motion due to bacgkround EM-field (orbit-following)
     * - Integrate scattering due to Coulomb collisions
     * - Check whether time step was accepted
     *   - NO:  revert to initial state and ignore the end of the loop
     *          (except CPU_TIME_MAX end condition if this is implemented)
     *   - YES: update particle time, clean redundant Wiener processes, and
     *          proceed
     * - Check for end condition(s)
     * - Update diagnostics
     */
    real *rnd = (real *)malloc(5 * vector_size * sizeof(real));
    while (n_running > 0)
    {

/* Store marker states in case time step will be rejected */
#pragma omp simd
        for (size_t i = 0; i < vector_size; i++)
        {
            marker_copy_gc(&p, i, &p0, i);
            hout_orb[i] = DUMMY_TIMESTEP_VAL;
            hout_col[i] = DUMMY_TIMESTEP_VAL;
            hout_rfof[i] = DUMMY_TIMESTEP_VAL;
            hnext[i] = DUMMY_TIMESTEP_VAL;
        }

/*************************** Physics **********************************/

/* Set time-step negative if tracing backwards in time */
#pragma omp simd
        for (size_t i = 0; i < vector_size; i++)
        {
            if (sim->options->reverse_time)
            {
                hin[i] = -hin[i];
            }
        }

        /* Cash-Karp method for orbit-following */
        if (sim->options->enable_orbit_following)
        {
            if (sim->options->enable_mhd)
            {
                step_gc_cashkarp_mhd(
                    &p, hin, hout_orb, tol_orb, &sim->bfield, &sim->efield,
                    sim->boozer, &sim->mhd, sim->options->enable_aldforce);
            }
            else
            {
                step_gc_cashkarp(
                    &p, hin, hout_orb, tol_orb, &sim->bfield, &sim->efield,
                    sim->options->enable_aldforce);
            }
/* Check whether time step was rejected */
#pragma omp simd
            for (size_t i = 0; i < vector_size; i++)
            {
                /* Switch sign of the time-step again if it was reverted earlier
                 */
                if (sim->options->reverse_time)
                {
                    hout_orb[i] = -hout_orb[i];
                    hin[i] = -hin[i];
                }
                if (p.running[i] && hout_orb[i] < 0)
                {
                    p.running[i] = 0;
                    hnext[i] = hout_orb[i];
                }
            }
        }

        /* Milstein method for collisions */
        if (sim->options->enable_coulomb_collisions)
        {
            random_normal_simd(sim->random_data, 5 * vector_size, rnd);
            mccc_gc_milstein(
                &p, hin, hout_col, tol_col, wienarr, &sim->bfield, &sim->plasma,
                sim->mccc_data, rnd);

/* Check whether time step was rejected */
#pragma omp simd
            for (size_t i = 0; i < vector_size; i++)
            {
                if (p.running[i] && hout_col[i] < 0)
                {
                    p.running[i] = 0;
                    hnext[i] = hout_col[i];
                }
            }
        }

        /* Performs the ICRH kick if in resonance. */
        if (sim->options->enable_icrh)
        {
            rfof_resonance_check_and_kick_gc(
                &p, hin, hout_rfof, &rfof_mrk, &sim->rfof, &sim->bfield);

/* Check whether time step was rejected */
#pragma omp simd
            for (size_t i = 0; i < vector_size; i++)
            {
                if (p.running[i] && hout_rfof[i] < 0)
                {
                    p.running[i] = 0;
                    hnext[i] = hout_rfof[i];
                }
            }
        }

        /**********************************************************************/

        cputime = A5_WTIME;
#pragma omp simd
        for (size_t i = 0; i < vector_size; i++)
        {
            if (p.id[i] > 0 && !p.err[i])
            {
                /* Check other time step limitations */
                if (hnext[i] > 0)
                {
                    real dphi = fabs(p0.phi[i] - p.phi[i]) /
                                sim->options->adaptive_max_dphi;
                    real drho = fabs(p0.rho[i] - p.rho[i]) /
                                sim->options->adaptive_max_drho;

                    if (dphi > 1 && dphi > drho)
                    {
                        hnext[i] = -hin[i] / dphi;
                    }
                    else if (drho > 1 && drho > dphi)
                    {
                        hnext[i] = -hin[i] / drho;
                    }
                }

                /* Retrieve marker states in case time step was rejected      */
                if (hnext[i] < 0)
                {
                    marker_copy_gc(&p0, i, &p, i);
                }
                if (p.running[i])
                {

                    /* Advance time (if time step was accepted) and determine
                       next time step */
                    if (hnext[i] < 0)
                    {
                        /* if hnext < 0, you screwed up and had to copy the
                        previous state. Therefore, let us use the suggestion
                        given by the integrator when retaking the failed step.*/
                        hin[i] = -hnext[i];
                    }
                    else
                    {
                        p.time[i] +=
                            (1.0 - 2.0 * (sim->options->reverse_time > 0)) *
                            hin[i];
                        p.mileage[i] += hin[i];
                        /* In case the time step was succesful, pick the
                        smallest recommended value for the next step */
                        if (hnext[i] > hout_orb[i])
                        {
                            /* Use time step suggested by the orbit-following
                               integrator */
                            hnext[i] = hout_orb[i];
                        }
                        if (hnext[i] > hout_col[i])
                        {
                            /* Use time step suggested by the collision
                               integrator */
                            hnext[i] = hout_col[i];
                        }
                        if (hnext[i] > hout_rfof[i])
                        {
                            /* Use time step suggested by RFOF */
                            hnext[i] = hout_rfof[i];
                        }
                        if (hnext[i] == 1.0)
                        {
                            /* Time step is unchanged (happens when no physics
                               are enabled) */
                            hnext[i] = hin[i];
                        }
                        hin[i] = hnext[i];
                        if (sim->options->enable_coulomb_collisions)
                        {
                            /* Clear wiener processes */
                            mccc_wiener_clean(&(wienarr[i]), p.time[i]);
                        }
                    }

                    p.cputime[i] += cputime - cputime_last;
                }
            }
        }
        cputime_last = cputime;

        /* Check possible end conditions */
        endcond_check_gc(&p, &p0, sim);

        /* Update diagnostics */
        diag_update_gc(sim->diagnostics, sim->options, &sim->bfield, &p, &p0);

        /* Update number of running particles */
        n_running = marker_cycle_gc(pq, &p, &sim->bfield, cycle);

/* Determine simulation time-step for new particles */
#pragma omp simd
        for (size_t i = 0; i < vector_size; i++)
        {
            if (cycle[i] > 0)
            {
                hin[i] = simulate_gc_adaptive_inidt(sim, &p, i);
                if (sim->options->enable_coulomb_collisions)
                {
                    /* Re-allocate array storing the Wiener processes */
                    mccc_wiener_initialize(&(wienarr[i]), p.time[i]);
                }
                if (sim->options->enable_icrh)
                {
                    /* Reset icrh (rfof) resonance memory matrix. */
                    rfof_clear_history(&rfof_mrk, i);
                }
            }
        }
    }

    /* All markers simulated! */
    free(hin);
    free(hout_orb);
    free(hout_col);
    free(hout_rfof);
    free(hnext);
    free(cycle);
    free(rnd);

    /* Deallocate rfof structs */
    if (sim->options->enable_icrh)
    {
        rfof_tear_down(&rfof_mrk);
    }
}

/**
 * @brief Calculates time step value
 *
 * The returned time step is either directly user-defined, 1/100th of collision
 * frequency or user-defined fraction of gyro-motion.
 *
 * @param sim pointer to simulation data struct
 * @param p SIMD array of markers
 * @param i index of marker for which time step is assessed
 *
 * @return Calculated time step
 */
real simulate_gc_adaptive_inidt(Simulation *sim, MarkerGuidingCenter *p, int i)
{
    /* Just use some large value if no physics are defined */
    real h = DUMMY_TIMESTEP_VAL;

    /* Value defined directly by user */
    if (sim->options->use_explicit_fixedstep)
    {
        h = sim->options->explicit_fixedstep;
    }
    else
    {
        /* Value calculated from gyrotime */
        if (sim->options->enable_orbit_following)
        {
            real Bnorm = math_normc(p->B_r[i], p->B_phi[i], p->B_z[i]);
            real gyrotime = CONST_2PI / phys_gyrofreq_ppar(
                                            p->mass[i], p->charge[i], p->mu[i],
                                            p->ppar[i], Bnorm);
            if (h > gyrotime)
            {
                h = gyrotime;
            }
        }

        /* Value calculated from collision frequency */
        if (sim->options->enable_coulomb_collisions)
        {
            real nu = 1;
            /*mccc_collfreq_gc(p, &sim->B_data, &sim->plasma_data,
                sim->coldata, &nu, i); */

            /* Only small angle collisions so divide this by 100 */
            real colltime = 1 / (100 * nu);
            if (h > colltime)
            {
                h = colltime;
            }
        }
    }
    return h;
}
