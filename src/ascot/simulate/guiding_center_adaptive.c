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

/**
 * A simple struct to keep track of when a marker crosses OMP twice in same
 * direction.
 */
typedef struct
{
    unsigned int crossed_once : 1;  /**< Flag for first crossing.             */
    unsigned int crossed_twice : 1; /**< Flag for second crossing.            */
    unsigned int first_ppar : 1;    /**< Direction of the first crossing.     */
} Crossing;

/**
 * Keeps track of acceleration factor and data necessary to update it after each
 * poloidal orbit.
 */
typedef struct
{
    real *acc;       /**< Acceleration factor.                                */
    real *orbittime; /**< Passed orbit time [s].                              */
    real *collfreq;  /**< Average collision frequency along the orbit [1/s].  */
    Crossing *cross; /**< Stores information about OMP crossings.             */
} Acceleration;

/**
 * Recalculate acceleration factor.
 *
 * The acceleration is updated when crossing OMP. During the first crossing,
 * the counter for orbit time is started. For the next crossing, we check if
 * OMP was crossed in the same direction as first. If not, the crossing is
 * ignored and for the third crossing we check the direction again. If this is
 * not in the same direction as the first, the counters are nullified,
 * acceleration is set to one, and the process is started again. This way we can
 * account both for passing and banana particles, and for the cases where
 * collisions have changed the orbit topology.
 *
 * When we have two suitable crossings, the acceleration factor is updated and
 * the counter for the orbit time and crossings are nullified.
 *
 * @param acc acceleration struct
 * @param sim simulation struct
 * @param p current marker
 * @param p0 previous marker
 */
void recalculate_acceleration(
    Acceleration *acc, Simulation *sim, MarkerGuidingCenter *p,
    MarkerGuidingCenter *p0);

/**
 * Allocates struct representing acceleration struct.
 *
 * Size used for memory allocation is NSIMD for CPU run and the total number
 * of particles for GPU.
 *
 * @param acceleration struct to allocate
 * @param vector_size the number of markers that the struct represents
 */
void acceleration_allocate(Acceleration *acceleration, size_t vector_size);

/**
 * Offload acceleration struct to GPU.
 *
 * @param acceleration pointer to the acceleration struct to be offloaded.
 * @param vector_size the number of markers that the struct represents
 */
void acceleration_offload(Acceleration *acceleration, size_t vector_size);

/** Dummy time step value [s]. */
#define DUMMY_TIMESTEP_VAL 1.0

/**
 * Replace markers in the simulation vector with new ones from the queue.
 *
 * A marker is replaced if it is no longer running or if it is a dummy marker.
 *
 * @param vector_size Number of markers in the simulation vector.
 * @param queue Marker queue.
 * @param p_current Current marker vector.
 * @param sim Simulation data.
 * @param time_step Time step for each marker.
 *        This is set to the initial step size for newly added markers.
 */
static size_t cycle_markers(
    size_t vector_size, MarkerQueue *queue, MarkerGuidingCenter *p_current,
    Simulation *sim, real time_step[vector_size])
{
    size_t start = 0;
    while (start < vector_size)
    {
        size_t next_in_queue;
        size_t idx = MarkerQueue_cycle(
            &next_in_queue, queue, vector_size, start, p_current->id,
            p_current->running);
        if (idx == vector_size)
            break;
        if (p_current->id[idx] != 0)
            MarkerGuidingCenter_to_queue(queue, p_current, idx, &sim->bfield);
        p_current->id[idx] = 0;
        p_current->running[idx] = 0;

        if (next_in_queue < queue->n)
        {
            if (MarkerGuidingCenter_from_queue(
                    p_current, queue, idx, next_in_queue, &sim->bfield))
            {
                p_current->id[idx] = 0;
                p_current->running[idx] = 0;
            }
            time_step[idx] = sim->options->timestep;
        }
        start = idx;
    }

    size_t n_running = 0;
#pragma omp simd reduction(+ : n_running)
    for (size_t i = 0; i < vector_size; i++)
        n_running += p_current->running[i];

    return n_running;
}

int simulate_gc_adaptive(Simulation *sim, MarkerQueue *pq, size_t vector_size)
{

    /* Wiener arrays needed for the adaptive time step */
    mccc_wienarr* wienarr = (mccc_wienarr*) malloc(vector_size*sizeof(mccc_wienarr));
    Acceleration acceleration;
    acceleration_allocate(&acceleration, vector_size);

    /* Current time step, suggestions for the next time step and next time
     * step                                                                */
    int err = 0;
    real *hin = (real *)xmalloc(&err, vector_size * sizeof(real));
    if (err)
        return 1;
    real *hout_orb = (real *)xmalloc(&err, vector_size * sizeof(real));
    if (err)
        return 1;
    real *hout_col = (real *)xmalloc(&err, vector_size * sizeof(real));
    if (err)
        return 1;
    real *hout_rfof = (real *)xmalloc(&err, vector_size * sizeof(real));
    if (err)
        return 1;
    real *hnext = (real *)xmalloc(&err, vector_size * sizeof(real));
    if (err)
        return 1;
    real *rnd = (real *)xmalloc(&err, 5 * vector_size * sizeof(real));
    if (err)
        return 1;

    /* Flag indicating whether a new marker was initialized */
    size_t *cycle = (size_t *)xmalloc(&err, vector_size * sizeof(size_t));
    if (err)
        return 1;

    real tol_col = sim->options->adaptive_tolerance_collisions;
    real tol_orb = sim->options->adaptive_tolerance_orbit;

    real cputime, cputime_last; // Global cpu time: recent and previous record

    MarkerGuidingCenter p;  // This array holds current states
    MarkerGuidingCenter p0; // This array stores previous states
    if (MarkerGuidingCenter_allocate(&p, vector_size))
        return 1;
    if (MarkerGuidingCenter_allocate(&p0, vector_size))
        return 1;

    rfof_marker rfof_mrk;

    for (size_t i = 0; i < vector_size; i++)
    {
        p.id[i] = 0;
        p.running[i] = 0;
        acceleration.acc[i] = 1.0;
        acceleration.orbittime[i] = -1;
        acceleration.cross[i].crossed_once = 0;
    }

    /* Initialize running particles */
    size_t n_running = cycle_markers(vector_size, pq, &p, sim, hin);

    if (sim->options->enable_icrh)
        rfof_set_up(&rfof_mrk, sim->rfof);

    GPU_PARALLEL_LOOP_ALL_LEVELS
    for (size_t i = 0; i < vector_size; i++)
    {
        if (cycle[i] > 0)
        {
            hin[i] = sim->options->timestep;
            if (sim->options->enable_coulomb_collisions)
                mccc_wiener_initialize(&(wienarr[i]), p.time[i]);
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
    MarkerGuidingCenter_offload(&p);
    MarkerGuidingCenter_offload(&p0);
    acceleration_offload(&acceleration, vector_size);
    GPU_MAP_TO_DEVICE(
        hin [0:vector_size], rnd [0:5 * vector_size], hout_orb [0:vector_size],
        hout_col [0:vector_size], hout_rfof [0:vector_size],
        hnext [0:vector_size], cycle [0:vector_size])
    mccc_wiener_offload(wienarr, vector_size);
    while (n_running > 0)
    {

        /* Store marker states in case time step will be rejected */
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for (size_t i = 0; i < vector_size; i++)
        {
            MarkerGuidingCenter_copy(&p0, &p, i);
            hout_orb[i] = DUMMY_TIMESTEP_VAL;
            hout_col[i] = DUMMY_TIMESTEP_VAL;
            hout_rfof[i] = DUMMY_TIMESTEP_VAL;
            hnext[i] = DUMMY_TIMESTEP_VAL;
        }

        /*************************** Physics **********************************/

        if (sim->options->enable_orbit_following)
        {

            GPU_PARALLEL_LOOP_ALL_LEVELS
            for (size_t i = 0; i < vector_size; i++)
                hin[i] = (1 - 2*(sim->options->reverse_time)) * hin[i];

            if (sim->options->enable_mhd)
                step_gc_cashkarp_mhd(
                    &p, hin, hout_orb, tol_orb, &sim->bfield, &sim->efield,
                    sim->boozer, &sim->mhd, sim->options->enable_aldforce);
            else
                step_gc_cashkarp(
                    &p, hin, hout_orb, tol_orb, &sim->bfield, &sim->efield,
                    sim->options->enable_aldforce);

            /* Check whether time step was rejected */
            GPU_PARALLEL_LOOP_ALL_LEVELS
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
            random_normal_simd(sim->random_data, 5 * p.size, rnd);
            mccc_gc_milstein(
                &p, hin, acceleration.acc, acceleration.collfreq, hout_col,
                tol_col, wienarr, &sim->bfield, &sim->plasma, sim->mccc_data,
                rnd);

            /* Check whether time step was rejected */
            GPU_PARALLEL_LOOP_ALL_LEVELS
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
                &p, hin, hout_rfof, &rfof_mrk, sim->rfof, &sim->bfield);

            /* Check whether time step was rejected */
            GPU_PARALLEL_LOOP_ALL_LEVELS
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
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for (size_t i = 0; i < vector_size; i++)
        {
            if (p.id[i] > 0 && !p.err[i])
            {

                /* Retrieve marker states in case time step was rejected */
                if (hnext[i] < 0)
                {
                    MarkerGuidingCenter_copy(&p, &p0, i);
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
                            hin[i] * acceleration.acc[i];
                        p.mileage[i] += hin[i] * acceleration.acc[i];
                        if (acceleration.orbittime[i] >= 0)
                            acceleration.orbittime[i] +=
                                hin[i] * acceleration.acc[i];
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
        if (sim->options->enable_adaptive > 1)
        {
            recalculate_acceleration(&acceleration, sim, &p, &p0);
        }
        cputime_last = cputime;
        endcond_check_gc(&p, &p0, sim);
        Diag_update_gc(&sim->diagnostics, &sim->bfield, &p, &p0);
#ifdef GPU
        n_running = 0;
        GPU_PARALLEL_LOOP_ALL_LEVELS_REDUCTION(n_running)
        for (size_t i = 0; i < p.size; i++)
        {
            if (p.running[i] > 0)
                n_running++;
        }
#else
        n_running = cycle_markers(vector_size, pq, &p, sim, hin);
#endif
        /* Determine simulation time-step for new particles */
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for (size_t i = 0; i < vector_size; i++)
        {
            if (cycle[i] > 0)
            {
                acceleration.acc[i] = 1.0;
                acceleration.orbittime[i] = -1;
                acceleration.cross[i].crossed_once = 0;
                hin[i] = sim->options->timestep;
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
    MarkerGuidingCenter_onload(&p);
    MarkerGuidingCenter_onload(&p0);
#ifdef GPU
    GPU_MAP_FROM_DEVICE(sim [0:1])
    n_running = cycle_markers(vector_size, pq, &p, sim, hin);
#endif

    /* All markers simulated! */
    free(hin);
    free(hout_orb);
    free(hout_col);
    free(hout_rfof);
    free(hnext);
    free(cycle);
    free(rnd);

    if (sim->options->enable_icrh)
        rfof_tear_down(&rfof_mrk);

    MarkerGuidingCenter_deallocate(&p0);
    MarkerGuidingCenter_deallocate(&p);
    return 0;
}

void recalculate_acceleration(
    Acceleration *acc, Simulation *sim, MarkerGuidingCenter *p,
    MarkerGuidingCenter *p0)
{
    real rz[2];
    real SAFETY_FACTOR = (float)sim->options->enable_adaptive / 1000.0;
    GPU_PARALLEL_LOOP_ALL_LEVELS
    for (size_t i = 0; i < p->size; i++)
    {
        Bfield_eval_axis_rz(rz, &sim->bfield, p->phi[i]);
        int omp_crossed =
            ((p->z[i] - rz[1]) * (p0->z[i] - rz[1]) < 0) && p->r[i] > rz[0];
        if (omp_crossed && acc->cross[i].crossed_twice)
        {
            if (((float)acc->cross[i].first_ppar - 0.5) * p->ppar[i] > 0)
            {
                acc->acc[i] = fmax(
                    1.0,
                    SAFETY_FACTOR / (acc->orbittime[i] * acc->collfreq[i]));
            }
            else
            {
                acc->acc[i] = 1;
            }
            acc->cross[i].crossed_once = 1;
            acc->cross[i].crossed_twice = 0;
            acc->cross[i].first_ppar = p->ppar[i] > 0;
            acc->orbittime[i] = 0;
        }
        else if (omp_crossed && acc->cross[i].crossed_once)
        {
            acc->cross[i].crossed_twice = 1;
            if (((float)acc->cross[i].first_ppar - 0.5) * p->ppar[i] > 0)
            {
                acc->acc[i] = fmax(
                    1.0,
                    SAFETY_FACTOR / (acc->orbittime[i] * acc->collfreq[i]));
                acc->cross[i].crossed_once = 1;
                acc->cross[i].crossed_twice = 0;
                acc->cross[i].first_ppar = p->ppar[i] > 0;
                acc->orbittime[i] = 0;
            }
        }
        else if (omp_crossed)
        {
            acc->cross[i].crossed_once = 1;
            acc->cross[i].first_ppar = p->ppar[i] > 0;
            acc->orbittime[i] = 0;
        }
    }
}

void acceleration_allocate(Acceleration *acceleration, size_t vector_size)
{
    acceleration->acc = malloc(vector_size * sizeof(acceleration->acc));
    acceleration->orbittime =
        malloc(vector_size * sizeof(acceleration->orbittime));
    acceleration->collfreq =
        malloc(vector_size * sizeof(acceleration->collfreq));
    acceleration->cross = malloc(vector_size * sizeof(acceleration->cross));
}

void acceleration_offload(Acceleration *acceleration, size_t vector_size)
{
    SUPPRESS_UNUSED_WARNING(acceleration);
    SUPPRESS_UNUSED_WARNING(vector_size);
    GPU_MAP_TO_DEVICE(
        acceleration [0:1], acceleration->acc [0:vector_size],
        acceleration->orbittime [0:vector_size],
        acceleration->collfreq [0:vector_size],
        acceleration->cross [0:vector_size])
}
