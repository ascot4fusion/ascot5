/**
 * Implements marker.h.
 */
#include "bfield.h"
#include "consts.h"
#include "defines.h"
#include "efield.h"
#include "utils/gctransform.h"
#include "marker.h"
#include "utils/mathlib.h"
#include "utils/physlib.h"
#include "datatypes.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * @brief Allocates struct representing particle markers
 *
 * Size used for memory allocation is NSIMD for CPU run and the total number
 * of particles for GPU.
 *
 * @param p_go marker struct to allocate
 * @param nmrk the number of markers that the struct represents
 */
void marker_allocate_go(MarkerGyroOrbit *mrk, size_t nmrk)
{
    mrk->r = malloc(nmrk * sizeof(mrk->r));
    mrk->phi = malloc(nmrk * sizeof(mrk->phi));
    mrk->z = malloc(nmrk * sizeof(mrk->z));
    mrk->p_r = malloc(nmrk * sizeof(mrk->p_r));
    mrk->p_phi = malloc(nmrk * sizeof(mrk->p_phi));
    mrk->p_z = malloc(nmrk * sizeof(mrk->p_z));
    mrk->mass = malloc(nmrk * sizeof(mrk->mass));
    mrk->charge = malloc(nmrk * sizeof(mrk->charge));
    mrk->time = malloc(nmrk * sizeof(mrk->time));
    mrk->znum = malloc(nmrk * sizeof(mrk->znum));
    mrk->anum = malloc(nmrk * sizeof(mrk->anum));
    mrk->B_r = malloc(nmrk * sizeof(mrk->B_r));
    mrk->B_phi = malloc(nmrk * sizeof(mrk->B_phi));
    mrk->B_z = malloc(nmrk * sizeof(mrk->B_z));
    mrk->B_r_dr = malloc(nmrk * sizeof(mrk->B_r_dr));
    mrk->B_phi_dr = malloc(nmrk * sizeof(mrk->B_phi_dr));
    mrk->B_z_dr = malloc(nmrk * sizeof(mrk->B_z_dr));
    mrk->B_r_dphi = malloc(nmrk * sizeof(mrk->B_r_dphi));
    mrk->B_phi_dphi = malloc(nmrk * sizeof(mrk->B_phi_dphi));
    mrk->B_z_dphi = malloc(nmrk * sizeof(mrk->B_z_dphi));
    mrk->B_r_dz = malloc(nmrk * sizeof(mrk->B_r_dz));
    mrk->B_phi_dz = malloc(nmrk * sizeof(mrk->B_phi_dz));
    mrk->B_z_dz = malloc(nmrk * sizeof(mrk->B_z_dz));
    mrk->bounces = malloc(nmrk * sizeof(mrk->bounces));
    mrk->weight = malloc(nmrk * sizeof(mrk->weight));
    mrk->cputime = malloc(nmrk * sizeof(mrk->cputime));
    mrk->rho = malloc(nmrk * sizeof(mrk->rho));
    mrk->theta = malloc(nmrk * sizeof(mrk->theta));
    mrk->id = malloc(nmrk * sizeof(mrk->id));
    mrk->endcond = malloc(nmrk * sizeof(mrk->endcond));
    mrk->walltile = malloc(nmrk * sizeof(mrk->walltile));
    mrk->mileage = malloc(nmrk * sizeof(mrk->mileage));
    mrk->running = malloc(nmrk * sizeof(mrk->running));
    mrk->err = malloc(nmrk * sizeof(mrk->err));
    mrk->index = malloc(nmrk * sizeof(mrk->index));
    mrk->n_mrk = nmrk;
}

/**
 * @brief Allocates struct representing particle markers
 *
 * Size used for memory allocation is NSIMD for CPU run and the total number
 * of particles for GPU.
 *
 * @param mrk marker struct to allocate
 * @param nmrk the number of markers that the struct represents
 */
void marker_allocate_gc(MarkerGuidingCenter *mrk, size_t nmrk)
{
    mrk->r = malloc(nmrk * sizeof(mrk->r));
    mrk->phi = malloc(nmrk * sizeof(mrk->phi));
    mrk->z = malloc(nmrk * sizeof(mrk->z));
    mrk->ppar = malloc(nmrk * sizeof(mrk->ppar));
    mrk->mu = malloc(nmrk * sizeof(mrk->mu));
    mrk->zeta = malloc(nmrk * sizeof(mrk->zeta));
    mrk->mass = malloc(nmrk * sizeof(mrk->mass));
    mrk->charge = malloc(nmrk * sizeof(mrk->charge));
    mrk->time = malloc(nmrk * sizeof(mrk->time));
    mrk->znum = malloc(nmrk * sizeof(mrk->znum));
    mrk->anum = malloc(nmrk * sizeof(mrk->anum));
    mrk->B_r = malloc(nmrk * sizeof(mrk->B_r));
    mrk->B_phi = malloc(nmrk * sizeof(mrk->B_phi));
    mrk->B_z = malloc(nmrk * sizeof(mrk->B_z));
    mrk->B_r_dr = malloc(nmrk * sizeof(mrk->B_r_dr));
    mrk->B_phi_dr = malloc(nmrk * sizeof(mrk->B_phi_dr));
    mrk->B_z_dr = malloc(nmrk * sizeof(mrk->B_z_dr));
    mrk->B_r_dphi = malloc(nmrk * sizeof(mrk->B_r_dphi));
    mrk->B_phi_dphi = malloc(nmrk * sizeof(mrk->B_phi_dphi));
    mrk->B_z_dphi = malloc(nmrk * sizeof(mrk->B_z_dphi));
    mrk->B_r_dz = malloc(nmrk * sizeof(mrk->B_r_dz));
    mrk->B_phi_dz = malloc(nmrk * sizeof(mrk->B_phi_dz));
    mrk->B_z_dz = malloc(nmrk * sizeof(mrk->B_z_dz));
    mrk->bounces = malloc(nmrk * sizeof(mrk->bounces));
    mrk->weight = malloc(nmrk * sizeof(mrk->weight));
    mrk->cputime = malloc(nmrk * sizeof(mrk->cputime));
    mrk->rho = malloc(nmrk * sizeof(mrk->rho));
    mrk->theta = malloc(nmrk * sizeof(mrk->theta));
    mrk->id = malloc(nmrk * sizeof(mrk->id));
    mrk->endcond = malloc(nmrk * sizeof(mrk->endcond));
    mrk->walltile = malloc(nmrk * sizeof(mrk->walltile));
    mrk->mileage = malloc(nmrk * sizeof(mrk->mileage));
    mrk->running = malloc(nmrk * sizeof(mrk->running));
    mrk->err = malloc(nmrk * sizeof(mrk->err));
    mrk->index = malloc(nmrk * sizeof(mrk->index));
    mrk->n_mrk = nmrk;
}

/**
 * @brief Allocates struct representing particle markers
 *
 * Size used for memory allocation is NSIMD for CPU run and the total number
 * of particles for GPU.
 *
 * @param mrk marker struct to allocate
 * @param nmrk the number of markers that the struct represents
 */
void marker_allocate_fl(MarkerFieldLine *mrk, size_t nmrk)
{
    mrk->r = malloc(nmrk * sizeof(mrk->r));
    mrk->phi = malloc(nmrk * sizeof(mrk->phi));
    mrk->z = malloc(nmrk * sizeof(mrk->z));
    mrk->pitch = malloc(nmrk * sizeof(mrk->pitch));
    mrk->B_r = malloc(nmrk * sizeof(mrk->B_r));
    mrk->B_phi = malloc(nmrk * sizeof(mrk->B_phi));
    mrk->B_z = malloc(nmrk * sizeof(mrk->B_z));
    mrk->B_r_dr = malloc(nmrk * sizeof(mrk->B_r_dr));
    mrk->B_phi_dr = malloc(nmrk * sizeof(mrk->B_phi_dr));
    mrk->B_z_dr = malloc(nmrk * sizeof(mrk->B_z_dr));
    mrk->B_r_dphi = malloc(nmrk * sizeof(mrk->B_r_dphi));
    mrk->B_phi_dphi = malloc(nmrk * sizeof(mrk->B_phi_dphi));
    mrk->B_z_dphi = malloc(nmrk * sizeof(mrk->B_z_dphi));
    mrk->B_r_dz = malloc(nmrk * sizeof(mrk->B_r_dz));
    mrk->B_phi_dz = malloc(nmrk * sizeof(mrk->B_phi_dz));
    mrk->B_z_dz = malloc(nmrk * sizeof(mrk->B_z_dz));
    mrk->bounces = malloc(nmrk * sizeof(mrk->bounces));
    mrk->cputime = malloc(nmrk * sizeof(mrk->cputime));
    mrk->rho = malloc(nmrk * sizeof(mrk->rho));
    mrk->theta = malloc(nmrk * sizeof(mrk->theta));
    mrk->id = malloc(nmrk * sizeof(mrk->id));
    mrk->endcond = malloc(nmrk * sizeof(mrk->endcond));
    mrk->walltile = malloc(nmrk * sizeof(mrk->walltile));
    mrk->mileage = malloc(nmrk * sizeof(mrk->mileage));
    mrk->running = malloc(nmrk * sizeof(mrk->running));
    mrk->err = malloc(nmrk * sizeof(mrk->err));
    mrk->index = malloc(nmrk * sizeof(mrk->index));
    mrk->n_mrk = nmrk;
}

/**
 * @brief Makes a dummy FO simulation marker
 *
 * A dummy marker is a marker whose fields have possible but unrealistic
 * values. It is intended to act as a filler for SIMD arrays that would
 * otherwise be left uninitialized.
 *
 * A dummy marker has negative one (-1) ID and negative one (-1) queue index,
 * and it is not running (0). Dummy marker should be ignored when encountered
 * during the simulation.
 *
 * @param p_go  pointer to MarkerGyroOrbit array where dummy is placed
 * @param j     index of the new dummy in the SIMD arrays
 */
void marker_to_go_dummy(MarkerGyroOrbit *p_go, size_t j)
{
    p_go->r[j] = 1;
    p_go->phi[j] = 1;
    p_go->z[j] = 1;
    p_go->p_r[j] = 1;
    p_go->p_phi[j] = 1;
    p_go->p_z[j] = 1;
    p_go->mass[j] = 1;
    p_go->charge[j] = 1;
    p_go->znum[j] = 1;
    p_go->anum[j] = 1;
    p_go->bounces[j] = 0;
    p_go->weight[j] = 0;
    p_go->time[j] = 0;
    p_go->id[j] = 0;
    p_go->mileage[j] = 0;
    p_go->running[j] = 0;
    p_go->endcond[j] = 0;
    p_go->walltile[j] = 0;
    p_go->B_r[j] = 1;
    p_go->B_phi[j] = 1;
    p_go->B_z[j] = 1;
    p_go->index[j] = -1;
    p_go->err[j] = 0;
}

/**
 * @brief Makes a dummy GC simulation marker
 *
 * A dummy marker is a marker whose fields have possible but unrealistic
 * values. It is intended to act as a filler for SIMD arrays that would
 * otherwise be left uninitialized.
 *
 * A dummy marker has negative one (-1) ID and negative one (-1) queue index,
 * and it is not running (0). Dummy marker should be ignored when encountered
 * during the simulation.
 *
 * @param p_gc  pointer to MarkerGuidingCenter array where dummy is placed
 * @param j     index of the new dummy in the SIMD arrays
 */
void marker_to_gc_dummy(MarkerGuidingCenter *p_gc, size_t j)
{
    p_gc->r[j] = 1;
    p_gc->phi[j] = 1;
    p_gc->z[j] = 1;
    p_gc->ppar[j] = 1;
    p_gc->mu[j] = 1;
    p_gc->zeta[j] = 1;
    p_gc->mass[j] = 1;
    p_gc->charge[j] = 1;
    p_gc->time[j] = 0;
    p_gc->bounces[j] = 0;
    p_gc->weight[j] = 0;
    p_gc->id[j] = 0;
    p_gc->mileage[j] = 0;
    p_gc->B_r[j] = 1;
    p_gc->B_r_dr[j] = 1;
    p_gc->B_r_dphi[j] = 1;
    p_gc->B_r_dz[j] = 1;

    p_gc->B_phi[j] = 1;
    p_gc->B_phi_dr[j] = 1;
    p_gc->B_phi_dphi[j] = 1;
    p_gc->B_phi_dz[j] = 1;

    p_gc->B_z[j] = 1;
    p_gc->B_z_dr[j] = 1;
    p_gc->B_z_dphi[j] = 1;
    p_gc->B_z_dz[j] = 1;
    p_gc->index[j] = -1;
    p_gc->err[j] = 0;
}

/**
 * @brief Makes a dummy ML simulation marker
 *
 * A dummy marker is a marker whose fields have possible but unrealistic
 * values. It is intended to act as a filler for SIMD arrays that would
 * otherwise be left uninitialized.
 *
 * A dummy marker has negative one (-1) ID and negative one (-1) queue index,
 * and it is not running (0). Dummy marker should be ignored when encountered
 * during the simulation.
 *
 * @param p_fl  pointer to MarkerFieldLine array where dummy is placed
 * @param j     index of the new dummy in the SIMD arrays
 */
void marker_to_fl_dummy(MarkerFieldLine *p_fl, size_t j)
{
    p_fl->r[j] = 1;
    p_fl->phi[j] = 1;
    p_fl->z[j] = 1;
    p_fl->time[j] = 0;
    p_fl->id[j] = 0;
    p_fl->mileage[j] = 0;
    p_fl->B_r[j] = 1;
    p_fl->B_r_dr[j] = 1;
    p_fl->B_r_dphi[j] = 1;
    p_fl->B_r_dz[j] = 1;

    p_fl->B_phi[j] = 1;
    p_fl->B_phi_dr[j] = 1;
    p_fl->B_phi_dphi[j] = 1;
    p_fl->B_phi_dz[j] = 1;

    p_fl->B_z[j] = 1;
    p_fl->B_z_dr[j] = 1;
    p_fl->B_z_dphi[j] = 1;
    p_fl->B_z_dz[j] = 1;
    p_fl->index[j] = -1;
    p_fl->err[j] = 0;
}

/**
 * @brief Replace finished FO markers with new ones or dummies
 *
 * A marker has finished simulation when marker.running = 0. If queue has
 * unsimulated markers, one is picked to replace the finished marker. If not,
 * a dummy marker is used as a replacement instead. Finished marker is converted
 * to marker state and stored in the queue.
 *
 * This function updates queue.next and queue.finished fields when a marker has
 * finished simulation. This is done thread-safe.
 *
 * This function returns values indicating what was done for each marker in a
 * SIMD array:
 *   0 : Nothing
 *  -1 : Finished marker replaced with a dummy (queue is empty)
 *   1 : Finished marker replaced with a fresh one
 *
 * @param q pointer to marker queue
 * @param p pointer to SIMD structure of markers
 * @param bfield pointer to magnetic field data
 * @param cycle pointer to integer array where what was done for each marker
 *              is stored
 *
 * @return Number of markers within the SIMD structure that are still running
 */
size_t marker_cycle_go(
    MarkerQueue *q, MarkerGyroOrbit *p, Bfield *bfield, size_t *cycle)
{

    /* Loop over markers.
     * A SIMD loop is not possible as we modify the queue. */
    for (size_t i = 0; i < p->n_mrk; i++)
    {

        /* First check whether we should pick a new marker */
        int newmarker = 0;

        /* 1. There are markers in queue and this marker is dummy */
        if (p->id[i] == 0 && q->next < q->n)
        {
            newmarker = 1;
        }

        /* 2. This marker has finished simulation */
        if (!p->running[i] && p->id[i] > 0)
        {
            /* Convert finished marker to state
             * and store it back to the queue */
            marker_go_to_state(p, i, q->p[p->index[i]], bfield);
            newmarker = 1;
#pragma omp critical
            q->finished++;
        }

        /* Init a new marker if one is needed */
        cycle[i] = 0;
        while (newmarker)
        {
            /* Get the next unsimulated marker from the queue */
            size_t i_prt;
#pragma omp critical
            i_prt = q->next++;

            if (i_prt >= q->n)
            {
                /* The queue is empty, place a dummy marker here */
                p->running[i] = 0;
                p->id[i] = -1;
                cycle[i] = -1;
                newmarker = 0;
            }
            else if (q->p[i_prt]->endcond)
            {
/* This marker already has an active end condition. Try next. */
#pragma omp critical
                q->finished++;
            }
            else
            {
                /* Try to convert marker state to a simulation marker */
                err_t err = state_to_go(q->p[i_prt], i_prt, p, i, bfield);
                if (err)
                {
                    /* Failed! Mark the marker candidate as finished
                     * and try with a new marker state*/
                    q->p[i_prt]->err = err;
#pragma omp critical
                    q->finished++;
                }
                else
                {
                    /* Success! We are good to go. */
                    cycle[i] = 1;
                    newmarker = 0;
                }
            }
        }
    }

    size_t n_running = 0;
#pragma omp simd reduction(+ : n_running)
    for (size_t i = 0; i < p->n_mrk; i++)
    {
        n_running += p->running[i];
    }

    return n_running;
}

/**
 * @brief Replace finished GC markers with new ones or dummies
 *
 * A marker has finished simulation when marker.running = 0. If queue has
 * unsimulated markers, one is picked to replace the finished marker. If not,
 * a dummy marker is used as a replacement instead. Finished marker is converted
 * to marker state and stored in the queue.
 *
 * This function updates queue.next and queue.finished fields when a marker has
 * finished simulation. This is done thread-safe.
 *
 * This function returns values indicating what was done for each marker in a
 * SIMD array:
 *   0 : Nothing
 *  -1 : Finished marker replaced with a dummy (queue is empty)
 *   1 : Finished marker replaced with a fresh one
 *
 * @param q pointer to marker queue
 * @param p pointer to SIMD structure of markers
 * @param bfield pointer to magnetic field data
 * @param cycle pointer to integer array where what was done for each marker
 *              is stored
 *
 * @return Number of markers within the SIMD structure that are still running
 */
size_t marker_cycle_gc(
    MarkerQueue *q, MarkerGuidingCenter *p, Bfield *bfield, size_t *cycle)
{

    /* Loop over markers.
     * A SIMD loop is not possible as we modify the queue. */
    for (size_t i = 0; i < NSIMD; i++)
    {

        /* First check whether we should pick a new marker */
        int newmarker = 0;

        /* 1. There are markers in queue and this marker is dummy */
        if (p->id[i] == 0 && q->next < q->n)
        {
            newmarker = 1;
        }

        /* 2. This marker has finished simulation */
        if (!p->running[i] && p->id[i] > 0)
        {
            /* Convert finished marker to state
             * and store it back to the queue */
            marker_gc_to_state(p, i, q->p[p->index[i]], bfield);
            newmarker = 1;
#pragma omp critical
            q->finished++;
        }

        /* Init a new marker if one is needed */
        cycle[i] = 0;
        while (newmarker)
        {
            /* Get the next unsimulated marker from the queue */
            size_t i_prt;
#pragma omp critical
            i_prt = q->next++;

            if (i_prt >= q->n)
            {
                /* The queue is empty, place a dummy marker here */
                p->running[i] = 0;
                p->id[i] = -1;
                cycle[i] = -1;
                newmarker = 0;
            }
            else
            {
                /* Try to convert marker state to a simulation marker */
                err_t err = state_to_gc(q->p[i_prt], i_prt, p, i, bfield);
                if (err)
                {
                    /* Failed! Mark the marker candidate as finished
                     * and try with a new marker state*/
                    q->p[i_prt]->err = err;
#pragma omp critical
                    q->finished++;
                }
                else
                {
                    /* Success! We are good to go. */
                    cycle[i] = 1;
                    newmarker = 0;
                }
            }
        }
    }

    size_t n_running = 0;
#pragma omp simd reduction(+ : n_running)
    for (size_t i = 0; i < NSIMD; i++)
    {
        n_running += p->running[i];
    }

    return n_running;
}

/**
 * @brief Replace finished ML markers with new ones or dummies
 *
 * A marker has finished simulation when marker.running = 0. If queue has
 * unsimulated markers, one is picked to replace the finished marker. If not,
 * a dummy marker is used as a replacement instead. Finished marker is converted
 * to marker state and stored in the queue.
 *
 * This function updates queue.next and queue.finished fields when a marker has
 * finished simulation. This is done thread-safe.
 *
 * This function returns values indicating what was done for each marker in a
 * SIMD array:
 *   0 : Nothing
 *  -1 : Finished marker replaced with a dummy (queue is empty)
 *   1 : Finished marker replaced with a fresh one
 *
 * @param q pointer to particle queue
 * @param p pointer to SIMD structure of markers
 * @param bfield pointer to magnetic field data
 * @param cycle pointer to integer array where what was done for each marker
 *              is stored
 *
 * @return Number of markers within the SIMD structure that are still running
 */
size_t marker_cycle_fl(
    MarkerQueue *q, MarkerFieldLine *p, Bfield *bfield, size_t *cycle)
{

    /* Loop over markers.
     * A SIMD loop is not possible as we modify the queue. */
    for (size_t i = 0; i < NSIMD; i++)
    {

        /* First check whether we should pick a new marker */
        int newmarker = 0;

        /* 1. There are markers in queue and this marker is dummy */
        if (p->id[i] == 0 && q->next < q->n)
        {
            newmarker = 1;
        }

        /* 2. This marker has finished simulation */
        if (!p->running[i] && p->id[i] > 0)
        {
            /* Convert finished marker to state
             * and store it back to the queue */
            marker_fl_to_state(p, i, q->p[p->index[i]]);
            newmarker = 1;
#pragma omp critical
            q->finished++;
        }

        /* Init a new marker if one is needed */
        cycle[i] = 0;
        while (newmarker)
        {
            /* Get the next unsimulated marker from the queue */
            size_t i_prt;
#pragma omp critical
            i_prt = q->next++;

            if (i_prt >= q->n)
            {
                /* The queue is empty, place a dummy marker here */
                p->running[i] = 0;
                p->id[i] = -1;
                cycle[i] = -1;
                newmarker = 0;
            }
            else
            {
                /* Try to convert marker state to a simulation marker */
                err_t err = state_to_fl(q->p[i_prt], i_prt, p, i, bfield);
                if (err)
                {
                    /* Failed! Mark the marker candidate as finished
                     * and try with a new marker state*/
                    q->p[i_prt]->err = err;
#pragma omp critical
                    q->finished++;
                }
                else
                {
                    /* Success! We are good to go. */
                    cycle[i] = 1;
                    newmarker = 0;
                }
            }
        }
    }

    size_t n_running = 0;
#pragma omp simd reduction(+ : n_running)
    for (size_t i = 0; i < NSIMD; i++)
    {
        n_running += p->running[i];
    }

    return n_running;
}

/**
 * @brief Convert state into a FO SIMD marker
 *
 * This function assumes markers are drawn from a marker queue where they will
 * be returned at the same position once simulation for this marker has ended,
 *
 * State is converted into a FO marker by copying all the fields, except for
 * magnetic field (and rho) which is evaluated here at FO position and then
 * stored.
 *
 * The simulation marker is set as running unless there is an error or the state
 * already have active end condition.
 *
 * If state marker already has an error flag, this error flag is returned and
 * nothing is done. If an error occurs within this function, the conversion
 * is terminated and the error is returned. It is the responsibility of the
 * caller to store the error in correct struct.
 *
 * This is a SIMD function.
 *
 * @param p pointer to a state being converted
 * @param i index of this state in the marker queue
 * @param p_go pointer to SIMD structure where marker is being stored
 * @param j index where in the SIMD structure marker is stored
 * @param bfield pointer to magnetic field data
 */
err_t state_to_go(State *p, size_t i, MarkerGyroOrbit *p_go, size_t j, Bfield *bfield)
{
    err_t err = p->err;

    if (!err)
    {
        p_go->r[j] = p->rprt;
        p_go->phi[j] = p->phiprt;
        p_go->z[j] = p->zprt;
        p_go->p_r[j] = p->pr;
        p_go->p_phi[j] = p->pphi;
        p_go->p_z[j] = p->pz;

        p_go->mass[j] = p->mass;
        p_go->charge[j] = p->charge;
        p_go->znum[j] = p->znum;
        p_go->anum[j] = p->anum;
        p_go->bounces[j] = 0;
        p_go->weight[j] = p->weight;
        p_go->time[j] = p->time;
        p_go->theta[j] = p->theta;
        p_go->id[j] = p->id;
        p_go->endcond[j] = p->endcond;
        p_go->walltile[j] = p->walltile;
        p_go->mileage[j] = p->mileage;
    }

    /* Magnetic field stored in state is for the gc position */
    real B_dB[15], psi[1], rho[2];
    if (!err)
    {
        err = Bfield_eval_b_db(
            B_dB, p->rprt, p->phiprt, p->zprt, p->time, bfield);
    }
    if (!err)
    {
        err =
            Bfield_eval_psi(psi, p->rprt, p->phiprt, p->zprt, p->time, bfield);
    }
    if (!err)
    {
        err = Bfield_eval_rho(rho, psi[0], bfield);
    }

    if (!err)
    {
        p_go->rho[j] = rho[0];

        p_go->B_r[j] = B_dB[0];
        p_go->B_r_dr[j] = B_dB[1];
        p_go->B_r_dphi[j] = B_dB[2];
        p_go->B_r_dz[j] = B_dB[3];

        p_go->B_phi[j] = B_dB[4];
        p_go->B_phi_dr[j] = B_dB[5];
        p_go->B_phi_dphi[j] = B_dB[6];
        p_go->B_phi_dz[j] = B_dB[7];

        p_go->B_z[j] = B_dB[8];
        p_go->B_z_dr[j] = B_dB[9];
        p_go->B_z_dphi[j] = B_dB[10];
        p_go->B_z_dz[j] = B_dB[11];

        p_go->running[j] = 1;
        if (p->endcond)
        {
            p_go->running[j] = 0;
        }
        p_go->cputime[j] = p->cputime;
        p_go->index[j] = i;

        p_go->err[j] = 0;
    }

    return err;
}

/**
 * @brief Convert FO to state
 *
 * This function converts FO simulation marker to a marker state. No new state
 * structure is initialized, but an existing one is filled with parameters
 * corresponding to the converted marker.
 *
 * This is a SIMD function.
 *
 * @param p_go pointer to SIMD structure
 * @param j SIMD index of the marker to be converted
 * @param p pointer to state whose fields are used to store the converted marker
 * @param bfield pointer to magnetic field data
 */
void marker_go_to_state(
    MarkerGyroOrbit *p_go, size_t j, State *p, Bfield *bfield)
{
    err_t err = 0;

    p->rprt = p_go->r[j];
    p->phiprt = p_go->phi[j];
    p->zprt = p_go->z[j];
    p->pr = p_go->p_r[j];
    p->pphi = p_go->p_phi[j];
    p->pz = p_go->p_z[j];

    p->mass = p_go->mass[j];
    p->charge = p_go->charge[j];
    p->znum = p_go->znum[j];
    p->anum = p_go->anum[j];
    p->weight = p_go->weight[j];
    p->time = p_go->time[j];
    p->theta = p_go->theta[j];
    p->id = p_go->id[j];
    p->endcond = p_go->endcond[j];
    p->walltile = p_go->walltile[j];
    p->cputime = p_go->cputime[j];
    p->mileage = p_go->mileage[j];

    /* Particle to guiding center */
    real B_dB[15], psi[1], rho[2];
    rho[0] = p_go->rho[j];
    B_dB[0] = p_go->B_r[j];
    B_dB[1] = p_go->B_r_dr[j];
    B_dB[2] = p_go->B_r_dphi[j];
    B_dB[3] = p_go->B_r_dz[j];
    B_dB[4] = p_go->B_phi[j];
    B_dB[5] = p_go->B_phi_dr[j];
    B_dB[6] = p_go->B_phi_dphi[j];
    B_dB[7] = p_go->B_phi_dz[j];
    B_dB[8] = p_go->B_z[j];
    B_dB[9] = p_go->B_z_dr[j];
    B_dB[10] = p_go->B_z_dphi[j];
    B_dB[11] = p_go->B_z_dz[j];

    /* Guiding center transformation */
    real ppar, mu;
    if (!err)
    {
        real pr = p->pr;
        real pphi = p->pphi;
        real pz = p->pz;

        gctransform_particle2guidingcenter(
            p->mass, p->charge, B_dB, p->rprt, p->phiprt, p->zprt, pr, pphi, pz,
            &p->r, &p->phi, &p->z, &ppar, &mu, &p->zeta);
    }
    if (!err && p->r <= 0)
    {
        err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);
    }
    if (!err && mu < 0)
    {
        err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);
    }

    if (!err)
    {
        err = Bfield_eval_b_db(B_dB, p->r, p->phi, p->z, p->time, bfield);
    }
    if (!err)
    {
        err = Bfield_eval_psi(psi, p->r, p->phi, p->z, p->time, bfield);
    }
    if (!err)
    {
        err = Bfield_eval_rho(rho, psi[0], bfield);
    }

    real Bnorm = math_normc(B_dB[0], B_dB[4], B_dB[8]);
    p->ekin = physlib_Ekin_ppar(p->mass, mu, ppar, Bnorm);
    p->pitch = physlib_gc_xi(p->mass, mu, ppar, Bnorm);

    /* If marker already has error flag, make sure it is not overwritten here */
    err_t simerr = p_go->err[j];
    if (simerr)
    {
        p->err = simerr;
    }
    else
    {
        p->err = err;
    }
}

/**
 * @brief Convert state into a GC SIMD marker
 *
 * This function assumes markers are drawn from a marker queue where they will
 * be returned at the same position once simulation for this marker has ended,
 *
 * State is converted into a GC marker by simply copying all the necessary
 * fields.
 *
 * The simulation marker is set as running unless there is an error or the state
 * already have active end condition.
 *
 * If state marker already has an error flag, this error flag is returned and
 * nothing is done. If an error occurs within this function, the conversion
 * is terminated and the error is returned. It is the responsibility of the
 * caller to store the error in correct struct.
 *
 * This is a SIMD function.
 *
 * @param p pointer to a state being converted
 * @param i index of this state in the marker queue
 * @param p_gc pointer to SIMD structure where marker is being stored
 * @param j index where in the SIMD structure marker is stored
 * @param bfield pointer to magnetic field data
 */
err_t state_to_gc(
    State *p, size_t i, MarkerGuidingCenter *p_gc, size_t j, Bfield *bfield)
{
    err_t err = p->err;

    if (!err)
    {
        p_gc->r[j] = p->r;
        p_gc->phi[j] = p->phi;
        p_gc->z[j] = p->z;
        p_gc->zeta[j] = p->zeta;

        p_gc->mass[j] = p->mass;
        p_gc->charge[j] = p->charge;
        p_gc->time[j] = p->time;
        p_gc->bounces[j] = 0;
        p_gc->weight[j] = p->weight;
        p_gc->theta[j] = p->theta;
        p_gc->id[j] = p->id;
        p_gc->endcond[j] = p->endcond;
        p_gc->walltile[j] = p->walltile;
        p_gc->mileage[j] = p->mileage;

        real B_dB[15], psi[1], rho[2];
        if (!err)
        {
            err = Bfield_eval_b_db(B_dB, p->r, p->phi, p->z, p->time, bfield);
        }
        if (!err)
        {
            err = Bfield_eval_psi(psi, p->r, p->phi, p->z, p->time, bfield);
        }
        if (!err)
        {
            err = Bfield_eval_rho(rho, psi[0], bfield);
        }
        p_gc->rho[j] = rho[0];
        p_gc->B_r[j] = B_dB[0];
        p_gc->B_r_dr[j] = B_dB[1];
        p_gc->B_r_dphi[j] = B_dB[2];
        p_gc->B_r_dz[j] = B_dB[3];
        p_gc->B_phi[j] = B_dB[4];
        p_gc->B_phi_dr[j] = B_dB[5];
        p_gc->B_phi_dphi[j] = B_dB[6];
        p_gc->B_phi_dz[j] = B_dB[7];
        p_gc->B_z[j] = B_dB[8];
        p_gc->B_z_dr[j] = B_dB[9];
        p_gc->B_z_dphi[j] = B_dB[10];
        p_gc->B_z_dz[j] = B_dB[11];

        real Bnorm = math_normc(B_dB[0], B_dB[4], B_dB[8]);
        real gamma = physlib_gamma_Ekin(p->mass, p->ekin);
        real pnorm = physlib_pnorm_gamma(p->mass, gamma);
        p_gc->mu[j] = physlib_gc_mu(p->mass, pnorm, p->pitch, Bnorm);
        p_gc->ppar[j] = phys_ppar_Ekin(p->mass, p->ekin, p_gc->mu[j], Bnorm);

        p_gc->running[j] = 1;
        if (p->endcond)
        {
            p_gc->running[j] = 0;
        }
        p_gc->cputime[j] = p->cputime;
        p_gc->index[j] = i;
        p_gc->err[j] = 0;
    }

    return err;
}

/**
 * @brief Convert GC to state
 *
 * This function converts GC simulation marker to a marker state. No new state
 * structure is initialized, but an existing one is filled with parameters
 * corresponding to the converted marker.
 *
 * This is a SIMD function.
 *
 * @param p_gc pointer to SIMD structure
 * @param j SIMD index of the marker to be converted
 * @param p pointer to state whose fields are used to store the converted marker
 * @param bfield pointer to magnetic field data
 */
void marker_gc_to_state(
    MarkerGuidingCenter *p_gc, size_t j, State *p, Bfield *bfield)
{
    err_t err = 0;

    p->r = p_gc->r[j];
    p->phi = p_gc->phi[j];
    p->z = p_gc->z[j];

    p->mass = p_gc->mass[j];
    p->charge = p_gc->charge[j];
    p->time = p_gc->time[j];
    p->weight = p_gc->weight[j];
    p->id = p_gc->id[j];
    p->cputime = p_gc->cputime[j];
    p->theta = p_gc->theta[j];
    p->endcond = p_gc->endcond[j];
    p->walltile = p_gc->walltile[j];
    p->mileage = p_gc->mileage[j];

    /* Guiding center to particle transformation */
    real B_dB[15];
    B_dB[0] = p_gc->B_r[j];
    B_dB[1] = p_gc->B_r_dr[j];
    B_dB[2] = p_gc->B_r_dphi[j];
    B_dB[3] = p_gc->B_r_dz[j];
    B_dB[4] = p_gc->B_phi[j];
    B_dB[5] = p_gc->B_phi_dr[j];
    B_dB[6] = p_gc->B_phi_dphi[j];
    B_dB[7] = p_gc->B_phi_dz[j];
    B_dB[8] = p_gc->B_z[j];
    B_dB[9] = p_gc->B_z_dr[j];
    B_dB[10] = p_gc->B_z_dphi[j];
    B_dB[11] = p_gc->B_z_dz[j];

    real Bnorm = math_normc(p_gc->B_r[j], p_gc->B_phi[j], p_gc->B_z[j]);
    p->ekin =
        physlib_Ekin_ppar(p_gc->mass[j], p_gc->mu[j], p_gc->ppar[j], Bnorm);
    p->pitch = physlib_gc_xi(p_gc->mass[j], p_gc->mu[j], p_gc->ppar[j], Bnorm);
    p->zeta = p_gc->zeta[j];

    real pr, pphi, pz;
    if (!err)
    {
        real pparprt, muprt, zetaprt;
        gctransform_guidingcenter2particle(
            p->mass, p->charge, B_dB, p->r, p->phi, p->z, p_gc->ppar[j],
            p_gc->mu[j], p->zeta, &p->rprt, &p->phiprt, &p->zprt, &pparprt,
            &muprt, &zetaprt);

        Bfield_eval_b_db(B_dB, p->rprt, p->phiprt, p->zprt, p->time, bfield);

        gctransform_pparmuzeta2prpphipz(
            p->mass, p->charge, B_dB, p->phiprt, pparprt, muprt, zetaprt, &pr,
            &pphi, &pz);
    }
    if (!err && p->rprt <= 0)
    {
        err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);
    }

    if (!err)
    {
        p->pr = pr;
        p->pphi = pphi;
        p->pz = pz;
    }

    /* If marker already has error flag, make sure it is not overwritten here */
    err_t simerr = p_gc->err[j];
    if (simerr)
    {
        p->err = simerr;
    }
    else
    {
        p->err = err;
    }
}

/**
 * @brief Convert state to a ML SIMD marker
 *
 * This function assumes markers are drawn from a marker queue where they will
 * be returned at the same position once simulation for this marker has ended,
 *
 * State is converted into a ML marker by using particle position as
 * initial position. Other fields are simply copied exceot pitch, which is
 * calculated from vpar as pitch = 2*(vpar >= 0) - 1, i.e., sign of vpar
 * determines sign of pitch.
 *
 * The simulation marker is set as running unless there is an error or the state
 * already have active end condition.
 *
 * If state marker already has an error flag, this error flag is returned and
 * nothing is done. If an error occurs within this function, the conversion
 * is terminated and the error is returned. It is the responsibility of the
 * caller to store the error in correct struct.
 *
 * This is a SIMD function.
 *
 * @todo A minor thing but it would be better if guiding center position were
 *       used instead. This would require evaluation of the magnetic field.
 *
 * @param p pointer to a state being converted
 * @param i index of this state in the marker queue
 * @param p_fl pointer to SIMD structure where marker is being stored
 * @param j index where in the SIMD structure marker is stored
 * @param bfield pointer to magnetic field data
 */
err_t state_to_fl(State *p, size_t i, MarkerFieldLine *p_fl, size_t j, Bfield *bfield)
{
    err_t err = p->err;

    if (!err)
    {
        p_fl->r[j] = p->r;
        p_fl->phi[j] = p->phi;
        p_fl->z[j] = p->z;

        p_fl->pitch[j] = 2 * (p->pitch >= 0) - 1.0;
        p_fl->time[j] = p->time;
        p_fl->id[j] = p->id;
        p_fl->cputime[j] = p->cputime;
        p_fl->theta[j] = p->theta;
        p_fl->endcond[j] = p->endcond;
        p_fl->walltile[j] = p->walltile;
        p_fl->mileage[j] = p->mileage;

        real B_dB[15], psi[1], rho[2];
        if (!err)
        {
            err = Bfield_eval_b_db(B_dB, p->r, p->phi, p->z, p->time, bfield);
        }
        if (!err)
        {
            err = Bfield_eval_psi(psi, p->r, p->phi, p->z, p->time, bfield);
        }
        if (!err)
        {
            err = Bfield_eval_rho(rho, psi[0], bfield);
        }
        p_fl->rho[j] = rho[0];
        p_fl->B_r[j] = B_dB[0];
        p_fl->B_r_dr[j] = B_dB[1];
        p_fl->B_r_dphi[j] = B_dB[2];
        p_fl->B_r_dz[j] = B_dB[3];
        p_fl->B_phi[j] = B_dB[4];
        p_fl->B_phi_dr[j] = B_dB[5];
        p_fl->B_phi_dphi[j] = B_dB[6];
        p_fl->B_phi_dz[j] = B_dB[7];
        p_fl->B_z[j] = B_dB[8];
        p_fl->B_z_dr[j] = B_dB[9];
        p_fl->B_z_dphi[j] = B_dB[10];
        p_fl->B_z_dz[j] = B_dB[11];

        p_fl->running[j] = 1;
        if (p->endcond)
        {
            p_fl->running[j] = 0;
        }
        p_fl->index[j] = i;

        p_fl->err[j] = 0;
    }

    return err;
}

/**
 * @brief Convert ML to state
 *
 * This function converts ML simulation marker to a marker state. No new state
 * structure is initialized, but an existing one is filled with parameters
 * corresponding to the converted marker.
 *
 * Field line position is stored in both particle and guiding center position
 * fields. Direction is stored in vpar field. rdot, phidot, zdot, mu, zeta,
 * mass, and charge are left empty.
 *
 * This is a SIMD function.
 *
 * @param p_fl pointer to SIMD structure
 * @param j SIMD index of the marker to be converted
 * @param p pointer to state whose fields are used to store the converted marker
 * @param bfield pointer to magnetic field data
 */
void marker_fl_to_state(MarkerFieldLine *p_fl, size_t j, State *p)
{
    err_t err = 0;

    p->rprt = p_fl->r[j];
    p->phiprt = p_fl->phi[j];
    p->zprt = p_fl->z[j];
    p->pr = 0;
    p->pphi = 0;
    p->pz = 0;

    p->r = p_fl->r[j];
    p->phi = p_fl->phi[j];
    p->z = p_fl->z[j];
    p->ekin = 0;
    p->pitch = p_fl->pitch[j];
    p->zeta = 0;
    p->mass = 0;
    p->charge = 0;
    p->time = p_fl->time[j];
    p->id = p_fl->id[j];
    p->cputime = p_fl->cputime[j];
    p->theta = p_fl->theta[j];
    p->endcond = p_fl->endcond[j];
    p->walltile = p_fl->walltile[j];
    p->mileage = p_fl->mileage[j];
    p->err = p_fl->err[j];

    /* If marker already has error flag, make sure it is not overwritten here */
    err_t simerr = p_fl->err[j];
    if (simerr)
    {
        p->err = simerr;
    }
    else
    {
        p->err = err;
    }
    if (!simerr && err)
    {
        err = err;
    }
    p->err = err;
}

/**
 * @brief Convert FO struct into a GC struct
 *
 * @deprecated Do we need this as same can be accomplished by turning
 *             FO to state and state to GC?
 *
 * @param p_go  fo SIMD structure being transformed
 * @param j     index where in the SIMD structure marker is stored
 * @param p_gc  gc SIMD structure where marker is transformed
 * @param bfield pointer to magnetic field data
 */
int marker_go_to_gc(
    MarkerGyroOrbit *p_go, size_t j, MarkerGuidingCenter *p_gc, Bfield *bfield)
{
    err_t err = p_go->err[j];
    real axisrz[2];
    int simerr = 0; /* Error has already occurred */
    if (err)
    {
        simerr = 1;
    }
    p_gc->id[j] = p_go->id[j];
    p_gc->index[j] = p_go->index[j];

    real r, phi, z, ppar, mu, zeta, B_dB[15];
    if (!err)
    {
        real Rprt = p_go->r[j];
        real phiprt = p_go->phi[j];
        real zprt = p_go->z[j];
        real pr = p_go->p_r[j];
        real pphi = p_go->p_phi[j];
        real pz = p_go->p_z[j];
        real mass = p_go->mass[j];
        real charge = p_go->charge[j];

        p_gc->mass[j] = p_go->mass[j];
        p_gc->charge[j] = p_go->charge[j];
        p_gc->weight[j] = p_go->weight[j];
        p_gc->time[j] = p_go->time[j];
        p_gc->mileage[j] = p_go->mileage[j];
        p_gc->endcond[j] = p_go->endcond[j];
        p_gc->running[j] = p_go->running[j];
        p_gc->walltile[j] = p_go->walltile[j];
        p_gc->cputime[j] = p_go->cputime[j];

        B_dB[0] = p_go->B_r[j];
        B_dB[1] = p_go->B_r_dr[j];
        B_dB[2] = p_go->B_r_dphi[j];
        B_dB[3] = p_go->B_r_dz[j];
        B_dB[4] = p_go->B_phi[j];
        B_dB[5] = p_go->B_phi_dr[j];
        B_dB[6] = p_go->B_phi_dphi[j];
        B_dB[7] = p_go->B_phi_dz[j];
        B_dB[8] = p_go->B_z[j];
        B_dB[9] = p_go->B_z_dr[j];
        B_dB[10] = p_go->B_z_dphi[j];
        B_dB[11] = p_go->B_z_dz[j];

        /* Guiding center transformation */
        gctransform_particle2guidingcenter(
            mass, charge, B_dB, Rprt, phiprt, zprt, pr, pphi, pz, &r, &phi, &z,
            &ppar, &mu, &zeta);
    }
    if (!err && r <= 0)
    {
        err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);
    }
    if (!err && mu < 0)
    {
        err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);
    }

    real psi[1], rho[2];
    if (!err)
    {
        err = Bfield_eval_b_db(B_dB, r, phi, z, p_go->time[j], bfield);
    }
    if (!err)
    {
        err = Bfield_eval_psi(psi, r, phi, z, p_go->time[j], bfield);
    }
    if (!err)
    {
        err = Bfield_eval_rho(rho, psi[0], bfield);
    }
    if (!err)
    {
        err = Bfield_eval_axis_rz(axisrz, bfield, p_gc->phi[j]);
    }

    if (!err)
    {
        p_gc->r[j] = r;
        p_gc->phi[j] = phi;
        p_gc->z[j] = z;
        p_gc->mu[j] = mu;
        p_gc->zeta[j] = zeta;
        p_gc->ppar[j] = ppar;
        p_gc->rho[j] = rho[0];

        /* Evaluate pol angle so that it is cumulative and at gc position */
        p_gc->theta[j] = p_go->theta[j];
        p_gc->theta[j] += atan2(
            (p_go->r[j] - axisrz[0]) * (p_gc->z[j] - axisrz[1]) -
                (p_go->z[j] - axisrz[1]) * (p_gc->r[j] - axisrz[0]),
            (p_go->r[j] - axisrz[0]) * (p_gc->r[j] - axisrz[0]) +
                (p_go->z[j] - axisrz[1]) * (p_gc->z[j] - axisrz[1]));

        p_gc->B_r[j] = B_dB[0];
        p_gc->B_r_dr[j] = B_dB[1];
        p_gc->B_r_dphi[j] = B_dB[2];
        p_gc->B_r_dz[j] = B_dB[3];

        p_gc->B_phi[j] = B_dB[4];
        p_gc->B_phi_dr[j] = B_dB[5];
        p_gc->B_phi_dphi[j] = B_dB[6];
        p_gc->B_phi_dz[j] = B_dB[7];

        p_gc->B_z[j] = B_dB[8];
        p_gc->B_z_dr[j] = B_dB[9];
        p_gc->B_z_dphi[j] = B_dB[10];
        p_gc->B_z_dz[j] = B_dB[11];
    }
    if (!simerr)
    {
        err = err;
    }
    p_gc->err[j] = err;
    if (p_gc->err[j])
    {
        p_gc->running[j] = 0;
        p_gc->endcond[j] = 0;
    }

    return err > 0;
}

/**
 * @brief Copy FO struct
 *
 * @param p1 SIMD structure for input
 * @param i  index for the copied input
 * @param p2 SIMD structure for output
 * @param j  index for the output slot
 */
void marker_copy_go(MarkerGyroOrbit *p1, size_t i, MarkerGyroOrbit *p2, size_t j)
{
    p2->r[j] = p1->r[i];
    p2->phi[j] = p1->phi[i];
    p2->z[j] = p1->z[i];
    p2->p_r[j] = p1->p_r[i];
    p2->p_phi[j] = p1->p_phi[i];
    p2->p_z[j] = p1->p_z[i];

    p2->time[j] = p1->time[i];
    p2->mileage[j] = p1->mileage[i];
    p2->cputime[j] = p1->cputime[i];
    p2->rho[j] = p1->rho[i];
    p2->weight[j] = p1->weight[i];
    p2->cputime[j] = p1->cputime[i];
    p2->rho[j] = p1->rho[i];
    p2->theta[j] = p1->theta[i];

    p2->mass[j] = p1->mass[i];
    p2->charge[j] = p1->charge[i];
    p2->znum[j] = p1->znum[i];
    p2->anum[j] = p1->anum[i];

    p2->id[j] = p1->id[i];
    p2->bounces[j] = p1->bounces[i];
    p2->running[j] = p1->running[i];
    p2->endcond[j] = p1->endcond[i];
    p2->walltile[j] = p1->walltile[i];

    p2->B_r[j] = p1->B_r[i];
    p2->B_phi[j] = p1->B_phi[i];
    p2->B_z[j] = p1->B_z[i];

    p2->B_r_dr[j] = p1->B_r_dr[i];
    p2->B_r_dphi[j] = p1->B_r_dphi[i];
    p2->B_r_dz[j] = p1->B_r_dz[i];

    p2->B_phi_dr[j] = p1->B_phi_dr[i];
    p2->B_phi_dphi[j] = p1->B_phi_dphi[i];
    p2->B_phi_dz[j] = p1->B_phi_dz[i];

    p2->B_z_dr[j] = p1->B_z_dr[i];
    p2->B_z_dphi[j] = p1->B_z_dphi[i];
    p2->B_z_dz[j] = p1->B_z_dz[i];
}

/**
 * @brief Copy GC struct
 *
 * @param p1 SIMD structure for input
 * @param i  index for the copied input
 * @param p2 SIMD structure for output
 * @param j  index for the output slot
 */
void marker_copy_gc(
    MarkerGuidingCenter *p1, size_t i, MarkerGuidingCenter *p2, size_t j)
{
    p2->r[j] = p1->r[i];
    p2->phi[j] = p1->phi[i];
    p2->z[j] = p1->z[i];
    p2->ppar[j] = p1->ppar[i];
    p2->mu[j] = p1->mu[i];
    p2->zeta[j] = p1->zeta[i];

    p2->time[j] = p1->time[i];
    p2->mileage[j] = p1->mileage[i];
    p2->weight[j] = p1->weight[i];
    p2->cputime[j] = p1->cputime[i];
    p2->rho[j] = p1->rho[i];
    p2->theta[j] = p1->theta[i];

    p2->mass[j] = p1->mass[i];
    p2->charge[j] = p1->charge[i];

    p2->id[j] = p1->id[i];
    p2->bounces[j] = p1->bounces[i];
    p2->running[j] = p1->running[i];
    p2->endcond[j] = p1->endcond[i];
    p2->walltile[j] = p1->walltile[i];

    p2->B_r[j] = p1->B_r[i];
    p2->B_phi[j] = p1->B_phi[i];
    p2->B_z[j] = p1->B_z[i];

    p2->B_r_dr[j] = p1->B_r_dr[i];
    p2->B_r_dphi[j] = p1->B_r_dphi[i];
    p2->B_r_dz[j] = p1->B_r_dz[i];

    p2->B_phi_dr[j] = p1->B_phi_dr[i];
    p2->B_phi_dphi[j] = p1->B_phi_dphi[i];
    p2->B_phi_dz[j] = p1->B_phi_dz[i];

    p2->B_z_dr[j] = p1->B_z_dr[i];
    p2->B_z_dphi[j] = p1->B_z_dphi[i];
    p2->B_z_dz[j] = p1->B_z_dz[i];
}

/**
 * @brief Copy ML struct
 *
 * @param p1 SIMD structure for input
 * @param i  index for the copied input
 * @param p2 SIMD structure for output
 * @param j  index for the output slot
 */
void marker_copy_fl(MarkerFieldLine *p1, size_t i, MarkerFieldLine *p2, size_t j)
{
    p2->r[j] = p1->r[i];
    p2->phi[j] = p1->phi[i];
    p2->z[j] = p1->z[i];
    p2->pitch[j] = p1->pitch[i];

    p2->time[j] = p1->time[i];
    p2->mileage[j] = p1->mileage[i];
    p2->cputime[j] = p1->cputime[i];
    p2->rho[j] = p1->rho[i];
    p2->theta[j] = p1->theta[i];

    p2->id[j] = p1->id[i];
    p2->running[j] = p1->running[i];
    p2->endcond[j] = p1->endcond[i];
    p2->walltile[j] = p1->walltile[i];

    p2->B_r[j] = p1->B_r[i];
    p2->B_phi[j] = p1->B_phi[i];
    p2->B_z[j] = p1->B_z[i];

    p2->B_r_dr[j] = p1->B_r_dr[i];
    p2->B_r_dphi[j] = p1->B_r_dphi[i];
    p2->B_r_dz[j] = p1->B_r_dz[i];

    p2->B_phi_dr[j] = p1->B_phi_dr[i];
    p2->B_phi_dphi[j] = p1->B_phi_dphi[i];
    p2->B_phi_dz[j] = p1->B_phi_dz[i];

    p2->B_z_dr[j] = p1->B_z_dr[i];
    p2->B_z_dphi[j] = p1->B_z_dphi[i];
    p2->B_z_dz[j] = p1->B_z_dz[i];
}

/**
 * @brief Offload particle struct to GPU.
 *
 * @param p pointer to the particle struct to be offloaded.
 */
void marker_offload_go(MarkerGyroOrbit *p)
{
    SUPPRESS_UNUSED_WARNING(p);
    GPU_MAP_TO_DEVICE(
        p [0:1], p->running [0:p->n_mrk], p->r [0:p->n_mrk],
        p->phi [0:p->n_mrk], p->p_r [0:p->n_mrk], p->p_phi [0:p->n_mrk],
        p->p_z [0:p->n_mrk], p->mileage [0:p->n_mrk], p->z [0:p->n_mrk],
        p->charge [0:p->n_mrk], p->mass [0:p->n_mrk], p->B_r [0:p->n_mrk],
        p->B_r_dr [0:p->n_mrk], p->B_r_dphi [0:p->n_mrk],
        p->B_r_dz [0:p->n_mrk], p->B_phi [0:p->n_mrk], p->B_phi_dr [0:p->n_mrk],
        p->B_phi_dphi [0:p->n_mrk], p->B_phi_dz [0:p->n_mrk],
        p->B_z [0:p->n_mrk], p->B_z_dr [0:p->n_mrk], p->B_z_dphi [0:p->n_mrk],
        p->B_z_dz [0:p->n_mrk], p->rho [0:p->n_mrk], p->theta [0:p->n_mrk],
        p->err [0:p->n_mrk], p->time [0:p->n_mrk], p->weight [0:p->n_mrk],
        p->cputime [0:p->n_mrk], p->id [0:p->n_mrk], p->endcond [0:p->n_mrk],
        p->walltile [0:p->n_mrk], p->index [0:p->n_mrk], p->znum [0:p->n_mrk],
        p->anum [0:p->n_mrk], p->bounces [0:p->n_mrk])
}

/**
 * @brief Onload particle struct from the GPU.
 *
 * @param p pointer to the particle struct to be onloaded.
 */
void marker_onload_go(MarkerGyroOrbit *p)
{
    SUPPRESS_UNUSED_WARNING(p);
    GPU_UPDATE_FROM_DEVICE(
        p->running [0:p->n_mrk], p->r [0:p->n_mrk], p->phi [0:p->n_mrk],
        p->p_r [0:p->n_mrk], p->p_phi [0:p->n_mrk], p->p_z [0:p->n_mrk],
        p->mileage [0:p->n_mrk], p->z [0:p->n_mrk], p->charge [0:p->n_mrk],
        p->mass [0:p->n_mrk], p->B_r [0:p->n_mrk], p->B_r_dr [0:p->n_mrk],
        p->B_r_dphi [0:p->n_mrk], p->B_r_dz [0:p->n_mrk], p->B_phi [0:p->n_mrk],
        p->B_phi_dr [0:p->n_mrk], p->B_phi_dphi [0:p->n_mrk],
        p->B_phi_dz [0:p->n_mrk], p->B_z [0:p->n_mrk], p->B_z_dr [0:p->n_mrk],
        p->B_z_dphi [0:p->n_mrk], p->B_z_dz [0:p->n_mrk], p->rho [0:p->n_mrk],
        p->theta [0:p->n_mrk], p->err [0:p->n_mrk], p->time [0:p->n_mrk],
        p->weight [0:p->n_mrk], p->cputime [0:p->n_mrk], p->id [0:p->n_mrk],
        p->endcond [0:p->n_mrk], p->walltile [0:p->n_mrk],
        p->index [0:p->n_mrk], p->znum [0:p->n_mrk], p->anum [0:p->n_mrk],
        p->bounces [0:p->n_mrk])
}
