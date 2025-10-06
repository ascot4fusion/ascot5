/**
 * Implements marker.h.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defines.h"
#include "errors.h"
#include "consts.h"
#include "mathlib.h"
#include "physlib.h"
#include "gctransform.h"
#include "particle.h"
#include "bfield.h"
#include "efield.h"

/**
 * @brief Allocates struct representing particle markers
 *
 * Size used for memory allocation is NSIMD for CPU run and the total number
 * of particles for GPU.
 *
 * @param p_fo marker struct to allocate
 * @param nmrk the number of markers that the struct represents
 */
void particle_allocate_fo(particle_simd_fo* p_fo, int nmrk){
    p_fo->r      = malloc(nmrk * sizeof(p_fo->r)     );
    p_fo->phi    = malloc(nmrk * sizeof(p_fo->phi)   );
    p_fo->z      = malloc(nmrk * sizeof(p_fo->z)     );
    p_fo->p_r    = malloc(nmrk * sizeof(p_fo->p_r)   );
    p_fo->p_phi  = malloc(nmrk * sizeof(p_fo->p_phi) );
    p_fo->p_z    = malloc(nmrk * sizeof(p_fo->p_z)   );
    p_fo->mass   = malloc(nmrk * sizeof(p_fo->mass)  );
    p_fo->charge = malloc(nmrk * sizeof(p_fo->charge));
    p_fo->time   = malloc(nmrk * sizeof(p_fo->time)  );
    p_fo->znum   = malloc(nmrk * sizeof(p_fo->znum)  );
    p_fo->anum   = malloc(nmrk * sizeof(p_fo->anum)  );

    /* Magnetic field data */
    p_fo->B_r        = malloc(nmrk * sizeof(p_fo->B_r)       );
    p_fo->B_phi      = malloc(nmrk * sizeof(p_fo->B_phi)     );
    p_fo->B_z        = malloc(nmrk * sizeof(p_fo->B_z)       );
    p_fo->B_r_dr     = malloc(nmrk * sizeof(p_fo->B_r_dr)    );
    p_fo->B_phi_dr   = malloc(nmrk * sizeof(p_fo->B_phi_dr)  );
    p_fo->B_z_dr     = malloc(nmrk * sizeof(p_fo->B_z_dr)    );
    p_fo->B_r_dphi   = malloc(nmrk * sizeof(p_fo->B_r_dphi)  );
    p_fo->B_phi_dphi = malloc(nmrk * sizeof(p_fo->B_phi_dphi));
    p_fo->B_z_dphi   = malloc(nmrk * sizeof(p_fo->B_z_dphi)  );
    p_fo->B_r_dz     = malloc(nmrk * sizeof(p_fo->B_r_dz)    );
    p_fo->B_phi_dz   = malloc(nmrk * sizeof(p_fo->B_phi_dz)  );
    p_fo->B_z_dz     = malloc(nmrk * sizeof(p_fo->B_z_dz)    );

    /* Quantities used in diagnostics */
    p_fo->bounces = malloc(nmrk * sizeof(p_fo->bounces));
    p_fo->weight  = malloc(nmrk * sizeof(p_fo->weight) );
    p_fo->cputime = malloc(nmrk * sizeof(p_fo->cputime));
    p_fo->rho     = malloc(nmrk * sizeof(p_fo->rho)    );
    p_fo->theta   = malloc(nmrk * sizeof(p_fo->theta)  );

    p_fo->id       = malloc(nmrk * sizeof(p_fo->id)      );
    p_fo->endcond  = malloc(nmrk * sizeof(p_fo->endcond) );
    p_fo->walltile = malloc(nmrk * sizeof(p_fo->walltile));

    /* Meta data */
    p_fo->mileage = malloc(nmrk * sizeof(p_fo->mileage));

    p_fo->running = malloc(nmrk * sizeof(p_fo->running));

    p_fo->err   = malloc(nmrk * sizeof(p_fo->err)  );
    p_fo->index = malloc(nmrk * sizeof(p_fo->index));
    p_fo->n_mrk = nmrk;
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
 * @param p_fo  pointer to particle_simd_fo array where dummy is placed
 * @param j     index of the new dummy in the SIMD arrays
 */
void particle_to_fo_dummy(particle_simd_fo* p_fo, int j){
    p_fo->r[j]        = 1;
    p_fo->phi[j]      = 1;
    p_fo->z[j]        = 1;
    p_fo->p_r[j]      = 1;
    p_fo->p_phi[j]    = 1;
    p_fo->p_z[j]      = 1;
    p_fo->mass[j]     = 1;
    p_fo->charge[j]   = 1;
    p_fo->znum[j]     = 1;
    p_fo->anum[j]     = 1;
    p_fo->bounces[j]  = 0;
    p_fo->weight[j]   = 0;
    p_fo->time[j]     = 0;
    p_fo->id[j]       = -1;
    p_fo->mileage[j]  = 0;
    p_fo->running[j]  = 0;
    p_fo->endcond[j]  = 0;
    p_fo->walltile[j] = 0;
    p_fo->B_r[j]      = 1;
    p_fo->B_phi[j]    = 1;
    p_fo->B_z[j]      = 1;
    p_fo->index[j]    = -1;
    p_fo->err[j]      = 0;
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
 * @param p_gc  pointer to particle_simd_gc array where dummy is placed
 * @param j     index of the new dummy in the SIMD arrays
 */
void particle_to_gc_dummy(particle_simd_gc* p_gc, int j) {
    p_gc->r[j]          = 1;
    p_gc->phi[j]        = 1;
    p_gc->z[j]          = 1;
    p_gc->ppar[j]       = 1;
    p_gc->mu[j]         = 1;
    p_gc->zeta[j]       = 1;
    p_gc->mass[j]       = 1;
    p_gc->charge[j]     = 1;
    p_gc->time[j]       = 0;
    p_gc->bounces[j]    = 0;
    p_gc->weight[j]     = 0;
    p_gc->id[j]         = -1;
    p_gc->mileage[j]    = 0;
    p_gc->B_r[j]        = 1;
    p_gc->B_r_dr[j]     = 1;
    p_gc->B_r_dphi[j]   = 1;
    p_gc->B_r_dz[j]     = 1;

    p_gc->B_phi[j]      = 1;
    p_gc->B_phi_dr[j]   = 1;
    p_gc->B_phi_dphi[j] = 1;
    p_gc->B_phi_dz[j]   = 1;

    p_gc->B_z[j]        = 1;
    p_gc->B_z_dr[j]     = 1;
    p_gc->B_z_dphi[j]   = 1;
    p_gc->B_z_dz[j]     = 1;
    p_gc->index[j]      = -1;
    p_gc->err[j]        = 0;
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
 * @param p_ml  pointer to particle_simd_ml array where dummy is placed
 * @param j     index of the new dummy in the SIMD arrays
 */
void particle_to_ml_dummy(particle_simd_ml* p_ml, int j){
    p_ml->r[j]          = 1;
    p_ml->phi[j]        = 1;
    p_ml->z[j]          = 1;
    p_ml->time[j]       = 0;
    p_ml->id[j]         = -1;
    p_ml->mileage[j]    = 0;
    p_ml->B_r[j]        = 1;
    p_ml->B_r_dr[j]     = 1;
    p_ml->B_r_dphi[j]   = 1;
    p_ml->B_r_dz[j]     = 1;

    p_ml->B_phi[j]      = 1;
    p_ml->B_phi_dr[j]   = 1;
    p_ml->B_phi_dphi[j] = 1;
    p_ml->B_phi_dz[j]   = 1;

    p_ml->B_z[j]        = 1;
    p_ml->B_z_dr[j]     = 1;
    p_ml->B_z_dphi[j]   = 1;
    p_ml->B_z_dz[j]     = 1;
    p_ml->index[j]      = -1;
    p_ml->err[j]        = 0;
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
int particle_cycle_fo(particle_queue* q, particle_simd_fo* p,
                      Bfield *bfield, int* cycle) {

    /* Loop over markers.
     * A SIMD loop is not possible as we modify the queue. */
    for(size_t i = 0; i < p->n_mrk; i++) {

        /* First check whether we should pick a new marker */
        int newmarker = 0;

        /* 1. There are markers in queue and this marker is dummy */
        if(p->id[i] < 0 && q->next < q->n) {
            newmarker = 1;
        }

        /* 2. This marker has finished simulation */
        if(!p->running[i] && p->id[i] >= 0) {
            /* Convert finished marker to state
             * and store it back to the queue */
            particle_fo_to_state(p, i, q->p[p->index[i]], bfield);
            newmarker = 1;
            #pragma omp critical
            q->finished++;
        }

        /* Init a new marker if one is needed */
        cycle[i] = 0;
        while(newmarker) {
            /* Get the next unsimulated marker from the queue */
            int i_prt;
            #pragma omp critical
            i_prt = q->next++;

            if(i_prt >= q->n) {
                /* The queue is empty, place a dummy marker here */
                p->running[i] = 0;
                p->id[i] = -1;
                cycle[i] = -1;
                newmarker = 0;
            }
            else if(q->p[i_prt]->endcond) {
                /* This marker already has an active end condition. Try next. */
                #pragma omp critical
                q->finished++;
            }
            else {
                /* Try to convert marker state to a simulation marker */
                a5err err = particle_state_to_fo(
                    q->p[i_prt], i_prt, p, i, bfield);
                if(err) {
                    /* Failed! Mark the marker candidate as finished
                     * and try with a new marker state*/
                    q->p[i_prt]->err = err;
                    #pragma omp critical
                    q->finished++;
                }
                else {
                    /* Success! We are good to go. */
                    cycle[i] = 1;
                    newmarker = 0;
                }
            }
        }
    }

    int n_running = 0;
    #pragma omp simd reduction(+:n_running)
    for(size_t i = 0; i < p->n_mrk; i++) {
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
int particle_cycle_gc(particle_queue* q, particle_simd_gc* p,
                      Bfield *bfield, int* cycle) {

    /* Loop over markers.
     * A SIMD loop is not possible as we modify the queue. */
    for(int i = 0; i < NSIMD; i++) {

        /* First check whether we should pick a new marker */
        int newmarker = 0;

        /* 1. There are markers in queue and this marker is dummy */
        if(p->id[i] < 0 && q->next < q->n) {
            newmarker = 1;
        }

        /* 2. This marker has finished simulation */
        if(!p->running[i] && p->id[i] >= 0) {
            /* Convert finished marker to state
             * and store it back to the queue */
            particle_gc_to_state(p, i, q->p[p->index[i]], bfield);
            newmarker = 1;
            #pragma omp critical
            q->finished++;
        }

        /* Init a new marker if one is needed */
        cycle[i] = 0;
        while(newmarker) {
            /* Get the next unsimulated marker from the queue */
            int i_prt;
            #pragma omp critical
            i_prt = q->next++;

            if(i_prt >= q->n) {
                /* The queue is empty, place a dummy marker here */
                p->running[i] = 0;
                p->id[i] = -1;
                cycle[i] = -1;
                newmarker = 0;
            }
            else {
                /* Try to convert marker state to a simulation marker */
                a5err err = particle_state_to_gc(
                    q->p[i_prt], i_prt, p, i, bfield);
                if(err) {
                    /* Failed! Mark the marker candidate as finished
                     * and try with a new marker state*/
                    q->p[i_prt]->err = err;
                    #pragma omp critical
                    q->finished++;
                }
                else {
                    /* Success! We are good to go. */
                    cycle[i] = 1;
                    newmarker = 0;
                }
            }
        }
    }

    int n_running = 0;
    #pragma omp simd reduction(+:n_running)
    for(int i = 0; i < NSIMD; i++) {
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
int particle_cycle_ml(particle_queue* q, particle_simd_ml* p,
                      Bfield *bfield, int* cycle) {

    /* Loop over markers.
     * A SIMD loop is not possible as we modify the queue. */
    for(int i = 0; i < NSIMD; i++) {

        /* First check whether we should pick a new marker */
        int newmarker = 0;

        /* 1. There are markers in queue and this marker is dummy */
        if(p->id[i] < 0 && q->next < q->n) {
            newmarker = 1;
        }

        /* 2. This marker has finished simulation */
        if(!p->running[i] && p->id[i] >= 0) {
            /* Convert finished marker to state
             * and store it back to the queue */
            particle_ml_to_state(p, i, q->p[p->index[i]]);
            newmarker = 1;
            #pragma omp critical
            q->finished++;
        }

        /* Init a new marker if one is needed */
        cycle[i] = 0;
        while(newmarker) {
            /* Get the next unsimulated marker from the queue */
            int i_prt;
            #pragma omp critical
            i_prt = q->next++;

            if(i_prt >= q->n) {
                /* The queue is empty, place a dummy marker here */
                p->running[i] = 0;
                p->id[i] = -1;
                cycle[i] = -1;
                newmarker = 0;
            }
            else {
                /* Try to convert marker state to a simulation marker */
                a5err err = particle_state_to_ml(
                    q->p[i_prt], i_prt, p, i, bfield);
                if(err) {
                    /* Failed! Mark the marker candidate as finished
                     * and try with a new marker state*/
                    q->p[i_prt]->err = err;
                    #pragma omp critical
                    q->finished++;
                }
                else {
                    /* Success! We are good to go. */
                    cycle[i] = 1;
                    newmarker = 0;
                }
            }
        }
    }

    int n_running = 0;
    #pragma omp simd reduction(+:n_running)
    for(int i = 0; i < NSIMD; i++) {
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
 * @param p_fo pointer to SIMD structure where marker is being stored
 * @param j index where in the SIMD structure marker is stored
 * @param bfield pointer to magnetic field data
 */
a5err particle_state_to_fo(particle_state* p, int i, particle_simd_fo* p_fo,
                           int j, Bfield *bfield) {
    a5err err = p->err;

    if(!err) {
        p_fo->r[j]          = p->rprt;
        p_fo->phi[j]        = p->phiprt;
        p_fo->z[j]          = p->zprt;
        p_fo->p_r[j]        = p->pr;
        p_fo->p_phi[j]      = p->pphi;
        p_fo->p_z[j]        = p->pz;

        p_fo->mass[j]       = p->mass;
        p_fo->charge[j]     = p->charge;
        p_fo->znum[j]       = p->znum;
        p_fo->anum[j]       = p->anum;
        p_fo->bounces[j]    = 0;
        p_fo->weight[j]     = p->weight;
        p_fo->time[j]       = p->time;
        p_fo->theta[j]      = p->theta;
        p_fo->id[j]         = p->id;
        p_fo->endcond[j]    = p->endcond;
        p_fo->walltile[j]   = p->walltile;
        p_fo->mileage[j]    = p->mileage;
    }

    /* Magnetic field stored in state is for the gc position */
    real B_dB[15], psi[1], rho[2];
    if(!err) {
        err = B_field_eval_B_dB(B_dB, p->rprt, p->phiprt, p->zprt, p->time,
                                bfield);
    }
    if(!err) {
        err = B_field_eval_psi(psi, p->rprt, p->phiprt, p->zprt, p->time,
                               bfield);
    }
    if(!err) {
        err = B_field_eval_rho(rho, psi[0], bfield);
    }

    if(!err) {
        p_fo->rho[j]        = rho[0];

        p_fo->B_r[j]        = B_dB[0];
        p_fo->B_r_dr[j]     = B_dB[1];
        p_fo->B_r_dphi[j]   = B_dB[2];
        p_fo->B_r_dz[j]     = B_dB[3];

        p_fo->B_phi[j]      = B_dB[4];
        p_fo->B_phi_dr[j]   = B_dB[5];
        p_fo->B_phi_dphi[j] = B_dB[6];
        p_fo->B_phi_dz[j]   = B_dB[7];

        p_fo->B_z[j]        = B_dB[8];
        p_fo->B_z_dr[j]     = B_dB[9];
        p_fo->B_z_dphi[j]   = B_dB[10];
        p_fo->B_z_dz[j]     = B_dB[11];

        p_fo->running[j] = 1;
        if(p->endcond) {
            p_fo->running[j] = 0;
        }
        p_fo->cputime[j] = p->cputime;
        p_fo->index[j]   = i;

        p_fo->err[j] = 0;
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
 * @param p_fo pointer to SIMD structure
 * @param j SIMD index of the marker to be converted
 * @param p pointer to state whose fields are used to store the converted marker
 * @param bfield pointer to magnetic field data
 */
void particle_fo_to_state(particle_simd_fo* p_fo, int j, particle_state* p,
                          Bfield *bfield) {
    a5err err = 0;

    p->rprt       = p_fo->r[j];
    p->phiprt     = p_fo->phi[j];
    p->zprt       = p_fo->z[j];
    p->pr         = p_fo->p_r[j];
    p->pphi       = p_fo->p_phi[j];
    p->pz         = p_fo->p_z[j];

    p->mass       = p_fo->mass[j];
    p->charge     = p_fo->charge[j];
    p->znum       = p_fo->znum[j];
    p->anum       = p_fo->anum[j];
    p->weight     = p_fo->weight[j];
    p->time       = p_fo->time[j];
    p->theta      = p_fo->theta[j];
    p->id         = p_fo->id[j];
    p->endcond    = p_fo->endcond[j];
    p->walltile   = p_fo->walltile[j];
    p->cputime    = p_fo->cputime[j];
    p->mileage    = p_fo->mileage[j];

    /* Particle to guiding center */
    real B_dB[15], psi[1], rho[2];
    rho[0]        = p_fo->rho[j];
    B_dB[0]       = p_fo->B_r[j];
    B_dB[1]       = p_fo->B_r_dr[j];
    B_dB[2]       = p_fo->B_r_dphi[j];
    B_dB[3]       = p_fo->B_r_dz[j];
    B_dB[4]       = p_fo->B_phi[j];
    B_dB[5]       = p_fo->B_phi_dr[j];
    B_dB[6]       = p_fo->B_phi_dphi[j];
    B_dB[7]       = p_fo->B_phi_dz[j];
    B_dB[8]       = p_fo->B_z[j];
    B_dB[9]       = p_fo->B_z_dr[j];
    B_dB[10]      = p_fo->B_z_dphi[j];
    B_dB[11]      = p_fo->B_z_dz[j];

    /* Guiding center transformation */
    real ppar, mu;
    if(!err) {
        real pr   = p->pr;
        real pphi = p->pphi;
        real pz   = p->pz;

        gctransform_particle2guidingcenter(
            p->mass, p->charge, B_dB,
            p->rprt, p->phiprt, p->zprt, pr , pphi, pz,
            &p->r, &p->phi, &p->z, &ppar, &mu, &p->zeta);
    }
    if(!err && p->r <= 0)  {err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);}
    if(!err && mu < 0)  {err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);}

    if(!err) {
        err = B_field_eval_B_dB(B_dB, p->r, p->phi, p->z, p->time, bfield);
    }
    if(!err) {
        err = B_field_eval_psi(psi, p->r, p->phi, p->z, p->time, bfield);
    }
    if(!err) {
        err = B_field_eval_rho(rho, psi[0], bfield);
    }

    real Bnorm = math_normc(B_dB[0], B_dB[4], B_dB[8]);
    p->ekin  = physlib_Ekin_ppar(p->mass, mu, ppar, Bnorm);
    p->pitch = physlib_gc_xi(p->mass, mu, ppar, Bnorm);


    /* If marker already has error flag, make sure it is not overwritten here */
    a5err simerr  = p_fo->err[j];
    if(simerr) {
        p->err = simerr;
    }
    else {
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
a5err particle_state_to_gc(particle_state* p, int i, particle_simd_gc* p_gc,
                           int j, Bfield *bfield) {
    a5err err = p->err;

    if(!err) {
        p_gc->r[j]          = p->r;
        p_gc->phi[j]        = p->phi;
        p_gc->z[j]          = p->z;
        p_gc->zeta[j]       = p->zeta;

        p_gc->mass[j]       = p->mass;
        p_gc->charge[j]     = p->charge;
        p_gc->time[j]       = p->time;
        p_gc->bounces[j]    = 0;
        p_gc->weight[j]     = p->weight;
        p_gc->theta[j]      = p->theta;
        p_gc->id[j]         = p->id;
        p_gc->endcond[j]    = p->endcond;
        p_gc->walltile[j]   = p->walltile;
        p_gc->mileage[j]    = p->mileage;

        real B_dB[15], psi[1], rho[2];
        if(!err) {
            err = B_field_eval_B_dB(B_dB, p->r, p->phi, p->z, p->time, bfield);
        }
        if(!err) {
            err = B_field_eval_psi(psi, p->r, p->phi, p->z, p->time, bfield);
        }
        if(!err) {
            err = B_field_eval_rho(rho, psi[0], bfield);
        }
        p_gc->rho[j]        = rho[0];
        p_gc->B_r[j]        = B_dB[0];
        p_gc->B_r_dr[j]     = B_dB[1];
        p_gc->B_r_dphi[j]   = B_dB[2];
        p_gc->B_r_dz[j]     = B_dB[3];
        p_gc->B_phi[j]      = B_dB[4];
        p_gc->B_phi_dr[j]   = B_dB[5];
        p_gc->B_phi_dphi[j] = B_dB[6];
        p_gc->B_phi_dz[j]   = B_dB[7];
        p_gc->B_z[j]        = B_dB[8];
        p_gc->B_z_dr[j]     = B_dB[9];
        p_gc->B_z_dphi[j]   = B_dB[10];
        p_gc->B_z_dz[j]     = B_dB[11];


        real Bnorm = math_normc(B_dB[0], B_dB[4], B_dB[8]);
        real gamma = physlib_gamma_Ekin(p->mass, p->ekin);
        real pnorm = physlib_pnorm_gamma(p->mass, gamma);
        p_gc->mu[j] = physlib_gc_mu(p->mass, pnorm, p->pitch, Bnorm);
        p_gc->ppar[j] = phys_ppar_Ekin(p->mass, p->ekin, p_gc->mu[j], Bnorm);

        p_gc->running[j] = 1;
        if(p->endcond) {
            p_gc->running[j] = 0;
        }
        p_gc->cputime[j] = p->cputime;
        p_gc->index[j]   = i;
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
void particle_gc_to_state(particle_simd_gc* p_gc, int j, particle_state* p,
                          Bfield *bfield) {
    a5err err = 0;

    p->r          = p_gc->r[j];
    p->phi        = p_gc->phi[j];
    p->z          = p_gc->z[j];

    p->mass       = p_gc->mass[j];
    p->charge     = p_gc->charge[j];
    p->time       = p_gc->time[j];
    p->weight     = p_gc->weight[j];
    p->id         = p_gc->id[j];
    p->cputime    = p_gc->cputime[j];
    p->theta      = p_gc->theta[j];
    p->endcond    = p_gc->endcond[j];
    p->walltile   = p_gc->walltile[j];
    p->mileage    = p_gc->mileage[j];

    /* Guiding center to particle transformation */
    real B_dB[15];
    B_dB[0]       = p_gc->B_r[j];
    B_dB[1]       = p_gc->B_r_dr[j];
    B_dB[2]       = p_gc->B_r_dphi[j];
    B_dB[3]       = p_gc->B_r_dz[j];
    B_dB[4]       = p_gc->B_phi[j];
    B_dB[5]       = p_gc->B_phi_dr[j];
    B_dB[6]       = p_gc->B_phi_dphi[j];
    B_dB[7]       = p_gc->B_phi_dz[j];
    B_dB[8]       = p_gc->B_z[j];
    B_dB[9]       = p_gc->B_z_dr[j];
    B_dB[10]      = p_gc->B_z_dphi[j];
    B_dB[11]      = p_gc->B_z_dz[j];

    real Bnorm = math_normc(p_gc->B_r[j], p_gc->B_phi[j], p_gc->B_z[j]);
    p->ekin       = physlib_Ekin_ppar(p_gc->mass[j], p_gc->mu[j], p_gc->ppar[j], Bnorm);
    p->pitch      = physlib_gc_xi(p_gc->mass[j], p_gc->mu[j], p_gc->ppar[j], Bnorm);
    p->zeta       = p_gc->zeta[j];

    real pr, pphi, pz;
    if(!err) {
        real pparprt, muprt, zetaprt;
        gctransform_guidingcenter2particle(
            p->mass, p->charge, B_dB,
            p->r, p->phi, p->z, p_gc->ppar[j], p_gc->mu[j], p->zeta,
            &p->rprt, &p->phiprt, &p->zprt, &pparprt, &muprt, &zetaprt);

        B_field_eval_B_dB(B_dB, p->rprt, p->phiprt, p->zprt, p->time, bfield);

        gctransform_pparmuzeta2prpphipz(
            p->mass, p->charge, B_dB,
            p->phiprt, pparprt, muprt, zetaprt,
            &pr, &pphi, &pz);
    }
    if(!err && p->rprt <= 0) {err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);}

    if(!err) {
        p->pr   = pr;
        p->pphi = pphi;
        p->pz   = pz;
    }

    /* If marker already has error flag, make sure it is not overwritten here */
    a5err simerr  = p_gc->err[j];
    if(simerr) {
        p->err = simerr;
    }
    else {
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
 * @param p_ml pointer to SIMD structure where marker is being stored
 * @param j index where in the SIMD structure marker is stored
 * @param bfield pointer to magnetic field data
 */
a5err particle_state_to_ml(particle_state* p, int i, particle_simd_ml* p_ml,
                           int j, Bfield *bfield) {
    a5err err = p->err;

    if(!err) {
        p_ml->r[j]          = p->r;
        p_ml->phi[j]        = p->phi;
        p_ml->z[j]          = p->z;

        p_ml->pitch[j]      = 2*(p->pitch >= 0) - 1.0;
        p_ml->time[j]       = p->time;
        p_ml->weight[j]     = p->weight;
        p_ml->id[j]         = p->id;
        p_ml->cputime[j]    = p->cputime;
        p_ml->theta[j]      = p->theta;
        p_ml->endcond[j]    = p->endcond;
        p_ml->walltile[j]   = p->walltile;
        p_ml->mileage[j]    = p->mileage;

        real B_dB[15], psi[1], rho[2];
        if(!err) {
            err = B_field_eval_B_dB(B_dB, p->r, p->phi, p->z, p->time, bfield);
        }
        if(!err) {
            err = B_field_eval_psi(psi, p->r, p->phi, p->z, p->time, bfield);
        }
        if(!err) {
            err = B_field_eval_rho(rho, psi[0], bfield);
        }
        p_ml->rho[j]        = rho[0];
        p_ml->B_r[j]        = B_dB[0];
        p_ml->B_r_dr[j]     = B_dB[1];
        p_ml->B_r_dphi[j]   = B_dB[2];
        p_ml->B_r_dz[j]     = B_dB[3];
        p_ml->B_phi[j]      = B_dB[4];
        p_ml->B_phi_dr[j]   = B_dB[5];
        p_ml->B_phi_dphi[j] = B_dB[6];
        p_ml->B_phi_dz[j]   = B_dB[7];
        p_ml->B_z[j]        = B_dB[8];
        p_ml->B_z_dr[j]     = B_dB[9];
        p_ml->B_z_dphi[j]   = B_dB[10];
        p_ml->B_z_dz[j]     = B_dB[11];

        p_ml->running[j] = 1;
        if(p->endcond) {
            p_ml->running[j] = 0;
        }
        p_ml->index[j] = i;

        p_ml->err[j] = 0;
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
 * @param p_ml pointer to SIMD structure
 * @param j SIMD index of the marker to be converted
 * @param p pointer to state whose fields are used to store the converted marker
 * @param bfield pointer to magnetic field data
 */
void particle_ml_to_state(particle_simd_ml* p_ml, int j, particle_state* p) {
    a5err err = 0;

    p->rprt       = p_ml->r[j];
    p->phiprt     = p_ml->phi[j];
    p->zprt       = p_ml->z[j];
    p->pr         = 0;
    p->pphi       = 0;
    p->pz         = 0;

    p->r          = p_ml->r[j];
    p->phi        = p_ml->phi[j];
    p->z          = p_ml->z[j];
    p->ekin       = 0;
    p->pitch      = p_ml->pitch[j];
    p->zeta       = 0;
    p->mass       = 0;
    p->charge     = 0;
    p->time       = p_ml->time[j];
    p->weight     = p_ml->weight[j];
    p->id         = p_ml->id[j];
    p->cputime    = p_ml->cputime[j];
    p->theta      = p_ml->theta[j];
    p->endcond    = p_ml->endcond[j];
    p->walltile   = p_ml->walltile[j];
    p->mileage    = p_ml->mileage[j];
    p->err        = p_ml->err[j];

    /* If marker already has error flag, make sure it is not overwritten here */
    a5err simerr  = p_ml->err[j];
    if(simerr) {
        p->err = simerr;
    }
    else {
        p->err = err;
    }
    if(!simerr && err) {err = err;}
    p->err = err;
}

/**
 * @brief Convert FO struct into a GC struct
 *
 * @deprecated Do we need this as same can be accomplished by turning
 *             FO to state and state to GC?
 *
 * @param p_fo  fo SIMD structure being transformed
 * @param j     index where in the SIMD structure marker is stored
 * @param p_gc  gc SIMD structure where marker is transformed
 * @param bfield pointer to magnetic field data
 */
int particle_fo_to_gc(particle_simd_fo* p_fo, int j, particle_simd_gc* p_gc,
                      Bfield *bfield) {
    a5err err = p_fo->err[j];
    real axisrz[2];
    int simerr = 0; /* Error has already occurred */
    if(err) {simerr = 1;}
    p_gc->id[j]      = p_fo->id[j];
    p_gc->index[j]   = p_fo->index[j];

    real r, phi, z, ppar, mu, zeta, B_dB[15];
    if(!err) {
        real Rprt   = p_fo->r[j];
        real phiprt = p_fo->phi[j];
        real zprt   = p_fo->z[j];
        real pr     = p_fo->p_r[j];
        real pphi   = p_fo->p_phi[j];
        real pz     = p_fo->p_z[j];
        real mass   = p_fo->mass[j];
        real charge = p_fo->charge[j];

        p_gc->mass[j]     = p_fo->mass[j];
        p_gc->charge[j]   = p_fo->charge[j];
        p_gc->weight[j]   = p_fo->weight[j];
        p_gc->time[j]     = p_fo->time[j];
        p_gc->mileage[j]  = p_fo->mileage[j];
        p_gc->endcond[j]  = p_fo->endcond[j];
        p_gc->running[j]  = p_fo->running[j];
        p_gc->walltile[j] = p_fo->walltile[j];
        p_gc->cputime[j]  = p_fo->cputime[j];

        B_dB[0]       = p_fo->B_r[j];
        B_dB[1]       = p_fo->B_r_dr[j];
        B_dB[2]       = p_fo->B_r_dphi[j];
        B_dB[3]       = p_fo->B_r_dz[j];
        B_dB[4]       = p_fo->B_phi[j];
        B_dB[5]       = p_fo->B_phi_dr[j];
        B_dB[6]       = p_fo->B_phi_dphi[j];
        B_dB[7]       = p_fo->B_phi_dz[j];
        B_dB[8]       = p_fo->B_z[j];
        B_dB[9]       = p_fo->B_z_dr[j];
        B_dB[10]      = p_fo->B_z_dphi[j];
        B_dB[11]      = p_fo->B_z_dz[j];

        /* Guiding center transformation */
        gctransform_particle2guidingcenter(
            mass, charge, B_dB,
            Rprt, phiprt, zprt, pr , pphi, pz,
            &r, &phi, &z, &ppar, &mu, &zeta);
    }
    if(!err && r <= 0)  {err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);}
    if(!err && mu < 0)  {err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);}

    real psi[1], rho[2];
    if(!err) {
        err = B_field_eval_B_dB(
            B_dB, r, phi, z, p_fo->time[j], bfield);
    }
    if(!err) {
        err = B_field_eval_psi(
            psi, r, phi, z, p_fo->time[j], bfield);
    }
    if(!err) {
        err = B_field_eval_rho(rho, psi[0], bfield);
    }
    if(!err) {
        err = B_field_get_axis_rz(axisrz, bfield, p_gc->phi[j]);
    }

    if(!err) {
        p_gc->r[j]    = r;
        p_gc->phi[j]  = phi;
        p_gc->z[j]    = z;
        p_gc->mu[j]   = mu;
        p_gc->zeta[j] = zeta;
        p_gc->ppar[j] = ppar;
        p_gc->rho[j]  = rho[0];

        /* Evaluate pol angle so that it is cumulative and at gc position */
        p_gc->theta[j]  = p_fo->theta[j];
        p_gc->theta[j] += atan2( (p_fo->r[j]-axisrz[0]) * (p_gc->z[j]-axisrz[1])
                               - (p_fo->z[j]-axisrz[1]) * (p_gc->r[j]-axisrz[0]),
                                 (p_fo->r[j]-axisrz[0]) * (p_gc->r[j]-axisrz[0])
                               + (p_fo->z[j]-axisrz[1]) * (p_gc->z[j]-axisrz[1]) );

        p_gc->B_r[j]        = B_dB[0];
        p_gc->B_r_dr[j]     = B_dB[1];
        p_gc->B_r_dphi[j]   = B_dB[2];
        p_gc->B_r_dz[j]     = B_dB[3];

        p_gc->B_phi[j]      = B_dB[4];
        p_gc->B_phi_dr[j]   = B_dB[5];
        p_gc->B_phi_dphi[j] = B_dB[6];
        p_gc->B_phi_dz[j]   = B_dB[7];

        p_gc->B_z[j]        = B_dB[8];
        p_gc->B_z_dr[j]     = B_dB[9];
        p_gc->B_z_dphi[j]   = B_dB[10];
        p_gc->B_z_dz[j]     = B_dB[11];
    }
    if(!simerr) {err = err;}
    p_gc->err[j] = err;
    if(p_gc->err[j]) {
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
void particle_copy_fo(particle_simd_fo* p1, int i, particle_simd_fo* p2, int j) {
        p2->r[j]          = p1->r[i];
        p2->phi[j]        = p1->phi[i];
        p2->z[j]          = p1->z[i];
        p2->p_r[j]        = p1->p_r[i];
        p2->p_phi[j]      = p1->p_phi[i];
        p2->p_z[j]        = p1->p_z[i];

        p2->time[j]       = p1->time[i];
        p2->mileage[j]    = p1->mileage[i];
        p2->cputime[j]    = p1->cputime[i];
        p2->rho[j]        = p1->rho[i];
        p2->weight[j]     = p1->weight[i];
        p2->cputime[j]    = p1->cputime[i];
        p2->rho[j]        = p1->rho[i];
        p2->theta[j]      = p1->theta[i];

        p2->mass[j]       = p1->mass[i];
        p2->charge[j]     = p1->charge[i];
        p2->znum[j]       = p1->znum[i];
        p2->anum[j]       = p1->anum[i];

        p2->id[j]         = p1->id[i];
        p2->bounces[j]    = p1->bounces[i];
        p2->running[j]    = p1->running[i];
        p2->endcond[j]    = p1->endcond[i];
        p2->walltile[j]   = p1->walltile[i];

        p2->B_r[j]        = p1->B_r[i];
        p2->B_phi[j]      = p1->B_phi[i];
        p2->B_z[j]        = p1->B_z[i];

        p2->B_r_dr[j]     = p1->B_r_dr[i];
        p2->B_r_dphi[j]   = p1->B_r_dphi[i];
        p2->B_r_dz[j]     = p1->B_r_dz[i];

        p2->B_phi_dr[j]   = p1->B_phi_dr[i];
        p2->B_phi_dphi[j] = p1->B_phi_dphi[i];
        p2->B_phi_dz[j]   = p1->B_phi_dz[i];

        p2->B_z_dr[j]     = p1->B_z_dr[i];
        p2->B_z_dphi[j]   = p1->B_z_dphi[i];
        p2->B_z_dz[j]     = p1->B_z_dz[i];
}

/**
 * @brief Copy GC struct
 *
 * @param p1 SIMD structure for input
 * @param i  index for the copied input
 * @param p2 SIMD structure for output
 * @param j  index for the output slot
 */
void particle_copy_gc(particle_simd_gc* p1, int i, particle_simd_gc* p2, int j) {
    p2->r[j]          = p1->r[i];
    p2->phi[j]        = p1->phi[i];
    p2->z[j]          = p1->z[i];
    p2->ppar[j]       = p1->ppar[i];
    p2->mu[j]         = p1->mu[i];
    p2->zeta[j]       = p1->zeta[i];

    p2->time[j]       = p1->time[i];
    p2->mileage[j]    = p1->mileage[i];
    p2->weight[j]     = p1->weight[i];
    p2->cputime[j]    = p1->cputime[i];
    p2->rho[j]        = p1->rho[i];
    p2->theta[j]      = p1->theta[i];

    p2->mass[j]       = p1->mass[i];
    p2->charge[j]     = p1->charge[i];

    p2->id[j]         = p1->id[i];
    p2->bounces[j]    = p1->bounces[i];
    p2->running[j]    = p1->running[i];
    p2->endcond[j]    = p1->endcond[i];
    p2->walltile[j]   = p1->walltile[i];

    p2->B_r[j]        = p1->B_r[i];
    p2->B_phi[j]      = p1->B_phi[i];
    p2->B_z[j]        = p1->B_z[i];

    p2->B_r_dr[j]     = p1->B_r_dr[i];
    p2->B_r_dphi[j]   = p1->B_r_dphi[i];
    p2->B_r_dz[j]     = p1->B_r_dz[i];

    p2->B_phi_dr[j]   = p1->B_phi_dr[i];
    p2->B_phi_dphi[j] = p1->B_phi_dphi[i];
    p2->B_phi_dz[j]   = p1->B_phi_dz[i];

    p2->B_z_dr[j]     = p1->B_z_dr[i];
    p2->B_z_dphi[j]   = p1->B_z_dphi[i];
    p2->B_z_dz[j]     = p1->B_z_dz[i];
}

/**
 * @brief Copy ML struct
 *
 * @param p1 SIMD structure for input
 * @param i  index for the copied input
 * @param p2 SIMD structure for output
 * @param j  index for the output slot
 */
void particle_copy_ml(particle_simd_ml* p1, int i, particle_simd_ml* p2,
                      int j) {
    p2->r[j]          = p1->r[i];
    p2->phi[j]        = p1->phi[i];
    p2->z[j]          = p1->z[i];
    p2->pitch[j]      = p1->pitch[i];

    p2->time[j]       = p1->time[i];
    p2->mileage[j]    = p1->mileage[i];
    p2->cputime[j]    = p1->cputime[i];
    p2->rho[j]        = p1->rho[i];
    p2->weight[j]     = p1->weight[i];
    p2->theta[j]      = p1->theta[i];

    p2->id[j]         = p1->id[i];
    p2->running[j]    = p1->running[i];
    p2->endcond[j]    = p1->endcond[i];
    p2->walltile[j]   = p1->walltile[i];

    p2->B_r[j]        = p1->B_r[i];
    p2->B_phi[j]      = p1->B_phi[i];
    p2->B_z[j]        = p1->B_z[i];

    p2->B_r_dr[j]     = p1->B_r_dr[i];
    p2->B_r_dphi[j]   = p1->B_r_dphi[i];
    p2->B_r_dz[j]     = p1->B_r_dz[i];

    p2->B_phi_dr[j]   = p1->B_phi_dr[i];
    p2->B_phi_dphi[j] = p1->B_phi_dphi[i];
    p2->B_phi_dz[j]   = p1->B_phi_dz[i];

    p2->B_z_dr[j]     = p1->B_z_dr[i];
    p2->B_z_dphi[j]   = p1->B_z_dphi[i];
    p2->B_z_dz[j]     = p1->B_z_dz[i];
}

/**
 * @brief Offload particle struct to GPU.
 *
 * @param p pointer to the particle struct to be offloaded.
 */
void particle_offload_fo(particle_simd_fo* p) {
    SUPPRESS_UNUSED_WARNING(p);
    GPU_MAP_TO_DEVICE(
        p[0:1],\
        p->running   [0:p->n_mrk],\
        p->r         [0:p->n_mrk],\
        p->phi       [0:p->n_mrk],\
        p->p_r       [0:p->n_mrk],\
        p->p_phi     [0:p->n_mrk],\
        p->p_z       [0:p->n_mrk],\
        p->mileage   [0:p->n_mrk],\
        p->z         [0:p->n_mrk],\
        p->charge    [0:p->n_mrk],\
        p->mass      [0:p->n_mrk],\
        p->B_r       [0:p->n_mrk],\
        p->B_r_dr    [0:p->n_mrk],\
        p->B_r_dphi  [0:p->n_mrk],\
        p->B_r_dz    [0:p->n_mrk],\
        p->B_phi     [0:p->n_mrk],\
        p->B_phi_dr  [0:p->n_mrk],\
        p->B_phi_dphi[0:p->n_mrk],\
        p->B_phi_dz  [0:p->n_mrk],\
        p->B_z       [0:p->n_mrk],\
        p->B_z_dr    [0:p->n_mrk],\
        p->B_z_dphi  [0:p->n_mrk],\
        p->B_z_dz    [0:p->n_mrk],\
        p->rho       [0:p->n_mrk],\
        p->theta     [0:p->n_mrk],\
        p->err       [0:p->n_mrk],\
        p->time      [0:p->n_mrk],\
        p->weight    [0:p->n_mrk],\
        p->cputime   [0:p->n_mrk],\
        p->id        [0:p->n_mrk],\
        p->endcond   [0:p->n_mrk],\
        p->walltile  [0:p->n_mrk],\
        p->index     [0:p->n_mrk],\
        p->znum      [0:p->n_mrk],\
        p->anum      [0:p->n_mrk],\
        p->bounces   [0:p->n_mrk]
    )
}

/**
 * @brief Onload particle struct from the GPU.
 *
 * @param p pointer to the particle struct to be onloaded.
 */
void particle_onload_fo(particle_simd_fo* p) {
    SUPPRESS_UNUSED_WARNING(p);
    GPU_UPDATE_FROM_DEVICE(
        p->running   [0:p->n_mrk],\
        p->r         [0:p->n_mrk],\
        p->phi       [0:p->n_mrk],\
        p->p_r       [0:p->n_mrk],\
        p->p_phi     [0:p->n_mrk],\
        p->p_z       [0:p->n_mrk],\
        p->mileage   [0:p->n_mrk],\
        p->z         [0:p->n_mrk],\
        p->charge    [0:p->n_mrk],\
        p->mass      [0:p->n_mrk],\
        p->B_r       [0:p->n_mrk],\
        p->B_r_dr    [0:p->n_mrk],\
        p->B_r_dphi  [0:p->n_mrk],\
        p->B_r_dz    [0:p->n_mrk],\
        p->B_phi     [0:p->n_mrk],\
        p->B_phi_dr  [0:p->n_mrk],\
        p->B_phi_dphi[0:p->n_mrk],\
        p->B_phi_dz  [0:p->n_mrk],\
        p->B_z       [0:p->n_mrk],\
        p->B_z_dr    [0:p->n_mrk],\
        p->B_z_dphi  [0:p->n_mrk],\
        p->B_z_dz    [0:p->n_mrk],\
        p->rho       [0:p->n_mrk],\
        p->theta     [0:p->n_mrk],\
        p->err       [0:p->n_mrk],\
        p->time      [0:p->n_mrk],\
        p->weight    [0:p->n_mrk],\
        p->cputime   [0:p->n_mrk],\
        p->id        [0:p->n_mrk],\
        p->endcond   [0:p->n_mrk],\
        p->walltile  [0:p->n_mrk],\
        p->index     [0:p->n_mrk],\
        p->znum      [0:p->n_mrk],\
        p->anum      [0:p->n_mrk],\
        p->bounces   [0:p->n_mrk]
    )
}
