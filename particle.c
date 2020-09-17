/**
 * @file particle.c
 * @brief Marker structs and conversions between them
 *
 * ASCOT5 is a code that simulates markers.
 *
 * This file contains functions to generate dummy markers, fetching markers
 * from the marker queue, and handles conversions between various marker
 * structs.
 *
 * The relationship between the seven different marker structs is:
 *
 *    particle >--+            +--> particle_simd_fo
 *                |            |
 * particle_gc >--particle_state--> particle_simd_gc
 *                |            |
 * particle_ml >--+            +--> particle_simd_ml
 *
 * The structs on the left are input structs which are constructed when marker
 * data is read. These structs must be independend of any other input data. A
 * consistent representation of markers is provided by particle_state which
 * can be constructed from any type of input marker.
 *
 * Particle state can be used to obtain any information from marker physical
 * quantities such as guiding center energy, which requires one knows both
 * guiding center magnetic moment and local magnetic field value. All markers
 * are represented by particle state when simulation begins and again when
 * simulation ends.
 *
 * During simulation, markers are represented by a SIMD compatible structure
 * (all fields in these structs are NSIMD long arrays representing NSIMD
 * markers), and the type of SIMD structure depends on the simulation type:
 * whether simulation models markers as particles, guiding centers or field
 * lines. Each simulation struct can be constructed from particle state.
 *
 * So in principle, one can have any type of marker input and still use any kind
 * of simulation method, e.g. input can be guiding centers but these can be
 * modelled as particles, or to be precise, the particles are modelled whose
 * guiding centers are given as inputs.
 *
 * Exception to this rule: field line input can only be used to simulate field
 * lines as field lines have no mass.
 *
 * If one is to implement a new simulation mode, e.g. particle simulation in
 * Boozer coordinates, one is advised to implement a new marker struct for it.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ascot5.h"
#include "error.h"
#include "consts.h"
#include "math.h"
#include "physlib.h"
#include "gctransform.h"
#include "particle.h"
#include "B_field.h"
#include "E_field.h"


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
    p_fo->rdot[j]     = 1;
    p_fo->phidot[j]   = 1;
    p_fo->zdot[j]     = 1;
    p_fo->mass[j]     = 1;
    p_fo->charge[j]   = 1;
    p_fo->weight[j]   = 0;
    p_fo->time[j]     = 0;
    p_fo->id[j]       = -1;
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
    p_gc->vpar[j]       = 1;
    p_gc->mu[j]         = 1;
    p_gc->zeta[j]       = 1;
    p_gc->mass[j]       = 1;
    p_gc->charge[j]     = 1;
    p_gc->time[j]       = 0;
    p_gc->weight[j]     = 0;
    p_gc->id[j]         = -1;
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
 * @param Bdata pointer to magnetic field data
 * @param cycle pointer to integer array where what was done for each marker
 *              is stored
 *
 * @return Number of markers within the SIMD structure that are still running
 */
int particle_cycle_fo(particle_queue* q, particle_simd_fo* p,
                      B_field_data* Bdata, int* cycle) {

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
            particle_fo_to_state(p, i, q->p[p->index[i]], Bdata);
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
                a5err err = particle_state_to_fo(
                    q->p[i_prt], i_prt, p, i, Bdata);
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
 * @param Bdata pointer to magnetic field data
 * @param cycle pointer to integer array where what was done for each marker
 *              is stored
 *
 * @return Number of markers within the SIMD structure that are still running
 */
int particle_cycle_gc(particle_queue* q, particle_simd_gc* p,
                      B_field_data* Bdata, int* cycle) {

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
            particle_gc_to_state(p, i, q->p[p->index[i]], Bdata);
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
                    q->p[i_prt], i_prt, p, i, Bdata);
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
 * @param Bdata pointer to magnetic field data
 * @param cycle pointer to integer array where what was done for each marker
 *              is stored
 *
 * @return Number of markers within the SIMD structure that are still running
 */
int particle_cycle_ml(particle_queue* q, particle_simd_ml* p,
                      B_field_data* Bdata, int* cycle) {

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
            particle_ml_to_state(p, i, q->p[p->index[i]], Bdata);
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
                    q->p[i_prt], i_prt, p, i, Bdata);
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
 * @brief Converts input marker to a marker state
 *
 * When marker input is read, the data is stored in one of the input marker
 * structs (particle, particle_gc, particle_ml). This function converts
 * these structs into a marker state, which then can be used to create
 * simulation structs.
 *
 * If an error is encountered while converting, the error message is stored
 * to the state.
 *
 * @todo Make separate functions for each different cases (that are called from
 *       here) to make this function cleaner
 * @todo Parameter p can be a particle_state itself but this option is not dealt
 *       with
 * @todo This sets p->type = input_particle_type_s which is not in accordance
 *       with the documentation.
 *
 * @param p pointer to marker input
 * @param ps pointer to state where converted marker will be stored
 * @param Bdata pointer to magnetic field data
 */
void particle_input_to_state(input_particle* p, particle_state* ps,
                             B_field_data* Bdata) {
    a5err err = 0;
    integer id;

    if(p->type == input_particle_type_p) {
        /* Check that input is valid */
        if(!err && p->p.r <= 0) {
            err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);
        }
        if(!err && ( math_normc(p->p.v_r,p->p.v_phi,p->p.v_z) >= CONST_C2 )) {
            err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);
        }
        if(!err && p->p.mass <= 0) {
            err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);
        }
        if(!err && p->p.weight <= 0) {
            err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);
        }
        if(!err && p->p.id <= 0) {
            err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);
        }

        /* Particle to state */
        p->type = input_particle_type_s;
        id = p->p.id;

        if(!err) {
            ps->rprt   = p->p.r;
            ps->phiprt = p->p.phi;
            ps->zprt   = p->p.z;
            ps->rdot   = p->p.v_r;
            ps->phidot = p->p.v_phi/ps->rprt;
            ps->zdot   = p->p.v_z;
            ps->mass   = p->p.mass;
            ps->charge = p->p.charge;
            ps->anum   = p->p.anum;
            ps->znum   = p->p.znum;
            ps->weight = p->p.weight;
            ps->time   = p->p.time;
            ps->theta  = atan2(ps->zprt-B_field_get_axis_r(Bdata, ps->phiprt),
                               ps->rprt-B_field_get_axis_z(Bdata, ps->phiprt));
            ps->id       = id;
            ps->endcond  = 0;
            ps->walltile = 0;
            ps->cputime  = 0;
        }

        /* Guiding center transformation */
        real B_dB[15], r, phi, z, vpar, mu, zeta, psi[1], rho[1];
        if(!err) {
            err = B_field_eval_B_dB(B_dB, ps->rprt, ps->phiprt, ps->zprt,
                                    ps->time, Bdata);
        }

        if(!err) {
            gctransform_particle2guidingcenter(
                ps->mass, ps->charge, B_dB,
                ps->rprt, ps->phiprt, ps->zprt, p->p.v_r, p->p.v_phi, p->p.v_z,
                &r, &phi, &z, &vpar, &mu, &zeta);
        }
        if(!err && r <= 0)  {err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);}
        if(!err && mu < 0)  {err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);}

        /* Update magnetic field at gc position */
        if(!err) {
            err = B_field_eval_B_dB(B_dB, r, phi, z, ps->time, Bdata);
        }
        if(!err) {
            err = B_field_eval_psi(psi, r, phi, z, ps->time, Bdata);
        }
        if(!err) {
            err = B_field_eval_rho(rho, psi[0], Bdata);
        }
        if(!err && vpar >= CONST_C) {err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);}

        if(!err) {
            ps->r     = r;
            ps->phi   = phi;
            ps->z     = z;
            ps->vpar  = vpar;
            ps->mu    = mu;
            ps->zeta  = zeta;

            ps->rho        = rho[0];
            ps->B_r        = B_dB[0];
            ps->B_phi      = B_dB[4];
            ps->B_z        = B_dB[8];
            ps->B_r_dr     = B_dB[1];
            ps->B_phi_dr   = B_dB[5];
            ps->B_z_dr     = B_dB[9];
            ps->B_r_dphi   = B_dB[2];
            ps->B_phi_dphi = B_dB[6];
            ps->B_z_dphi   = B_dB[10];
            ps->B_r_dz     = B_dB[3];
            ps->B_phi_dz   = B_dB[7];
            ps->B_z_dz     = B_dB[11];

            ps->err = 0;
        }

    }
    else if(p->type == input_particle_type_gc) {
        /* Check that input is valid */
        if(!err && p->p_gc.r <= 0)          {err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);}
        if(!err && fabs(p->p_gc.pitch) > 1) {err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);}
        if(!err && p->p_gc.energy <= 0)     {err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);}
        if(!err && p->p_gc.mass <= 0)       {err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);}
        if(!err && p->p_gc.weight <= 0)     {err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);}
        if(!err && p->p_gc.id <= 0)         {err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);}

        /* Guiding center to state */
        p->type = input_particle_type_s;
        id = p->p_gc.id;

        real B_dB[15], psi[1], rho[1];
        if(!err) {
            err = B_field_eval_B_dB(B_dB, p->p_gc.r, p->p_gc.phi, p->p_gc.z,
                                    p->p_gc.time, Bdata);
        }
        if(!err) {
            err = B_field_eval_psi(psi, p->p_gc.r, p->p_gc.phi, p->p_gc.z,
                                   p->p_gc.time, Bdata);
        }
        if(!err) {
            err = B_field_eval_rho(rho, psi[0], Bdata);
        }

        if(!err) {
            ps->rho        = rho[0];
            ps->B_r        = B_dB[0];
            ps->B_phi      = B_dB[4];
            ps->B_z        = B_dB[8];
            ps->B_r_dr     = B_dB[1];
            ps->B_phi_dr   = B_dB[5];
            ps->B_z_dr     = B_dB[9];
            ps->B_r_dphi   = B_dB[2];
            ps->B_phi_dphi = B_dB[6];
            ps->B_z_dphi   = B_dB[10];
            ps->B_r_dz     = B_dB[3];
            ps->B_phi_dz   = B_dB[7];
            ps->B_z_dz     = B_dB[11];
        }

        /* Input is in (Ekin,xi) coordinates but state needs (mu,vpar) so we need to do that
         * transformation first. */
        real gamma, mu, vpar;
        if(!err) {
            /* From kinetic energy we get Lorentz factor as gamma = 1 + Ekin/mc^2 */
            gamma = 1 + p->p_gc.energy / (p->p_gc.mass * CONST_C2);

            /* And then we can use the usual formula for Lorentz factor to get total velocity */
            real v = sqrt(1 - 1.0 / (gamma * gamma)) * CONST_C;

            /* Now we can use library functions for transformation */
            real Bnorm = math_normc(B_dB[0], B_dB[4], B_dB[8]);
            vpar = physlib_gc_vpar(v, p->p_gc.pitch);
            mu   = physlib_gc_mu(p->p_gc.mass, v, p->p_gc.pitch, Bnorm);
        }
        if(!err && mu < 0)          {err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);}
        if(!err && vpar >= CONST_C) {err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);}

        if(!err) {
            ps->r        = p->p_gc.r;
            ps->phi      = p->p_gc.phi;
            ps->z        = p->p_gc.z;
            ps->mu       = mu;
            ps->vpar     = vpar;
            ps->zeta     = p->p_gc.zeta;
            ps->mass     = p->p_gc.mass;
            ps->charge   = p->p_gc.charge;
            ps->anum     = p->p_gc.anum;
            ps->znum     = p->p_gc.znum;
            ps->weight   = p->p_gc.weight;
            ps->time     = p->p_gc.time;
            ps->theta    = atan2(ps->z-B_field_get_axis_z(Bdata, ps->phi),
                                 ps->r-B_field_get_axis_r(Bdata, ps->phi));
            ps->id       = id;
            ps->endcond  = 0;
            ps->walltile = 0;
            ps->cputime  = 0;
        }

        /* Guiding center transformation to get particle coordinates */
        real rprt, phiprt, zprt, vR, vphi, vz;
        if(!err) {
            real vparprt, muprt, zetaprt;
            gctransform_guidingcenter2particle(
                ps->mass, ps->charge, B_dB,
                ps->r, ps->phi, ps->z, ps->vpar, ps->mu, ps->zeta,
                &rprt, &phiprt, &zprt, &vparprt, &muprt, &zetaprt);

            B_field_eval_B_dB(B_dB, rprt, phiprt, zprt, ps->time, Bdata);

            gctransform_vparmuzeta2vRvphivz(
                ps->mass, ps->charge, B_dB,
                phiprt, vparprt, muprt, zetaprt,
                &vR, &vphi, &vz);
        }
        if(!err && rprt <= 0) {err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);}

        if(!err) {
            ps->rprt   = rprt;
            ps->phiprt = phiprt;
            ps->zprt   = zprt;
            ps->rdot   = vR;
            ps->phidot = vphi / ps->rprt;
            ps->zdot   = vz;

            ps->err = 0;
        }
    }
    else if(p->type == input_particle_type_ml) {
        /* Check that input is valid */
        if(!err && p->p_ml.r <= 0)      {err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);}
        if(!err && p->p_ml.weight <= 0) {err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);}
        if(!err && p->p_ml.id <= 0)     {err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);} 

        /* Magnetic field line to state */
        p->type = input_particle_type_s;
        id = p->p_ml.id;

        real B_dB[15], psi[1], rho[1];
        if(!err) {
            err = B_field_eval_B_dB(B_dB, p->p_ml.r, p->p_ml.phi, p->p_ml.z,
                                    p->p_ml.time, Bdata);
        }
        if(!err) {
            err = B_field_eval_psi(psi, p->p_ml.r, p->p_ml.phi, p->p_ml.z,
                                   p->p_ml.time, Bdata);
        }
        if(!err) {
            err = B_field_eval_rho(rho, psi[0], Bdata);
        }

        if(!err) {
            ps->rprt       = p->p_ml.r;
            ps->phiprt     = p->p_ml.phi;
            ps->zprt       = p->p_ml.z;
            ps->rdot       = 0;
            ps->phidot     = 0;
            ps->zdot       = 0;

            ps->mass       = 0;
            ps->charge     = 0;
            ps->anum       = 0;
            ps->znum       = 0;
            ps->weight     = p->p_ml.weight;
            ps->time       = p->p_ml.time;
            ps->id         = id;
            ps->theta      = atan2(p->p_ml.z - B_field_get_axis_z(Bdata, ps->phiprt),
                                   p->p_ml.r - B_field_get_axis_r(Bdata, ps->phiprt));
            ps->endcond    = 0;
            ps->walltile   = 0;
            ps->cputime    = 0;

            ps->r          = p->p_ml.r;
            ps->phi        = p->p_ml.phi;
            ps->z          = p->p_ml.z;
            ps->vpar       = p->p_ml.pitch >= 0;
            ps->mu         = 0;
            ps->zeta       = 0;

            ps->rho        = rho[0];
            ps->B_r        = B_dB[0];
            ps->B_phi      = B_dB[4];
            ps->B_z        = B_dB[8];
            ps->B_r_dr     = B_dB[1];
            ps->B_phi_dr   = B_dB[5];
            ps->B_z_dr     = B_dB[9];
            ps->B_r_dphi   = B_dB[2];
            ps->B_phi_dphi = B_dB[6];
            ps->B_z_dphi   = B_dB[10];
            ps->B_r_dz     = B_dB[3];
            ps->B_phi_dz   = B_dB[7];
            ps->B_z_dz     = B_dB[11];

            ps->err = 0;
        }
    }

    /* If particle was rejected, ensure that these fields are filled */
    if(err) {
        ps->id      = id;
        ps->endcond = 0;
        ps->err     = err;
    }
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
 * @param Bdata pointer to magnetic field data
 */
a5err particle_state_to_fo(particle_state* p, int i, particle_simd_fo* p_fo,
                           int j, B_field_data* Bdata) {
    a5err err = p->err;

    if(!err) {
        p_fo->r[j]          = p->rprt;
        p_fo->phi[j]        = p->phiprt;
        p_fo->z[j]          = p->zprt;
        p_fo->rdot[j]       = p->rdot;
        p_fo->phidot[j]     = p->phidot;
        p_fo->zdot[j]       = p->zdot;

        p_fo->mass[j]       = p->mass;
        p_fo->charge[j]     = p->charge;
        p_fo->weight[j]     = p->weight;
        p_fo->time[j]       = p->time;
        p_fo->theta[j]      = p->theta;
        p_fo->id[j]         = p->id;
        p_fo->endcond[j]    = p->endcond;
        p_fo->walltile[j]   = p->walltile;
    }

    /* Magnetic field stored in state is for the gc position */
    real B_dB[15], psi[1], rho[1];
    if(!err) {
        err = B_field_eval_B_dB(B_dB, p->rprt, p->phiprt, p->zprt, p->time,
                                Bdata);
    }
    if(!err) {
        err = B_field_eval_psi(psi, p->rprt, p->phiprt, p->zprt, p->time,
                               Bdata);
    }
    if(!err) {
        err = B_field_eval_rho(rho, psi[0], Bdata);
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
 * @param Bdata pointer to magnetic field data
 */
void particle_fo_to_state(particle_simd_fo* p_fo, int j, particle_state* p,
                          B_field_data* Bdata) {
    a5err err = 0;

    p->rprt       = p_fo->r[j];
    p->phiprt     = p_fo->phi[j];
    p->zprt       = p_fo->z[j];
    p->rdot       = p_fo->rdot[j];
    p->phidot     = p_fo->phidot[j];
    p->zdot       = p_fo->zdot[j];

    p->mass       = p_fo->mass[j];
    p->charge     = p_fo->charge[j];
    p->weight     = p_fo->weight[j];
    p->time       = p_fo->time[j];
    p->theta      = p_fo->theta[j];
    p->id         = p_fo->id[j];
    p->endcond    = p_fo->endcond[j];
    p->walltile   = p_fo->walltile[j];
    p->cputime    = p_fo->cputime[j];

    /* Particle to guiding center */
    real B_dB[15], psi[1], rho[1];
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
    real vpar;
    if(!err) {
        real vR   = p->rdot;
        real vphi = p->phidot * p->rprt;
        real vz   = p->zdot;

        gctransform_particle2guidingcenter(
            p->mass, p->charge, B_dB,
            p->rprt, p->phiprt, p->zprt, vR , vphi, vz,
            &p->r, &p->phi, &p->z, &vpar, &p->mu, &p->zeta);
    }
    if(!err && p->r <= 0)  {err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);}
    if(!err && p->mu < 0)  {err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);}

    if(!err) {
        err = B_field_eval_B_dB(B_dB, p->r, p->phi, p->z, p->time, Bdata);
    }
    if(!err) {
        err = B_field_eval_psi(psi, p->r, p->phi, p->z, p->time, Bdata);
    }
    if(!err) {
        err = B_field_eval_rho(rho, psi[0], Bdata);
    }

    if(!err) {
        p->vpar = vpar;
    }
    if(!err && p->vpar >= CONST_C) {
        err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);
    }

    /* Normally magnetic field data at gc position is stored here
     * but, if gc transformation fails, field at particle position is
     * stored instead. */
    p->rho        = rho[0];

    p->B_r        = B_dB[0];
    p->B_r_dr     = B_dB[1];
    p->B_r_dphi   = B_dB[2];
    p->B_r_dz     = B_dB[3];

    p->B_phi      = B_dB[4];
    p->B_phi_dr   = B_dB[5];
    p->B_phi_dphi = B_dB[6];
    p->B_phi_dz   = B_dB[7];

    p->B_z        = B_dB[8];
    p->B_z_dr     = B_dB[9];
    p->B_z_dphi   = B_dB[10];
    p->B_z_dz     = B_dB[11];

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
 * @param Bdata pointer to magnetic field data
 */
a5err particle_state_to_gc(particle_state* p, int i, particle_simd_gc* p_gc,
                           int j, B_field_data* Bdata) {
    a5err err = p->err;

    if(!err) {
        p_gc->hermite_weights[j]          = p->hermite_weights;
        p_gc->hermite_knots[j]          = p->hermite_knots;
        p_gc->use_hermite[j]          = p->use_hermite;
        p_gc->r[j]          = p->r;
        p_gc->phi[j]        = p->phi;
        p_gc->z[j]          = p->z;
        p_gc->vpar[j]       = p->vpar;
        p_gc->mu[j]         = p->mu;
        p_gc->zeta[j]       = p->zeta;

        p_gc->mass[j]       = p->mass;
        p_gc->charge[j]     = p->charge;
        p_gc->time[j]       = p->time;
        p_gc->weight[j]     = p->weight;
        p_gc->rho[j]        = p->rho;
        p_gc->theta[j]      = p->theta;
        p_gc->id[j]         = p->id;
        p_gc->endcond[j]    = p->endcond;
        p_gc->walltile[j]   = p->walltile;

        p_gc->B_r[j]        = p->B_r;
        p_gc->B_r_dr[j]     = p->B_r_dr;
        p_gc->B_r_dphi[j]   = p->B_r_dphi;
        p_gc->B_r_dz[j]     = p->B_r_dz;

        p_gc->B_phi[j]      = p->B_phi;
        p_gc->B_phi_dr[j]   = p->B_phi_dr;
        p_gc->B_phi_dphi[j] = p->B_phi_dphi;
        p_gc->B_phi_dz[j]   = p->B_phi_dz;

        p_gc->B_z[j]        = p->B_z;
        p_gc->B_z_dr[j]     = p->B_z_dr;
        p_gc->B_z_dphi[j]   = p->B_z_dphi;
        p_gc->B_z_dz[j]     = p->B_z_dz;

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
 * @param Bdata pointer to magnetic field data
 */
void particle_gc_to_state(particle_simd_gc* p_gc, int j, particle_state* p,
                          B_field_data* Bdata) {
    a5err err = 0;

    p->hermite_weights          = p_gc->hermite_weights[j];
    p->hermite_knots          = p_gc->hermite_knots[j];
    p->use_hermite          = p_gc->use_hermite[j];
    p->r          = p_gc->r[j];
    p->phi        = p_gc->phi[j];
    p->z          = p_gc->z[j];
    p->vpar       = p_gc->vpar[j];
    p->mu         = p_gc->mu[j];
    p->zeta       = p_gc->zeta[j];

    p->mass       = p_gc->mass[j];
    p->charge     = p_gc->charge[j];
    p->time       = p_gc->time[j];
    p->weight     = p_gc->weight[j];
    p->id         = p_gc->id[j];
    p->cputime    = p_gc->cputime[j];
    p->rho        = p_gc->rho[j];
    p->theta      = p_gc->theta[j];
    p->endcond    = p_gc->endcond[j];
    p->walltile   = p_gc->walltile[j];

    p->B_r        = p_gc->B_r[j];
    p->B_r_dr     = p_gc->B_r_dr[j];
    p->B_r_dphi   = p_gc->B_r_dphi[j];
    p->B_r_dz     = p_gc->B_r_dz[j];

    p->B_phi      = p_gc->B_phi[j];
    p->B_phi_dr   = p_gc->B_phi_dr[j];
    p->B_phi_dphi = p_gc->B_phi_dphi[j];
    p->B_phi_dz   = p_gc->B_phi_dz[j];

    p->B_z        = p_gc->B_z[j];
    p->B_z_dr     = p_gc->B_z_dr[j];
    p->B_z_dphi   = p_gc->B_z_dphi[j];
    p->B_z_dz     = p_gc->B_z_dz[j];

    /* Guiding center to particle transformation */
    real B_dB[15];
    B_dB[0]       = p->B_r;
    B_dB[1]       = p->B_r_dr;
    B_dB[2]       = p->B_r_dphi;
    B_dB[3]       = p->B_r_dz;
    B_dB[4]       = p->B_phi;
    B_dB[5]       = p->B_phi_dr;
    B_dB[6]       = p->B_phi_dphi;
    B_dB[7]       = p->B_phi_dz;
    B_dB[8]       = p->B_z;
    B_dB[9]       = p->B_z_dr;
    B_dB[10]      = p->B_z_dphi;
    B_dB[11]      = p->B_z_dz;

    real vR, vphi, vz;
    if(!err) {
        real vparprt, muprt, zetaprt;
        gctransform_guidingcenter2particle(
            p->mass, p->charge, B_dB,
            p->r, p->phi, p->z, p->vpar, p->mu, p->zeta,
            &p->rprt, &p->phiprt, &p->zprt, &vparprt, &muprt, &zetaprt);

        B_field_eval_B_dB(B_dB, p->rprt, p->phiprt, p->zprt, p->time, Bdata);

        gctransform_vparmuzeta2vRvphivz(
            p->mass, p->charge, B_dB,
            p->phiprt, vparprt, muprt, zetaprt,
            &vR, &vphi, &vz);
    }
    if(!err && p->rprt <= 0) {err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);}

    if(!err) {
        p->rdot       = vR;
        p->phidot     = vphi/p->rprt;
        p->zdot       = vz;
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
 * @param Bdata pointer to magnetic field data
 */
a5err particle_state_to_ml(particle_state* p, int i, particle_simd_ml* p_ml,
                           int j, B_field_data* Bdata) {
    a5err err = p->err;

    if(!err) {
        p_ml->r[j]          = p->r;
        p_ml->phi[j]        = p->phi;
        p_ml->z[j]          = p->z;

        p_ml->pitch[j]      = 2*(p->vpar >= 0) - 1.0;
        p_ml->time[j]       = p->time;
        p_ml->weight[j]     = p->weight;
        p_ml->id[j]         = p->id;
        p_ml->cputime[j]    = p->cputime;
        p_ml->rho[j]        = p->rho;
        p_ml->theta[j]      = p->theta;
        p_ml->endcond[j]    = p->endcond;
        p_ml->walltile[j]   = p->walltile;

        p_ml->B_r[j]        = p->B_r;
        p_ml->B_r_dr[j]     = p->B_r_dr;
        p_ml->B_r_dphi[j]   = p->B_r_dphi;
        p_ml->B_r_dz[j]     = p->B_r_dz;

        p_ml->B_phi[j]      = p->B_phi;
        p_ml->B_phi_dr[j]   = p->B_phi_dr;
        p_ml->B_phi_dphi[j] = p->B_phi_dphi;
        p_ml->B_phi_dz[j]   = p->B_phi_dz;

        p_ml->B_z[j]        = p->B_z;
        p_ml->B_z_dr[j]     = p->B_z_dr;
        p_ml->B_z_dphi[j]   = p->B_z_dphi;
        p_ml->B_z_dz[j]     = p->B_z_dz;

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
 * @param Bdata pointer to magnetic field data
 */
void particle_ml_to_state(particle_simd_ml* p_ml, int j, particle_state* p,
                          B_field_data* Bdata) {
    a5err err = 0;

    p->rprt       = p_ml->r[j];
    p->phiprt     = p_ml->phi[j];
    p->zprt       = p_ml->z[j];
    p->rdot       = 0;
    p->phidot     = 0;
    p->zdot       = 0;

    p->r          = p_ml->r[j];
    p->phi        = p_ml->phi[j];
    p->z          = p_ml->z[j];
    p->vpar       = p_ml->pitch[j];
    p->mu         = 0;
    p->zeta       = 0;
    p->mass       = 0;
    p->charge     = 0;
    p->time       = p_ml->time[j];
    p->weight     = p_ml->weight[j];
    p->id         = p_ml->id[j];
    p->cputime    = p_ml->cputime[j];
    p->rho        = p_ml->rho[j];
    p->theta      = p_ml->theta[j];
    p->endcond    = p_ml->endcond[j];
    p->walltile   = p_ml->walltile[j];
    p->err        = p_ml->err[j];

    p->B_r        = p_ml->B_r[j];
    p->B_r_dr     = p_ml->B_r_dr[j];
    p->B_r_dphi   = p_ml->B_r_dphi[j];
    p->B_r_dz     = p_ml->B_r_dz[j];

    p->B_phi      = p_ml->B_phi[j];
    p->B_phi_dr   = p_ml->B_phi_dr[j];
    p->B_phi_dphi = p_ml->B_phi_dphi[j];
    p->B_phi_dz   = p_ml->B_phi_dz[j];

    p->B_z        = p_ml->B_z[j];
    p->B_z_dr     = p_ml->B_z_dr[j];
    p->B_z_dphi   = p_ml->B_z_dphi[j];
    p->B_z_dz     = p_ml->B_z_dz[j];

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
 * @param Bdata pointer to magnetic field data
 */
int particle_fo_to_gc(particle_simd_fo* p_fo, int j, particle_simd_gc* p_gc,
                      B_field_data* Bdata) {
    a5err err = p_fo->err[j];
    int simerr = 0; /* Error has already occurred */
    if(err) {simerr = 1;}
    p_gc->id[j]      = p_fo->id[j];
    p_gc->index[j]   = p_fo->index[j];

    real r, phi, z, vpar, mu, zeta, B_dB[15];
    if(!err) {
        real Rprt   = p_fo->r[j];
        real phiprt = p_fo->phi[j];
        real zprt   = p_fo->z[j];
        real vR     = p_fo->rdot[j];
        real vphi   = p_fo->phidot[j] * p_fo->r[j];
        real vz     = p_fo->zdot[j];
        real mass   = p_fo->mass[j];
        real charge = p_fo->charge[j];

        p_gc->mass[j]     = p_fo->mass[j];
        p_gc->charge[j]   = p_fo->charge[j];
        p_gc->weight[j]   = p_fo->weight[j];
        p_gc->time[j]     = p_fo->time[j];
        p_gc->endcond[j]  = p_fo->endcond[j];
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
            Rprt, phiprt, zprt, vR , vphi, vz,
            &r, &phi, &z, &vpar, &mu, &zeta);
    }
    if(!err && r <= 0)  {err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);}
    if(!err && mu < 0)  {err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);}

    real psi[1], rho[1];
    if(!err) {
        err = B_field_eval_B_dB(
            B_dB, r, phi, z, p_fo->time[j], Bdata);
    }
    if(!err) {
        err = B_field_eval_psi(
            psi, r, phi, z, p_fo->time[j], Bdata);
    }
    if(!err) {
        err = B_field_eval_rho(rho, psi[0], Bdata);
    }

    if(!err) {
        p_gc->r[j]          = r;
        p_gc->phi[j]        = phi;
        p_gc->z[j]          = z;
        p_gc->mu[j]         = mu;
        p_gc->zeta[j]       = zeta;
        p_gc->vpar[j]       = vpar;
        p_gc->rho[j]        = rho[0];

        /* Evaluate pol angle so that it is cumulative and at gc position */
        real axis_r = B_field_get_axis_r(Bdata, p_gc->phi[j]);
        real axis_z = B_field_get_axis_z(Bdata, p_gc->phi[j]);
        p_gc->theta[j]  = p_fo->theta[j];
        p_gc->theta[j] += atan2(   (p_fo->r[j]-axis_r) * (p_gc->z[j]-axis_z)
                                 - (p_fo->z[j]-axis_z) * (p_gc->r[j]-axis_r),
                                   (p_fo->r[j]-axis_r) * (p_gc->r[j]-axis_r)
                                 + (p_fo->z[j]-axis_z) * (p_gc->z[j]-axis_z) );

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

void particle_copy_fo(particle_simd_fo* p1, int i, particle_simd_fo* p2, int j) {
        p2->r[j]          = p1->r[i];
        p2->phi[j]        = p1->phi[i];
        p2->z[j]          = p1->z[i];
        p2->rdot[j]       = p1->rdot[i];
        p2->phidot[j]     = p1->phidot[i];
        p2->zdot[j]       = p1->zdot[i];

        p2->time[j]       = p1->time[i];
        p2->cputime[j]    = p1->cputime[i];
        p2->rho[j]        = p1->rho[i];
        p2->weight[j]     = p1->weight[i];
        p2->cputime[j]    = p1->cputime[i];
        p2->rho[j]        = p1->rho[i];
        p2->theta[j]      = p1->theta[i];

        p2->mass[j]       = p1->mass[i];
        p2->charge[j]     = p1->charge[i];

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

void particle_copy_gc(particle_simd_gc* p1, int i, particle_simd_gc* p2, int j) {
    p2->hermite_weights[j]          = p1->hermite_weights[i];
    p2->hermite_knots[j]          = p1->hermite_knots[i];
    p2->use_hermite[j]          = p1->use_hermite[i];
    p2->r[j]          = p1->r[i];
    p2->phi[j]        = p1->phi[i];
    p2->z[j]          = p1->z[i];
    p2->vpar[j]       = p1->vpar[i];
    p2->mu[j]         = p1->mu[i];
    p2->zeta[j]       = p1->zeta[i];

    p2->time[j]       = p1->time[i];
    p2->weight[j]     = p1->weight[i];
    p2->cputime[j]    = p1->cputime[i];
    p2->rho[j]        = p1->rho[i];
    p2->theta[j]      = p1->theta[i];

    p2->mass[j]       = p1->mass[i];
    p2->charge[j]     = p1->charge[i];

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

void particle_copy_ml(particle_simd_ml* p1, int i, particle_simd_ml* p2, int j) {
    p2->r[j]          = p1->r[i];
    p2->phi[j]        = p1->phi[i];
    p2->z[j]          = p1->z[i];
    p2->pitch[j]      = p1->pitch[i];

    p2->time[j]       = p1->time[i];
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
