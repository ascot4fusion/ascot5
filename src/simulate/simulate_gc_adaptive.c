/**
 * @file simulate_gc_adaptive.c
 * @brief Simulate guiding centers using adaptive time-step
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <math.h>
#include "../ascot5.h"
#include "../endcond.h"
#include "../math.h"
#include "../consts.h"
#include "../physlib.h"
#include "../simulate.h"
#include "../particle.h"
#include "../wall.h"
#include "../diag.h"
#include "../B_field.h"
#include "../E_field.h"
#include "../boozer.h"
#include "../mhd.h"
#include "../rfof.h"
#include "../plasma.h"
#include "../boschhale.h"
#include "simulate_gc_adaptive.h"
#include "step/step_gc_cashkarp.h"
#include "mccc/mccc.h"
#include "mccc/mccc_wiener.h"

DECLARE_TARGET_SIMD_UNIFORM(sim)
real simulate_gc_adaptive_inidt(sim_data* sim, particle_simd_gc* p, int i);

#define DUMMY_TIMESTEP_VAL 1.0 /**< Dummy time step value */

/**
 * @brief Simulates guiding centers using adaptive time-step
 *
 * The simulation includes:
 * - orbit-following with Cash-Karp method
 * - Coulomb collisions with Milstein method
 *
 * The simulation is carried until all marker have met some
 * end condition or are aborted/rejected. The final state of the
 * markers is stored in the given marker array. Other output
 * is stored in the diagnostic array.
 *
 * The adaptive time-step is determined by integrator error
 * tolerances as well as user-defined limits for how much
 * marker state can change during a single time-step.
 *
 * @param pq particles to be simulated
 * @param sim simulation data
 *
 */
void simulate_gc_adaptive(particle_queue* pq, sim_data* sim, int mrk_array_size) {

    /* Wiener arrays needed for the adaptive time step */
    mccc_wienarr* wienarr = (mccc_wienarr*) malloc(mrk_array_size*sizeof(mccc_wienarr));

    /* Current time step, suggestions for the next time step and next time
     * step                                                                */
    real* hin = (real*) malloc(mrk_array_size*sizeof(real));
    real* hout_orb = (real*) malloc(mrk_array_size*sizeof(real));
    real* hout_col = (real*) malloc(mrk_array_size*sizeof(real));
    real* hout_rfof = (real*) malloc(mrk_array_size*sizeof(real));
    real* hnext = (real*) malloc(mrk_array_size*sizeof(real));

    Acceleration acceleration;
    acceleration_allocate(&acceleration, mrk_array_size);
    /* Flag indicateing whether a new marker was initialized */
    int* cycle = (int*) malloc(mrk_array_size*sizeof(int));

    real tol_col = sim->ada_tol_clmbcol;
    real tol_orb = sim->ada_tol_orbfol;

    real cputime, cputime_last; // Global cpu time: recent and previous record

    particle_simd_gc p;  // This array holds current states
    particle_simd_gc p0; // This array stores previous states
    particle_allocate_gc(&p, mrk_array_size);
    particle_allocate_gc(&p0, mrk_array_size);
    
    rfof_marker rfof_mrk; // RFOF specific data

    for(int i=0; i< mrk_array_size; i++) {
        p.id[i] = -1;
        p.running[i] = 0;
        acceleration.acc[i] = 1.0;
        acceleration.orbittime[i] = -1;
        acceleration.cross[i].crossed_once = 0;
    }

    /* Initialize running particles */
    int n_running = particle_cycle_gc(pq, &p, &sim->B_data, cycle);

    if(sim->enable_icrh) {
        rfof_set_up(&rfof_mrk, &sim->rfof_data);
    }

    #pragma omp simd
    for(int i = 0; i < mrk_array_size; i++) {
        if(cycle[i] > 0) {
            /* Determine initial time-step */
            hin[i] = simulate_gc_adaptive_inidt(sim, &p, i);
            if(sim->enable_clmbcol) {
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
    real* rnd = (real*) malloc(5*mrk_array_size*sizeof(real));
    particle_offload_gc(&p);
    particle_offload_gc(&p0);
    acceleration_offload(&acceleration,mrk_array_size);
    GPU_MAP_TO_DEVICE(hin[0:mrk_array_size],rnd[0:5*mrk_array_size],hout_orb[0:mrk_array_size],hout_col[0:mrk_array_size],hout_rfof[0:mrk_array_size],hnext[0:mrk_array_size],cycle[0:mrk_array_size])
    mccc_wiener_offload(wienarr,mrk_array_size);   
    while(n_running > 0) {

        /* Store marker states in case time step will be rejected */
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for(int i = 0; i < p.n_mrk; i++) {
            particle_copy_gc(&p, i, &p0, i);
            hout_orb[i]  = DUMMY_TIMESTEP_VAL;
            hout_col[i]  = DUMMY_TIMESTEP_VAL;
            hout_rfof[i] = DUMMY_TIMESTEP_VAL;
            hnext[i]     = DUMMY_TIMESTEP_VAL;
        }

        /*************************** Physics **********************************/

        /* Set time-step negative if tracing backwards in time */
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for(int i = 0; i < p.n_mrk; i++) {
            if(sim->reverse_time) {
                hin[i]  = -hin[i];
            }
        }

        /* Cash-Karp method for orbit-following */
        if(sim->enable_orbfol) {
            if(sim->enable_mhd) {
                step_gc_cashkarp_mhd(
                    &p, hin, hout_orb, tol_orb, &sim->B_data, &sim->E_data,
                    &sim->boozer_data, &sim->mhd_data, sim->enable_aldforce);
            }
            else {
                step_gc_cashkarp(
                    &p, hin, hout_orb, tol_orb, &sim->B_data, &sim->E_data,
                    sim->enable_aldforce);
            }
            /* Check whether time step was rejected */
            GPU_PARALLEL_LOOP_ALL_LEVELS
            for(int i = 0; i < p.n_mrk; i++) {
                /* Switch sign of the time-step again if it was reverted earlier
                */
                if(sim->reverse_time) {
                    hout_orb[i] = -hout_orb[i];
                    hin[i]      = -hin[i];
                }
                if(p.running[i] && hout_orb[i] < 0){
                    p.running[i] = 0;
                    hnext[i] = hout_orb[i];
                }
            }
        }

        /* Milstein method for collisions */
        if(sim->enable_clmbcol) {
            random_normal_simd(sim->random_data, 5*p.n_mrk, rnd);
            mccc_gc_milstein(&p, hin, acceleration.acc, acceleration.collfreq,
			     hout_col, tol_col, wienarr, &sim->B_data,
                             &sim->plasma_data, &sim->mccc_data, rnd);

            /* Check whether time step was rejected */
            GPU_PARALLEL_LOOP_ALL_LEVELS
            for(int i = 0; i < p.n_mrk; i++) {
                if(p.running[i] && hout_col[i] < 0){
                    p.running[i] = 0;
                    hnext[i] =  hout_col[i];
                }
            }
        }

        /* Performs the ICRH kick if in resonance. */
        if(sim->enable_icrh) {
            rfof_resonance_check_and_kick_gc(
                &p, hin, hout_rfof, &rfof_mrk, &sim->rfof_data, &sim->B_data);

            /* Check whether time step was rejected */
            GPU_PARALLEL_LOOP_ALL_LEVELS
            for(int i = 0; i < p.n_mrk; i++) {
                if(p.running[i] && hout_rfof[i] < 0){
                    p.running[i] = 0;
                    hnext[i] =  hout_rfof[i];
                }
            }
        }

        // Reduce weight based on the fusion reactivity
        int n_plasma_species = plasma_get_n_species(&sim->plasma_data);
        //const int * plasma_Z = (int*) malloc(n_plasma_species*sizeof(int));
        //const real * plasma_A = (real*) malloc(n_plasma_species*sizeof(real));

        // Get the atomic numbers and masses of the plasma species
        const int* plasma_Z = plasma_get_species_znum(&sim->plasma_data);
        const real* plasma_A = plasma_get_species_mass(&sim->plasma_data);
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for(int i = 0; i < p.n_mrk; i++) {
            // for now, we assume that p.charge corresponds to Z (the marker is fully ionised)

            int Z = round(p.charge[i]/1.602202e-19);
            int A = round(p.mass[i]/1.660539e-27);

            // Loop through all plasma species and check if we find any of the available reaction

            Reaction reaction[10];
            int species_index_in_reaction[10];
            int n_reactions_found = 0;

            for (int j = 0; j < n_plasma_species; j++) {
                if (plasma_Z[j] == 1 && round(plasma_A[j]/1.660539e-27) == 2) {
                    // D background
                    if (Z == 1 && A == 2) {
                        // Marker is also D ==> add DD reactions to a list of reactions
                        reaction[n_reactions_found] = DD_Tp;
                        species_index_in_reaction[n_reactions_found] = j;
                        n_reactions_found++;
                        reaction[n_reactions_found] = DD_He3n;
                        species_index_in_reaction[n_reactions_found] = j;
                        n_reactions_found++;
                    }
                    if (Z == 2 && A == 3) {
                        // Marker is He3 ==> add DHe3_He4p reaction to a list of reactions
                        reaction[n_reactions_found] = DHe3_He4p;
                        species_index_in_reaction[n_reactions_found] = j;
                        n_reactions_found++;
                    }
                    if (Z == 1 && A == 3) {
                        // Marker is T ==> add DT_He4n reaction to a list of reactions
                        reaction[n_reactions_found] = DT_He4n;
                        species_index_in_reaction[n_reactions_found] = j;
                        n_reactions_found++;
                    }
                }
                if (plasma_Z[j] == 1 && round(plasma_A[j]/1.660539e-27) == 3) {
                    // T background
                    if (Z == 1 && A == 2) {
                        // Marker is D ==> add DT_He4n reaction to a list of reactions
                        reaction[n_reactions_found] = DT_He4n;
                        species_index_in_reaction[n_reactions_found] = j;
                        n_reactions_found++;
                    }
                    if (Z == 2 && A == 3) {
                        // Marker is He3 ==> add DHe3_He4p reaction to a list of reactions
                        reaction[n_reactions_found] = DHe3_He4p;
                        species_index_in_reaction[n_reactions_found] = j;
                        n_reactions_found++;
                    }
                }
            }

            for (int j = 0; j < n_reactions_found; j++) {

                // Get thermal speed for background species
                real ti;
                plasma_eval_temp(&ti, p.rho[i], p.r[i], p.phi[i], p.z[i],
                    p.time[i], species_index_in_reaction[j],
                    &sim->plasma_data);
                real mT = plasma_A[species_index_in_reaction[j]];
                real vt = sqrt(2*ti/mT);

                // Get the fast ion speed
                real bnorm = sqrt(p.B_phi[i] * p.B_phi[i] + p.B_z[i] * p.B_z[i] + p.B_r[i] * p.B_r[i]);
                real pnorm = physlib_gc_p(p.mass[i], p.mu[i], p.ppar[i], bnorm);
                real vf = physlib_vnorm_pnorm(p.mass[i], pnorm);

                // Get <sigma*v> for the correct reaction
                real sigmav = boschhale_sigmav_beam_bulk(
                    reaction[j],
                    vt,
                    vf,
                    1000   // Fixed number of summation intervals for trapz
                );

                // Get n_background of the corresponding reaction
                real bulk_density;
                plasma_eval_dens(&bulk_density, p.rho[i], p.r[i], p.phi[i],
                    p.z[i], p.time[i], species_index_in_reaction[j],
                    &sim->plasma_data);

                /* The number of fusion reactions that "eats" particles,
                i.e., decreases (weight) */
                real N_fusions = p.weight[i] * bulk_density * sigmav * hin[i];
                p.weight[i] -= N_fusions;
            }
        }

        /**********************************************************************/

        cputime = A5_WTIME;
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for(int i = 0; i < p.n_mrk; i++) {
            if(p.id[i] > 0 && !p.err[i]) {
                /* Check other time step limitations */
                if(hnext[i] > 0) {
                    real dphi = fabs(p0.phi[i]-p.phi[i]) / sim->ada_max_dphi;
                    real drho = fabs(p0.rho[i]-p.rho[i]) / sim->ada_max_drho;

                    if(dphi > 1 && dphi > drho) {
                        hnext[i] = -hin[i]/dphi;
                    }
                    else if(drho > 1 && drho > dphi) {
                        hnext[i] = -hin[i]/drho;
                    }
                }

                /* Retrieve marker states in case time step was rejected      */
                if(hnext[i] < 0) {
                    particle_copy_gc(&p0, i, &p, i);
                }
                if(p.running[i]){

                    /* Advance time (if time step was accepted) and determine
                       next time step */
                    if(hnext[i] < 0){
                        /* if hnext < 0, you screwed up and had to copy the
                        previous state. Therefore, let us use the suggestion
                        given by the integrator when retaking the failed step.*/
                        hin[i] = -hnext[i];
                    }
                    else {
                        p.time[i] += ( 1.0 - 2.0 * ( sim->reverse_time > 0 ) )
                            * hin[i] * acceleration.acc[i];
                        p.mileage[i] += hin[i] * acceleration.acc[i];
                        if(acceleration.orbittime[i] >= 0)
                            acceleration.orbittime[i] += hin[i] * acceleration.acc[i];
                        /* In case the time step was succesful, pick the
                        smallest recommended value for the next step */
                        if(hnext[i] > hout_orb[i]) {
                            /* Use time step suggested by the orbit-following
                               integrator */
                            hnext[i] = hout_orb[i];
                        }
                        if(hnext[i] > hout_col[i]) {
                            /* Use time step suggested by the collision
                               integrator */
                            hnext[i] = hout_col[i];
                        }
                        if(hnext[i] > hout_rfof[i]) {
                            /* Use time step suggested by RFOF */
                            hnext[i] = hout_rfof[i];
                        }
                        if(hnext[i] == 1.0) {
                            /* Time step is unchanged (happens when no physics
                               are enabled) */
                            hnext[i] = hin[i];
                        }
                        hin[i] = hnext[i];
                        if(sim->enable_clmbcol) {
                            /* Clear wiener processes */
                            mccc_wiener_clean(&(wienarr[i]), p.time[i]);
                        }
                    }

                    p.cputime[i] += cputime - cputime_last;
                }
            }
        }
        cputime_last = cputime;

        /* If OMP is crossed, adjust acceleration */
        if(sim->enable_ada > 1) {
            recalculate_acceleration(&acceleration, sim, &p, &p0);
        }

        /* Check possible end conditions */
        endcond_check_gc(&p, &p0, sim);

        /* Update diagnostics */
        diag_update_gc(&sim->diag_data, &sim->B_data, &p, &p0);

        /* Update number of running particles */
#ifdef GPU
        n_running = 0;
        GPU_PARALLEL_LOOP_ALL_LEVELS_REDUCTION(n_running)
        for(int i = 0; i < p.n_mrk; i++)
        {
            if(p.running[i] > 0) n_running++;
        }
#else
        n_running = particle_cycle_gc(pq, &p, &sim->B_data, cycle);

        /* Determine simulation time-step for new particles */
        #pragma omp simd
        for(int i = 0; i <p.n_mrk; i++) {
            if(cycle[i] > 0) {
                hin[i] = simulate_gc_adaptive_inidt(sim, &p, i);
                acceleration.acc[i] = 1.0;
                acceleration.orbittime[i] = -1;
                acceleration.cross[i].crossed_once = 0;
                if(sim->enable_clmbcol) {
                    /* Re-allocate array storing the Wiener processes */
                    mccc_wiener_initialize(&(wienarr[i]), p.time[i]);
                }
                if(sim->enable_icrh) {
                    /* Reset icrh (rfof) resonance memory matrix. */
                    rfof_clear_history(&rfof_mrk, i);
                }
            }
        }
#endif
    }

    /* All markers simulated! */
#ifdef GPU
    GPU_MAP_FROM_DEVICE(sim[0:1])
    particle_onload_gc(&p);
    n_running = particle_cycle_gc(pq, &p, &sim->B_data, cycle);
    mccc_wiener_onload(wienarr,mrk_array_size);   
#endif
    free(cycle);
    free(hin);
    free(rnd);
    free(hout_orb);
    free(hout_col);
    free(hout_rfof);
    free(hnext);
    /* Deallocate rfof structs */
    if(sim->enable_icrh) {
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
real simulate_gc_adaptive_inidt(sim_data* sim, particle_simd_gc* p, int i) {
    /* Just use some large value if no physics are defined */
    real h = DUMMY_TIMESTEP_VAL;

    /* Value defined directly by user */
    if(sim->fix_usrdef_use) {
        h =  sim->fix_usrdef_val;
    }
    else {
        /* Value calculated from gyrotime */
        if(sim->enable_orbfol) {
            real Bnorm = math_normc(p->B_r[i], p->B_phi[i], p->B_z[i]);
            real gyrotime = CONST_2PI /
                phys_gyrofreq_ppar(p->mass[i], p->charge[i], p->mu[i],
                                   p->ppar[i], Bnorm);
            if(h > gyrotime) {
                h = gyrotime;
            }
        }

        /* Value calculated from collision frequency */
        if(sim->enable_clmbcol) {
            real nu = 1;
            /*mccc_collfreq_gc(p, &sim->B_data, &sim->plasma_data,
                sim->coldata, &nu, i); */

            /* Only small angle collisions so divide this by 100 */
            real colltime = 1/(100*nu);
            if(h > colltime) {h=colltime;}
        }
    }
    return h;
}


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
 */
void recalculate_acceleration(
    Acceleration* acc, sim_data* sim, particle_simd_gc* p, particle_simd_gc* p0)
{
    real rz[2];
    real SAFETY_FACTOR = (float)sim->enable_ada / 1000.0;
    GPU_PARALLEL_LOOP_ALL_LEVELS
    for(int i = 0; i < p->n_mrk; i++) {
        B_field_get_axis_rz(rz, &sim->B_data, p->phi[i]);
        int omp_crossed = ((p->z[i] - rz[1]) * (p0->z[i] - rz[1]) < 0) &&
            p->r[i] > rz[0];
        if(omp_crossed && acc->cross[i].crossed_twice) {
            if( ((float)acc->cross[i].first_ppar - 0.5) * p->ppar[i] > 0 ) {
                acc->acc[i] = fmax(1.0, SAFETY_FACTOR / (acc->orbittime[i] * acc->collfreq[i]));
            }
            else {
                acc->acc[i] = 1;
            }
            acc->cross[i].crossed_once = 1;
            acc->cross[i].crossed_twice = 0;
            acc->cross[i].first_ppar = p->ppar[i] > 0;
            acc->orbittime[i] = 0;
        }
        else if(omp_crossed && acc->cross[i].crossed_once) {
            acc->cross[i].crossed_twice = 1;
            if( ((float)acc->cross[i].first_ppar - 0.5) * p->ppar[i] > 0 ) {
                acc->acc[i] = fmax(1.0, SAFETY_FACTOR / (acc->orbittime[i] * acc->collfreq[i]));
                acc->cross[i].crossed_once = 1;
                acc->cross[i].crossed_twice = 0;
                acc->cross[i].first_ppar = p->ppar[i] > 0;
                acc->orbittime[i] = 0;
            }
        }
        else if(omp_crossed) {
            acc->cross[i].crossed_once = 1;
            acc->cross[i].first_ppar = p->ppar[i] > 0;
            acc->orbittime[i] = 0;
        }
    }
}

/**
 * @brief Allocates struct representing acceleration struc
 *
 * Size used for memory allocation is NSIMD for CPU run and the total number
 * of particles for GPU.
 *
 * @param Acceleration struct to allocate
 * @param nmrk the number of markers that the struct represents
 */
void acceleration_allocate(Acceleration* acceleration, int nmrk){
  acceleration->acc       = malloc(nmrk * sizeof(acceleration->acc) );
  acceleration->orbittime = malloc(nmrk * sizeof(acceleration->orbittime) );
  acceleration->collfreq  = malloc(nmrk * sizeof(acceleration->collfreq) );
  acceleration->cross     = malloc(nmrk * sizeof(acceleration->cross) );
}
/**
 * @brief Offload acceleration struct to GPU.
 *
 * @param acceleration pointer to the acceleration struct to be offloaded.
 */
void acceleration_offload(Acceleration* acceleration, int mrk_array_size) {
    GPU_MAP_TO_DEVICE(
        acceleration[0:1],\
        acceleration->acc       [0:mrk_array_size],\
        acceleration->orbittime [0:mrk_array_size],\
        acceleration->collfreq  [0:mrk_array_size],\
        acceleration->cross     [0:mrk_array_size]
    )
}
