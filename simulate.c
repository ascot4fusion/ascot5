/**
 * @file simulate.c
 * @brief Simulation is initialized and run from here
 *
 * This module acts as an interface through which different types of simulations
 * are initialized and run. This module handles no IO operations (with the
 * exception of writing of progress update), no offloading (only unpacking and
 * initialization is done here).
 *
 * Thread level parallelisation is done here and the threads have shared access
 * on the data once it has been initialized. However, threads should only modify
 * marker and diagnostic data.
 */
#include <string.h>
#include <unistd.h>
#include "endcond.h"
#include "hdf5io/hdf5_orbits.h"
#include "offload.h"
#include "particle.h"
#include "plasma.h"
#include "random.h"
#include "simulate.h"
#include "simulate/simulate_ml_adaptive.h"
#include "simulate/simulate_gc_adaptive.h"
#include "simulate/simulate_gc_fixed.h"
#include "simulate/simulate_fo_fixed.h"
#include "simulate/mccc/mccc_coefs.h"

#pragma omp declare target
void sim_init(sim_data* sim, sim_offload_data* offload_data);
void sim_monitor(FILE* f, volatile int* n, volatile int* finished);
#pragma omp end declare target

/**
 * @brief Execute marker simulation
 *
 * This simulates markers using given inputs and options. All different types of
 * simulations are initialized and run via this function.
 *
 * This function proceeds as follows:
 *
 * 1. Input offload data is unpacked and initialized by calling respective init
 *    functions.
 *
 * 2. Meta data (e.g. random number generator) is initialized.
 *
 * 3. Markers are put into simulation queue.
 *
 * 4. Threads are spawned. One thread is dedicated for monitoring progress, if
 *    monitoring is active.
 *
 * 5. Other threads execute marker simulation using the mode the user has
 *    chosen.
 *
 * -  Process continues once all markers have been simulated and each thread has
 *    finished. Monitoring is also terminated.
 *
 * 6. (If hybrid mode is active) Markers with hybrid end condition active are
 *    placed on a new queue, and they have their end condition deactivated and
 *    they are simulated with simulate_fo_fixed.c until they have met some other
 *    end condition. Threads are spawned and progress is monitored as
 *    previously.
 *
 * 7. Simulation data is deallocated except for data that is mapped back to
 *    host.
 *
 * 8. Execution returns to host where this function was called.
 *
 * @param id target id where this function is executed, zero if on host
 * @param n_particles total number of markers to be simulated
 * @param p pointer to array storing all marker states to be simulated
 * @param sim_offload pointer to simulation offload data
 * @param offload_data pointer to the rest of the offload data
 * @param offload_array pointer to input data offload array
 * @param diag_offload_array pointer to diagnostics offload array
 *
 * @todo Reorganize this function so that it conforms to documentation.
 */
void simulate(int id, int n_particles, particle_state* p,
        sim_offload_data* sim_offload,
        offload_package* offload_data,
        real* offload_array,
        real* diag_offload_array) {

    char targetname[5];
    if(id == 0) {
        sprintf(targetname, "host");
    }
    else {
        sprintf(targetname, "mic%d", id-1);
    }
    /**************************************************************************/
    /* 1. Input offload data is unpacked and initialized by calling           */
    /*    respective init functions.                                          */
    /*                                                                        */
    /**************************************************************************/
    sim_data sim;
    sim_init(&sim, sim_offload);

    real* ptr;
    ptr = offload_unpack(offload_data, offload_array,
            sim_offload->B_offload_data.offload_array_length);
    B_field_init(&sim.B_data, &sim_offload->B_offload_data, ptr);

    ptr = offload_unpack(offload_data, offload_array,
            sim_offload->E_offload_data.offload_array_length);
    E_field_init(&sim.E_data, &sim_offload->E_offload_data, ptr);

    ptr = offload_unpack(offload_data, offload_array,
            sim_offload->plasma_offload_data.offload_array_length);
    plasma_init(&sim.plasma_data, &sim_offload->plasma_offload_data, ptr);

    ptr = offload_unpack(offload_data, offload_array,
            sim_offload->neutral_offload_data.offload_array_length);
    neutral_init(&sim.neutral_data, &sim_offload->neutral_offload_data, ptr);

    ptr = offload_unpack(offload_data, offload_array,
            sim_offload->wall_offload_data.offload_array_length);
    wall_init(&sim.wall_data, &sim_offload->wall_offload_data, ptr);

    diag_init(&sim.diag_data, &sim_offload->diag_offload_data,
            diag_offload_array);

    /**************************************************************************/
    /* 2. Meta data (e.g. random number generator) is initialized.            */
    /*                                                                        */
    /**************************************************************************/

    /* Initialize collision coefficients */
    sim.coldata = NULL;
    mccc_coefs_init(sim.coldata);

    /**************************************************************************/
    /* 3. Markers are put into simulation queue.                              */
    /*                                                                        */
    /**************************************************************************/
    particle_queue pq;
    particle_queue pq_hybrid;

    pq.n = 0;
    pq_hybrid.n = 0;
    for(int i = 0; i < n_particles; i++) {
        if(p[i].endcond == 0) {
            pq.n++;
        }
        else if(p[i].endcond == endcond_hybrid) {
            pq_hybrid.n++;
        }
    }

    pq.p = (particle_state**) malloc(pq.n * sizeof(particle_state*));
    pq_hybrid.p = (particle_state**)malloc(pq_hybrid.n*sizeof(particle_state*));

    pq.finished = 0;
    pq_hybrid.finished = 0;

    pq.next = 0;
    pq_hybrid.next = 0;
    for(int i = 0; i < n_particles; i++) {
        if(p[i].endcond == 0) {
            pq.p[pq.next++] = &p[i];
        }
        else if(p[i].endcond == endcond_hybrid) {
            pq_hybrid.p[pq_hybrid.next++] = &p[i];
        }

    }
    pq.next = 0;

    random_init(&sim.random_data, 0);

    #if VERBOSE >= 1
    printf("%s: All fields initialized. Simulation begins, %d threads.\n",
        targetname, omp_get_max_threads());
    #endif

#if VERBOSE > 1
    /* Open a file for writing simulation progress */
    char filename[256];
    sprintf(filename, "%s_%06d.stdout", sim_offload->outfn,
        sim_offload->mpi_rank);
    FILE *f = fopen(filename, "w");
    if (f == NULL) {
        printf("%s: Error opening stdout file.\n", targetname);
    }
#endif

    /**************************************************************************/
    /* 4. Threads are spawned. One thread is dedicated for monitoring         */
    /*    progress, if monitoring is active.                                  */
    /*                                                                        */
    /**************************************************************************/
    #pragma omp parallel sections num_threads(2)
    {
        #pragma omp section
        {
            /******************************************************************/
            /* 5. Other threads execute marker simulation using the mode the  */
            /*    user has chosen.                                            */
            /*                                                                */
            /******************************************************************/
            if(pq.n > 0 && (sim.sim_mode == simulate_mode_gc
                        || sim.sim_mode == simulate_mode_hybrid)) {
                sim.diag_data.orbits.type = diag_orb_type_gc;
                if(sim.enable_ada) {
                    #pragma omp parallel
                    simulate_gc_adaptive(&pq, &sim);
                }
                else {
                    #pragma omp parallel
                    simulate_gc_fixed(&pq, &sim);
                }
            }
            else if(pq.n > 0 && sim.sim_mode == simulate_mode_fo) {
                if(sim.record_GOasGC) {
                    sim.diag_data.orbits.type = diag_orb_type_gc;
                }
                else {
                    sim.diag_data.orbits.type = diag_orb_type_fo;
                }

                #pragma omp parallel
                simulate_fo_fixed(&pq, &sim);
            }
            else if(pq.n > 0 && sim.sim_mode == simulate_mode_ml) {
                sim.diag_data.orbits.type = diag_orb_type_ml;

                #pragma omp parallel
                simulate_ml_adaptive(&pq, &sim);
            }
        }

        #pragma omp section
        {
#if VERBOSE > 1
            /* Update progress until simulation is complete. */
            sim_monitor(f, &pq.n, &pq.finished);
#endif
        }
    }

    /**************************************************************************/
    /* 6. (If hybrid mode is active) Markers with hybrid end condition active */
    /*    are placed on a new queue, and they have their end condition        */
    /*    deactivated and they are simulated with simulate_fo_fixed.c until   */
    /*    they have met some other end condition. Threads are spawned and     */
    /*    progress is monitored as previously.                                */
    /*                                                                        */
    /**************************************************************************/
    if(sim.sim_mode == simulate_mode_hybrid) {

        /* Determine the number markers that should be run
         * in fo after previous gc simulation */
        int n_new = 0;
        for(int i = 0; i < pq.n; i++) {
            if(pq.p[i]->endcond == endcond_hybrid) {
                /* Check that there was no wall between when moving from
                   gc to fo */
                int tile = wall_hit_wall(pq.p[i]->r, pq.p[i]->phi, pq.p[i]->z,
                        pq.p[i]->rprt, pq.p[i]->phiprt, pq.p[i]->zprt,
                        &sim.wall_data);
                if(tile > 0) {
                    pq.p[i]->walltile = tile;
                    pq.p[i]->endcond |= endcond_wall;
                }
                else {
                    n_new++;
                }
            }
        }

        if(n_new > 0) {
        /* Reallocate and add "old" hybrid particles to the hybrid queue */
            particle_state** tmp = pq_hybrid.p;
            pq_hybrid.p = (particle_state**) malloc((pq_hybrid.n + n_new)
                    * sizeof(particle_state*));
            memcpy(pq_hybrid.p, tmp, pq_hybrid.n * sizeof(particle_state*));
            free(tmp);

            /* Add "new" hybrid particles and reset their end condition */
            pq_hybrid.n += n_new;
            for(int i = 0; i < pq.n; i++) {
                if(pq.p[i]->endcond == endcond_hybrid) {
                    pq.p[i]->endcond = 0;
                    pq_hybrid.p[pq_hybrid.next++] = pq.p[i];
                }
            }
        }
        pq_hybrid.next = 0;

        sim.record_GOasGC = 1;//Make sure we don't collect fos in gc diagnostics
        sim.diag_data.orbits.type = diag_orb_type_gc;

        #pragma omp parallel sections num_threads(2)
        {
            #pragma omp section
            {
                #pragma omp parallel
                simulate_fo_fixed(&pq_hybrid, &sim);
            }

            #pragma omp section
            {
#if VERBOSE > 1
                sim_monitor(f, &pq_hybrid.n, &pq_hybrid.finished);
#endif
            }
        }
    }

    /**************************************************************************/
    /* 7. Simulation data is deallocated except for data that is mapped back  */
    /*    to host.                                                            */
    /*                                                                        */
    /**************************************************************************/
#if VERBOSE > 1
    /* Close progress file*/
    fclose(f);
#endif

    free(sim.coldata);
    free(pq.p);
    free(pq_hybrid.p);

    // Temporary solution
#ifndef TARGET
    hdf5_orbits_write(&sim, sim_offload->hdf5_out, sim_offload->qid);
#endif

    diag_clean(&sim.diag_data);

    /**************************************************************************/
    /* 8. Execution returns to host where this function was called.           */
    /*                                                                        */
    /**************************************************************************/
    #if VERBOSE >= 1
    printf("%s: Simulation complete.\n", targetname);
    #endif
}

/**
 * @brief Initialize simulation data struct on target
 *
 * This function copies the simulation parameters from the offload struct
 * to the struct on the target.
 *
 * @param sim pointer to data struct on target
 * @param offload_data pointer to offload data struct
 */
void sim_init(sim_data* sim, sim_offload_data* offload_data) {
    sim->sim_mode             = offload_data->sim_mode;
    sim->enable_ada           = offload_data->enable_ada;
    sim->record_GOasGC        = offload_data->record_GOasGC;

    sim->fix_usrdef_use       = offload_data->fix_usrdef_use;
    sim->fix_usrdef_val       = offload_data->fix_usrdef_val;
    sim->fix_stepsPerGO       = offload_data->fix_stepsPerGO;

    sim->ada_tol_orbfol       = offload_data->ada_tol_orbfol;
    sim->ada_tol_clmbcol      = offload_data->ada_tol_clmbcol;
    sim->ada_max_drho         = offload_data->ada_max_drho;
    sim->ada_max_dphi         = offload_data->ada_max_dphi;
    sim->ada_max_acc          = offload_data->ada_max_acc;

    sim->enable_orbfol        = offload_data->enable_orbfol;
    sim->enable_clmbcol       = offload_data->enable_clmbcol;

    sim->endcond_active       = offload_data->endcond_active;
    sim->endcond_maxSimTime   = offload_data->endcond_maxSimTime;
    sim->endcond_maxCpuTime   = offload_data->endcond_maxCpuTime;
    sim->endcond_minRho       = offload_data->endcond_minRho;
    sim->endcond_maxRho       = offload_data->endcond_maxRho;
    sim->endcond_minEkin      = offload_data->endcond_minEkin;
    sim->endcond_minEkinPerTe = offload_data->endcond_minEkinPerTe;
    sim->endcond_maxTorOrb    = offload_data->endcond_maxTorOrb;
    sim->endcond_maxPolOrb    = offload_data->endcond_maxPolOrb;
}

/**
 * @brief Monitor simulation progress
 *
 * This function contains a loop that is repeated until all markers have
 * finished simulation. Loops are executed at interval defined by
 * A5_PRINTPROGRESSINTERVAL in ascot5.h.
 *
 * At each loop, number of markers that have finished simulation is written
 * to output file, along with time spent on simulation and estimated time
 * remaining for the simulation to finish.
 *
 * @param f pointer to file where progress is written. File is opened and closed
 *          outside this function
 * @param n pointer to number of total markers in simulation queue
 * @param finished pointer to number of finished markers in simulation queue
 */
void sim_monitor(FILE* f, volatile int* n, volatile int* finished) {
    real time_sim_started = A5_WTIME;
    while(f != NULL && *n > *finished) {
        real fracprog = ((real) *finished)/(*n);
        real timespent = (A5_WTIME)-time_sim_started;

        if(fracprog == 0) {
            fprintf(f, "No marker has finished simulation yet. "
                    "Time spent: %.2f h\n", timespent/3600);
        }
        else {
            fprintf(f, "Progress: %d/%d, %.2f %%. Time spent: %.2f h, "
                    "estimated time to finish: %.2f h\n", *finished, *n,
                    100*fracprog, timespent/3600, (1/fracprog-1)*timespent/3600);
        }
        fflush(f);
        sleep(A5_PRINTPROGRESSINTERVAL);
    }
}
