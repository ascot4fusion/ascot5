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
#include "particle.h"
#include "plasma.h"
#include "wall.h"
#include "boozer.h"
#include "mhd.h"
#include "neutral.h"
#include "B_field.h"
#include "E_field.h"
#include "random.h"
#include "simulate.h"
#include "print.h"
#include "simulate/simulate_ml_adaptive.h"
#include "simulate/simulate_gc_adaptive.h"
#include "simulate/simulate_gc_fixed.h"
#include "simulate/simulate_fo_fixed.h"
#include "simulate/mccc/mccc.h"
#include "gctransform.h"
#include "asigma.h"
#include "rfof_interface.h"

void sim_monitor(char* filename, volatile int* n, volatile int* finished);

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
 * @param n_particles total number of markers to be simulated
 * @param p pointer to array storing all marker states to be simulated
 * @param sim_offload pointer to simulation offload data
 * @param offload_data pointer to the rest of the offload data
 * @param offload_array pointer to input data offload array
 * @param int_offload_array pointer to input data int offload array
 * @param diag_offload_array pointer to diagnostics offload array
 *
 * @todo Reorganize this function so that it conforms to documentation.
 */
void simulate(int n_particles, particle_state* p, sim_data* sim) {

    // Size = NSIMD on CPU and Size = Total number of particles on GPU
    int n_queue_size;
#ifdef GPU
    n_queue_size = n_particles;
#else
    n_queue_size = NSIMD;
#endif
    /**************************************************************************/
    /* 1. Input offload data is unpacked and initialized by calling           */
    /*    respective init functions.                                          */
    /*                                                                        */
    /**************************************************************************/

    simulate_init(sim);

#ifdef GPU
    if(sim->sim_mode != 1) {
        print_err("Only GO mode ported to GPU. Please set SIM_MODE=1.");
        exit(1);
    }
    if(sim->record_mode) {
        print_err("RECORD_MODE=1 not ported to GPU. Please disable it.");
        exit(1);
    }
    if(sim->enable_atomic) {
        print_err("Atomic not yet ported to GPU. Please set ENABLE_ATOMIC=0.");
        exit(1);
    }
    if(sim->enable_mhd) {
        print_err("MHD not yet ported to GPU. Please set ENABLE_MHD=0.");
        exit(1);
    }
    if(sim->diag_data.diagorb_collect) {
        print_err(
            "ENABLE_ORBITWRITE=1 not ported to GPU. Please disable it.");
        exit(1);
    }
    if(sim->diag_data.diagtrcof_collect) {
        print_err(
            "ENABLE_TRANSCOEF=1 not ported to GPU. Please disable it.");
        exit(1);
    }
#endif

    if(sim->enable_icrh) {
        char *xml_filename = "rfof_codeparam.xml";
       rfof_interface_initev_excl_marker_stuff(xml_filename, &(sim->rfof_data));
    }

    diag_init(&sim->diag_data, n_particles);
    GPU_MAP_TO_DEVICE(sim[0:1])
    B_field_offload(&sim->B_data);
    E_field_offload(&sim->E_data);
    plasma_offload(&sim->plasma_data);
    neutral_offload(&sim->neutral_data);
    wall_offload(&sim->wall_data);
    boozer_offload(&sim->boozer_data);
    mhd_offload(&sim->mhd_data);
    asigma_offload(&sim->asigma_data);
    diag_offload(&sim->diag_data);

    /**************************************************************************/
    /* 2. Meta data (e.g. random number generator) is initialized.            */
    /*                                                                        */
    /**************************************************************************/
    random_init(&sim->random_data, 0);

    /**************************************************************************/
    /* 3. Markers are put into simulation queue.                              */
    /*                                                                        */
    /**************************************************************************/
    particle_queue pq;

    pq.n = 0;
    for(int i = 0; i < n_particles; i++) {
        pq.n++;
    }

    pq.p = (particle_state**) malloc(pq.n * sizeof(particle_state*));
    pq.finished = 0;

    pq.next = 0;
    for(int i = 0; i < n_particles; i++) {
        pq.p[pq.next++] = &p[i];

    }
    pq.next = 0;

    print_out(VERBOSE_NORMAL, "Simulation begins; %d threads.\n",
              omp_get_max_threads());
    fflush(stdout);

    /**************************************************************************/
    /* 4. Threads are spawned. One thread is dedicated for monitoring         */
    /*    progress, if monitoring is active.                                  */
    /*                                                                        */
    /**************************************************************************/
#ifndef GPU
    omp_set_max_active_levels(2);
#endif
#if !defined(GPU) && VERBOSE > 1
    #pragma omp parallel sections num_threads(2)
    {
        #pragma omp section
#endif
        {
            /******************************************************************/
            /* 5. Other threads execute marker simulation using the mode the  */
            /*    user has chosen.                                            */
            /*                                                                */
            /******************************************************************/
            if(pq.n > 0 && (sim->sim_mode == simulate_mode_gc
                        || sim->sim_mode == simulate_mode_hybrid)) {
                if(sim->enable_ada) {
                    OMP_PARALLEL_CPU_ONLY
                    simulate_gc_adaptive(&pq, sim);
                }
                else {
                    OMP_PARALLEL_CPU_ONLY
                    simulate_gc_fixed(&pq, sim);
                }
            }
            else if(pq.n > 0 && sim->sim_mode == simulate_mode_fo) {
                OMP_PARALLEL_CPU_ONLY
                simulate_fo_fixed(&pq, sim, n_queue_size);
            }
            else if(pq.n > 0 && sim->sim_mode == simulate_mode_ml) {
                OMP_PARALLEL_CPU_ONLY
                simulate_ml_adaptive(&pq, sim);
            }
        }
#if !defined(GPU) && VERBOSE > 1
        #pragma omp section
        {
            /* Update progress until simulation is complete.             */
            /* Trim .h5 from filename and replace it with _<QID>.stdout  */
            if(id == 0) {
                char filename[519], outfn[256];
                strcpy(outfn, sim->hdf5_out);
                outfn[strlen(outfn)-3] = '\0';
                sprintf(filename, "%s_%s.stdout", outfn, sim->qid);
                sim_monitor(filename, &pq.n, &pq.finished);
            }
        }
    }
#endif

    /**************************************************************************/
    /* 6. (If hybrid mode is active) Markers with hybrid end condition active */
    /*    are placed on a new queue, and they have their end condition        */
    /*    deactivated and they are simulated with simulate_fo_fixed.c until   */
    /*    they have met some other end condition. Threads are spawned and     */
    /*    progress is monitored as previously.                                */
    /*                                                                        */
    /**************************************************************************/
    int n_new = 0;
    if(sim->sim_mode == simulate_mode_hybrid) {

        /* Determine the number markers that should be run
         * in fo after previous gc simulation */
        for(int i = 0; i < pq.n; i++) {
            if(pq.p[i]->endcond == endcond_hybrid) {
                /* Check that there was no wall between when moving from
                   gc to fo */
                real w_coll;
                int tile = wall_hit_wall(pq.p[i]->r, pq.p[i]->phi, pq.p[i]->z,
                        pq.p[i]->rprt, pq.p[i]->phiprt, pq.p[i]->zprt,
                                         &sim->wall_data, &w_coll);
                if(tile > 0) {
                    pq.p[i]->walltile = tile;
                    pq.p[i]->endcond |= endcond_wall;
                }
                else {
                    n_new++;
                }
            }
        }
    }
    if(n_new > 0) {

        /* Reset hybrid marker end condition */
        for(int i = 0; i < pq.n; i++) {
            if(pq.p[i]->endcond & endcond_hybrid) {
                pq.p[i]->endcond ^= endcond_hybrid;
            }
        }
        pq.next = 0;
        pq.finished = 0;

#if !defined(GPU) && VERBOSE > 1
        #pragma omp parallel sections num_threads(2)
        {
            #pragma omp section
#endif
            {
                OMP_PARALLEL_CPU_ONLY
                simulate_fo_fixed(&pq, sim, n_queue_size);
            }
#if !defined(GPU) && VERBOSE > 1
            #pragma omp section
            {
                /* Trim .h5 from filename and replace it with _<qid>.stdout */
                if(id == 0) {
                    char filename[519], outfn[256];
                    strcpy(outfn, sim->hdf5_out);
                    outfn[strlen(outfn)-3] = '\0';
                    sprintf(filename, "%s_%s.stdout", outfn, sim->qid);
                    sim_monitor(filename, &pq.n, &pq.finished);
                }
            }
        }
#endif
    }

    /**************************************************************************/
    /* 7. Simulation data is deallocated.                                     */
    /**************************************************************************/
    free(pq.p);

    if(sim->enable_icrh) {
        rfof_interface_deallocate_rfof_input_param(
            &(sim->rfof_data.cptr_rfof_input_params));
        rfof_interface_deallocate_rfglobal(&(sim->rfof_data.cptr_rfglobal));
    }
    /**************************************************************************/
    /* 8. Execution returns to host where this function was called.           */
    /*                                                                        */
    /**************************************************************************/
    print_out(VERBOSE_NORMAL, "Simulation complete.\n");
}

/**
 * @brief Initialize simulation data struct
 *
 * @param sim pointer to data struct
 */
void simulate_init(sim_data* sim) {

    mccc_init(&sim->mccc_data, !sim->disable_energyccoll,
              !sim->disable_pitchccoll, !sim->disable_gcdiffccoll);

    if(sim->disable_gctransform) {
        gctransform_setorder(0);
    }
    asigma_extrapolate(sim->enable_atomic==2);

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
 * @param filename pointer to file where progress is written. File is opened and
 *        closed outside this function
 * @param n pointer to number of total markers in simulation queue
 * @param finished pointer to number of finished markers in simulation queue
 */
void sim_monitor(char* filename, volatile int* n, volatile int* finished) {
    /* Open a file for writing simulation progress */
    FILE *f = fopen(filename, "w");
    if (f == NULL) {
        print_out(VERBOSE_DEBUG,
                  "Warning. %s could not be opened for progress updates.\n",
                  filename);
        return;
    }

    real time_sim_started = A5_WTIME;
    int stopflag = 1; /* Ensures progress is written one last time at 100% */
    int n_temp, finished_temp; /* Use these to store volatile variables so that
                                  their value does not change during one loop */
    while(stopflag) {
        n_temp = *n;
        finished_temp = *finished;
        real fracprog = ((real) finished_temp)/n_temp;
        real timespent = (A5_WTIME)-time_sim_started;

        if(n_temp == finished_temp) {
            stopflag = 0;
        }

        if(fracprog == 0) {
            fprintf(f, "No marker has finished simulation yet. "
                    "Time spent: %.2f h\n", timespent/3600);
        }
        else {
            fprintf(f, "Progress: %d/%d, %.2f %%. Time spent: %.2f h, "
                    "estimated time to finish: %.2f h\n", finished_temp, n_temp,
                    100*fracprog, timespent/3600, (1/fracprog-1)*timespent/3600);
        }
        fflush(f);
        //sleep(A5_PRINTPROGRESSINTERVAL);
    }

    fprintf(f, "Simulation finished.\n");
    fclose(f);
}
