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
#include "offload.h"
#include "particle.h"
#include "plasma.h"
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

#pragma omp declare target
void sim_monitor(char* filename, volatile int* n, volatile int* finished);
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
 * @param int_offload_array pointer to input data int offload array
 * @param diag_offload_array pointer to diagnostics offload array
 *
 * @todo Reorganize this function so that it conforms to documentation.
 */
void simulate(
    int id, int n_particles, particle_state* p, sim_offload_data* sim_offload,
    offload_package* offload_data, real* offload_array, int* int_offload_array,
    real* diag_offload_array) {

    char targetname[5];
    if(id == 0) {
        sprintf(targetname, "host");
    }
    else {
        sprintf(targetname, "mic%hu", (unsigned short)(id-1));
    }
    /**************************************************************************/
    /* 1. Input offload data is unpacked and initialized by calling           */
    /*    respective init functions.                                          */
    /*                                                                        */
    /**************************************************************************/
    sim_data sim;
    sim_init(&sim, sim_offload);

    real* ptr; int* ptrint;
    offload_unpack(offload_data, offload_array,
                   sim_offload->B_offload_data.offload_array_length,
                   NULL, 0, &ptr, &ptrint);
    B_field_init(&sim.B_data, &sim_offload->B_offload_data, ptr);

    offload_unpack(offload_data, offload_array,
                   sim_offload->E_offload_data.offload_array_length,
                   NULL, 0, &ptr, &ptrint);
    E_field_init(&sim.E_data, &sim_offload->E_offload_data, ptr);

    offload_unpack(offload_data, offload_array,
                   sim_offload->plasma_offload_data.offload_array_length,
                   NULL, 0, &ptr, &ptrint);
    plasma_init(&sim.plasma_data, &sim_offload->plasma_offload_data, ptr);

    offload_unpack(offload_data, offload_array,
                   sim_offload->neutral_offload_data.offload_array_length,
                   NULL, 0, &ptr, &ptrint);
    neutral_init(&sim.neutral_data, &sim_offload->neutral_offload_data, ptr);

    offload_unpack(offload_data, offload_array,
                   sim_offload->wall_offload_data.offload_array_length,
                   int_offload_array,
                   sim_offload->wall_offload_data.int_offload_array_length,
                   &ptr, &ptrint);
    wall_init(&sim.wall_data, &sim_offload->wall_offload_data, ptr, ptrint);

    offload_unpack(offload_data, offload_array,
                   sim_offload->boozer_offload_data.offload_array_length,
                   NULL, 0, &ptr, &ptrint);
    boozer_init(&sim.boozer_data, &sim_offload->boozer_offload_data, ptr);

    offload_unpack(offload_data, offload_array,
                   sim_offload->mhd_offload_data.offload_array_length,
                   NULL, 0, &ptr, &ptrint);
    mhd_init(&sim.mhd_data, &sim_offload->mhd_offload_data, ptr);

    offload_unpack(offload_data, offload_array,
                   sim_offload->asigma_offload_data.offload_array_length,
                   NULL, 0, &ptr, &ptrint);
    asigma_init(&sim.asigma_data, &sim_offload->asigma_offload_data, ptr);

    /* Offload complete. Reset struct so it can be reused. */
    offload_data->unpack_pos     = 0;
    offload_data->int_unpack_pos = 0;

    diag_init(&sim.diag_data, &sim_offload->diag_offload_data,
              diag_offload_array);

    /**************************************************************************/
    /* 2. Meta data (e.g. random number generator) is initialized.            */
    /*                                                                        */
    /**************************************************************************/
    random_init(&sim.random_data, 0);

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

    print_out(VERBOSE_NORMAL,
              "%s: All fields initialized. Simulation begins, %d threads.\n",
              targetname, omp_get_max_threads());

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

                #pragma omp parallel
                simulate_fo_fixed(&pq, &sim);
            }
            else if(pq.n > 0 && sim.sim_mode == simulate_mode_ml) {

                #pragma omp parallel
                simulate_ml_adaptive(&pq, &sim);
            }
        }

        #pragma omp section
        {
#if VERBOSE > 1
            /* Update progress until simulation is complete.             */
            /* Trim .h5 from filename and replace it with _<QID>.stdout  */
            if(id == 0) {
                char filename[519], outfn[256];
                strcpy(outfn, sim_offload->hdf5_out);
                outfn[strlen(outfn)-3] = '\0';
                sprintf(filename, "%s_%s.stdout", outfn, sim_offload->qid);
                sim_monitor(filename, &pq.n, &pq.finished);
            }
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
    int n_new = 0;
    if(sim.sim_mode == simulate_mode_hybrid) {

        /* Determine the number markers that should be run
         * in fo after previous gc simulation */
        for(int i = 0; i < pq.n; i++) {
            if(pq.p[i]->endcond == endcond_hybrid) {
                /* Check that there was no wall between when moving from
                   gc to fo */
                real w_coll;
                int tile = wall_hit_wall(pq.p[i]->r, pq.p[i]->phi, pq.p[i]->z,
                        pq.p[i]->rprt, pq.p[i]->phiprt, pq.p[i]->zprt,
                                         &sim.wall_data, &w_coll);
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

        #pragma omp parallel sections num_threads(2)
        {
            #pragma omp section
            {
                #pragma omp parallel
                simulate_fo_fixed(&pq, &sim);
            }

            #pragma omp section
            {
#if VERBOSE > 1
                /* Trim .h5 from filename and replace it with _<qid>.stdout */
                if(id == 0) {
                    char filename[519], outfn[256];
                    strcpy(outfn, sim_offload->hdf5_out);
                    outfn[strlen(outfn)-3] = '\0';
                    sprintf(filename, "%s_%s.stdout", outfn, sim_offload->qid);
                    sim_monitor(filename, &pq.n, &pq.finished);
                }
#endif
            }
        }
    }

    /**************************************************************************/
    /* 7. Simulation data is deallocated except for data that is mapped back  */
    /*    to host.                                                            */
    /*                                                                        */
    /**************************************************************************/
    free(pq.p);
    diag_free(&sim.diag_data);

    /**************************************************************************/
    /* 8. Execution returns to host where this function was called.           */
    /*                                                                        */
    /**************************************************************************/
    print_out(VERBOSE_NORMAL, "%s: Simulation complete.\n", targetname);
}

/**
 * @brief Initializes simulation settings
 *
 * This function adjusts simulation settings, e.g. how physics are included,
 * according to the given simulation data. This function should only be called
 * once right after input data has been read.
 *
 * @param sim simulation offload struct which has all fields initialized
 */
void simulate_init_offload(sim_offload_data* sim) {
    if(sim->disable_gctransform) {
        gctransform_setorder(0);
    }
    asigma_extrapolate(sim->enable_atomic==2);
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
    sim->record_mode          = offload_data->record_mode;

    sim->fix_usrdef_use       = offload_data->fix_usrdef_use;
    sim->fix_usrdef_val       = offload_data->fix_usrdef_val;
    sim->fix_gyrodef_nstep    = offload_data->fix_gyrodef_nstep;

    sim->ada_tol_orbfol       = offload_data->ada_tol_orbfol;
    sim->ada_tol_clmbcol      = offload_data->ada_tol_clmbcol;
    sim->ada_max_drho         = offload_data->ada_max_drho;
    sim->ada_max_dphi         = offload_data->ada_max_dphi;

    sim->enable_orbfol        = offload_data->enable_orbfol;
    sim->enable_clmbcol       = offload_data->enable_clmbcol;
    sim->enable_mhd           = offload_data->enable_mhd;
    sim->enable_atomic        = offload_data->enable_atomic;
    sim->disable_gctransform  = offload_data->disable_gctransform;
    sim->disable_energyccoll  = offload_data->disable_energyccoll;
    sim->disable_pitchccoll   = offload_data->disable_pitchccoll;
    sim->disable_gcdiffccoll  = offload_data->disable_gcdiffccoll;
    sim->reverse_time         = offload_data->reverse_time;

    sim->endcond_active       = offload_data->endcond_active;
    sim->endcond_lim_simtime  = offload_data->endcond_lim_simtime;
    sim->endcond_max_mileage  = offload_data->endcond_max_mileage;
    sim->endcond_max_cputime  = offload_data->endcond_max_cputime;
    sim->endcond_min_rho      = offload_data->endcond_min_rho;
    sim->endcond_max_rho      = offload_data->endcond_max_rho;
    sim->endcond_min_ekin     = offload_data->endcond_min_ekin;
    sim->endcond_min_thermal  = offload_data->endcond_min_thermal;
    sim->endcond_max_tororb   = offload_data->endcond_max_tororb;
    sim->endcond_max_polorb   = offload_data->endcond_max_polorb;
    sim->endcond_torandpol    = offload_data->endcond_torandpol;

    sim->bmc_timedependent   = offload_data->bmc_timedependent;
    sim->bmc_orbit_subcycles = offload_data->bmc_orbit_subcycles;
    sim->bmc_timestep        = offload_data->bmc_timestep;
    sim->bmc_tstart          = offload_data->bmc_tstart;
    sim->bmc_tstop           = offload_data->bmc_tstop;
    sim->bmc_mass            = offload_data->bmc_mass;
    sim->bmc_charge          = offload_data->bmc_charge;
    sim->bmc_anum            = offload_data->bmc_anum;
    sim->bmc_znum            = offload_data->bmc_znum;

    mccc_init(&sim->mccc_data, !sim->disable_energyccoll,
              !sim->disable_pitchccoll, !sim->disable_gcdiffccoll);

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
        sleep(A5_PRINTPROGRESSINTERVAL);
    }

    fprintf(f, "Simulation finished.\n");
    fclose(f);
}
