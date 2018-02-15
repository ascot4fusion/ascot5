#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "endcond.h"
#include "hdf5io/hdf5_orbits.h"
#include "offload.h"
#include "particle.h"
#include "simulate.h"
#include "simulate/simulate_ml_adaptive.h"
#include "simulate/simulate_gc_adaptive.h"
#include "simulate/simulate_gc_fixed.h"
#include "simulate/simulate_fo_fixed.h"

#pragma omp declare target
void sim_init(sim_data* sim, sim_offload_data* offload_data);
void sim_monitor(FILE* f, int* n, int finished);
#pragma omp end declare target

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
    plasma_1d_init(&sim.plasma_data, &sim_offload->plasma_offload_data, ptr);

    ptr = offload_unpack(offload_data, offload_array,
            sim_offload->wall_offload_data.offload_array_length);
    wall_init(&sim.wall_data, &sim_offload->wall_offload_data, ptr);

    diag_init(&sim.diag_data, &sim_offload->diag_offload_data,
            diag_offload_array);

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

    /* Spawn two threads: one for the actual simulation and one worker
     * for updating process. */
    #pragma omp parallel sections num_threads(2) 
    {
        #pragma omp section
        {
        /* Choose simulation mode and spawn more threads that each begins
           simulation. */
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
            sim_monitor(f, &pq.n, pq.finished);
#endif
        }
    }

    /* Finish simulating hybrid particles with fo */
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
                sim_monitor(f, &pq_hybrid.n, pq_hybrid.finished);
#endif
            }
        }
    }

#if VERBOSE > 1
    /* Close progress file*/
    fclose(f);
#endif

    free(pq.p);
    free(pq_hybrid.p);

    // Temporary solution
#ifdef NOTARGET
    hdf5_orbits_write(&sim, sim_offload->hdf5_out);
#endif

    diag_clean(&sim.diag_data);

    #if VERBOSE >= 1
    printf("%s: Simulation complete.\n", targetname);
    #endif
}

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

void sim_monitor(FILE* f, int* n, int finished) {
    real timer = A5_WTIME;
    while(f != NULL && *n > finished) {
        real epsilon = 1e-10;// To avoid division by zero
        real fracprog = ((real) finished)/(*n);
        real timespent = (A5_WTIME)-timer;
        fprintf(f, "Progress: %d/%d, %.2f %%. Time spent: %.2f h, "
            "estimated time to finish: %.2f h\n", finished, *n, 100*fracprog,
            timespent/3600, (1/(fracprog+epsilon)-1)*timespent/3600);
        fflush(f);
        sleep(A5_PRINTPROGRESSINTERVAL);
    }
}
