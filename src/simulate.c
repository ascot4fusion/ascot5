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
#include "random.h"
#include "simulate.h"
#include "particle.h"
#include "gctransform.h"
#include "simulate/mccc/mccc.h"
#include "simulate/simulate_ml_adaptive.h"
#include "simulate/simulate_gc_adaptive.h"
#include "simulate/simulate_gc_fixed.h"
#include "simulate/simulate_fo_fixed.h"
#include "B_field.h"
#include "E_field.h"
#include "plasma.h"
#include "neutral.h"
#include "wall.h"
#include "boozer.h"
#include "mhd.h"
#include "diag.h"
#include "asigma.h"
#include "rfof.h"


/**
 * @brief Simulate given markers
 *
 * This simulates markers using given inputs and options. All different types of
 * simulations are initialized and run via this function.
 *
 * @param nmarkers number of markers to be simulated by this process.
 * @param p markers to be simulated by this process.
 * @param bfield magnetic field data.
 * @param efield electric field data.
 * @param plasma plasma data.
 * @param neutral neutral data.
 * @param wall wall data.
 * @param boozer boozer data.
 * @param mhd mhd data.
 * @param atomic atomic data.
 */
void simulate(int nmarkers, particle_state* p, sim_data* sim) {

    // Size = NSIMD on CPU and Size = Total number of particles on GPU
    int n_queue_size;
#ifdef GPU
    n_queue_size = n_particles;
#else
    n_queue_size = NSIMD;
#endif

#ifdef GPU
    if(sim->params->simulation_mode != 1) {
        print_err("Only GO mode ported to GPU. Please set SIM_MODE=1.");
        exit(1);
    }
    if(sim->params->record_mode) {
        print_err("RECORD_MODE=1 not ported to GPU. Please disable it.");
        exit(1);
    }
    if(sim->params->enable_atomic) {
        print_err("Atomic not yet ported to GPU. Please set ENABLE_ATOMIC=0.");
        exit(1);
    }
    if(sim->params->enable_mhd) {
        print_err("MHD not yet ported to GPU. Please set ENABLE_MHD=0.");
        exit(1);
    }
    if(sim->params->collect_orbit) {
        print_err(
            "ENABLE_ORBITWRITE=1 not ported to GPU. Please disable it.");
        exit(1);
    }
    if(sim->params->collect_transport_coefficient) {
        print_err(
            "ENABLE_TRANSCOEF=1 not ported to GPU. Please disable it.");
        exit(1);
    }
#endif

    GPU_MAP_TO_DEVICE(sim[0:1])
    B_field_offload(&sim->B_data);
    E_field_offload(&sim->E_data);
    plasma_offload(&sim->plasma_data);
    neutral_offload(&sim->neutral_data);
    wall_offload(&sim->wall_data);
    boozer_offload(sim->boozer_data);
    mhd_offload(&sim->mhd_data);
    asigma_offload(&sim->asigma_data);
    random_init(&sim->random_data, 0);

    /* Tabulated data not yet implemented */
    //sim->mccc_data->usetabulated = 0;

    if(sim->params->disable_first_order_gctransformation) {
        gctransform_setorder(0);
    }
    asigma_extrapolate(sim->params->enable_atomic==2);

    particle_queue pq;

    pq.n = 0;
    for(int i = 0; i < nmarkers; i++) {
        pq.n++;
    }

    pq.p = (particle_state**) malloc(pq.n * sizeof(particle_state*));
    pq.finished = 0;

    pq.next = 0;
    for(int i = 0; i < nmarkers; i++) {
        pq.p[pq.next++] = &p[i];

    }
    pq.next = 0;
    if(pq.n > 0 &&
        (sim->params->simulation_mode == simulate_mode_gc ||
         sim->params->simulation_mode == simulate_mode_hybrid)) {
        if(sim->params->enable_adaptive) {
            OMP_PARALLEL_CPU_ONLY
            simulate_gc_adaptive(&pq, sim);
        }
        else {
            OMP_PARALLEL_CPU_ONLY
            simulate_gc_fixed(&pq, sim);
        }
    }
    else if(pq.n > 0 && sim->params->simulation_mode == simulate_mode_fo) {
        OMP_PARALLEL_CPU_ONLY
        simulate_fo_fixed(&pq, sim, n_queue_size);
    }
    else if(pq.n > 0 && sim->params->simulation_mode == simulate_mode_ml) {
        simulate_ml_adaptive(&pq, sim);
    }

    int n_new = 0;
    if(sim->params->simulation_mode == simulate_mode_hybrid) {

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
        OMP_PARALLEL_CPU_ONLY
        simulate_fo_fixed(&pq, sim, n_queue_size);
    }
    free(pq.p);
}
