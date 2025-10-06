/**
 * Simulation is initialized and run from here
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
#include "endcond.h"
#include "gctransform.h"
#include "particle.h"
#include "random.h"
#include "simulate.h"
#include <string.h>
#include <unistd.h>
// #include "simulate/mccc/mccc.h"
#include "atomic.h"
#include "bfield.h"
#include "boozer.h"
#include "diag.h"
#include "efield.h"
#include "fusion_source.h"
#include "mathlib.h"
#include "mhd.h"
#include "nbi_source.h"
#include "neutral.h"
#include "plasma.h"
#include "rfof.h"
#include "wall.h"

/**
 * Trace provided markers.
 *
 * This simulates markers using given inputs and options. All different types of
 * simulations are initialized and run via this function.
 *
 * @param nmarkers number of markers to be simulated by this process.
 * @param p markers to be simulated by this process.
 */
void ascot_solve_distribution(int nmarkers, particle_state *p, sim_data *sim)
{

    // Size = NSIMD on CPU and Size = Total number of particles on GPU
    int n_queue_size;
#ifdef GPU
    n_queue_size = n_particles;
#else
    n_queue_size = NSIMD;
#endif

#ifdef GPU
    if (sim->params->simulation_mode != 1)
    {
        print_err("Only GO mode ported to GPU. Please set SIM_MODE=1.");
        exit(1);
    }
    if (sim->params->record_mode)
    {
        print_err("RECORD_MODE=1 not ported to GPU. Please disable it.");
        exit(1);
    }
    if (sim->params->enable_atomic)
    {
        print_err("Atomic not yet ported to GPU. Please set ENABLE_ATOMIC=0.");
        exit(1);
    }
    if (sim->params->enable_mhd)
    {
        print_err("MHD not yet ported to GPU. Please set ENABLE_MHD=0.");
        exit(1);
    }
    if (sim->params->collect_orbit)
    {
        print_err("ENABLE_ORBITWRITE=1 not ported to GPU. Please disable it.");
        exit(1);
    }
    if (sim->params->collect_transport_coefficient)
    {
        print_err("ENABLE_TRANSCOEF=1 not ported to GPU. Please disable it.");
        exit(1);
    }
#endif

    GPU_MAP_TO_DEVICE(sim [0:1])
    B_field_offload(&sim->bfield);
    E_field_offload(&sim->efield);
    plasma_offload(&sim->plasma);
    neutral_offload(&sim->neutral);
    wall_offload(&sim->wall);
    boozer_offload(sim->boozer);
    mhd_offload(&sim->mhd);
    asigma_offload(&sim->atomic);
    random_init(&sim->random_data, 0);

    /* Tabulated data not yet implemented */
    // sim->mccc_data->usetabulated = 0;

    if (sim->params->disable_first_order_gctransformation)
    {
        gctransform_setorder(0);
    }
    asigma_extrapolate(sim->params->enable_atomic == 2);

    particle_queue pq;

    pq.n = 0;
    for (int i = 0; i < nmarkers; i++)
    {
        pq.n++;
    }

    pq.p = (particle_state **)malloc(pq.n * sizeof(particle_state *));
    pq.finished = 0;

    pq.next = 0;
    for (int i = 0; i < nmarkers; i++)
    {
        pq.p[pq.next++] = &p[i];
    }
    pq.next = 0;
    if (pq.n > 0 && (sim->params->simulation_mode == simulate_mode_gc ||
                     sim->params->simulation_mode == simulate_mode_hybrid))
    {
        if (sim->params->enable_adaptive)
        {
            OMP_PARALLEL_CPU_ONLY
            simulate_gc_adaptive(&pq, sim);
        }
        else
        {
            OMP_PARALLEL_CPU_ONLY
            simulate_gc_fixed(&pq, sim);
        }
    }
    else if (pq.n > 0 && sim->params->simulation_mode == simulate_mode_fo)
    {
        OMP_PARALLEL_CPU_ONLY
        simulate_fo_fixed(&pq, sim, n_queue_size);
    }
    else if (pq.n > 0 && sim->params->simulation_mode == simulate_mode_ml)
    {
        simulate_fl_adaptive(&pq, sim);
    }

    int n_new = 0;
    if (sim->params->simulation_mode == simulate_mode_hybrid)
    {

        /* Determine the number markers that should be run
         * in fo after previous gc simulation */
        for (int i = 0; i < pq.n; i++)
        {
            if (pq.p[i]->endcond == endcond_hybrid)
            {
                /* Check that there was no wall between when moving from
                   gc to fo */
                real w_coll;
                int tile = wall_hit_wall(
                    pq.p[i]->r, pq.p[i]->phi, pq.p[i]->z, pq.p[i]->rprt,
                    pq.p[i]->phiprt, pq.p[i]->zprt, &w_coll, &sim->wall);
                if (tile > 0)
                {
                    pq.p[i]->walltile = tile;
                    pq.p[i]->endcond |= endcond_wall;
                }
                else
                {
                    n_new++;
                }
            }
        }
    }
    if (n_new > 0)
    {

        /* Reset hybrid marker end condition */
        for (int i = 0; i < pq.n; i++)
        {
            if (pq.p[i]->endcond & endcond_hybrid)
            {
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

void ascot_solve_fusion(
    sim_data *sim, afsi_data *afsi, size_t n, histogram *prod1,
    histogram *prod2)
{

    random_init(&rdata, time((NULL)));

    real m1, q1, m2, q2, mprod1, qprod1, mprod2, qprod2, Q;
    boschhale_reaction(
        afsi->reaction, &m1, &q1, &m2, &q2, &mprod1, &qprod1, &mprod2, &qprod2,
        &Q);

    int prod_mom_space;
    int i0coord, i1coord, i2coord, p1coord, p2coord;
    if (prod1->axes[0].n)
    {
        i0coord = 0;
        i1coord = 1;
        i2coord = 2;
    }
    else if (prod1->axes[3].n)
    {
        i0coord = 3;
        i1coord = 4;
        i2coord = 1;
    }
    else
    {
        return;
    }
    if (prod1->axes[5].n)
    {
        p1coord = 5;
        p2coord = 6;
        prod_mom_space = PPARPPERP;
    }
    else if (prod1->axes[10].n)
    {
        p1coord = 10;
        p2coord = 11;
        prod_mom_space = EKINXI;
    }
    else
    {
        return;
    }

    real time = 0.0;
#pragma omp parallel for
    for (size_t i0 = 0; i0 < afsi->volshape[0]; i0++)
    {
        real *ppara1 = (real *)malloc(n * sizeof(real));
        real *pperp1 = (real *)malloc(n * sizeof(real));
        real *ppara2 = (real *)malloc(n * sizeof(real));
        real *pperp2 = (real *)malloc(n * sizeof(real));
        real *prod1_p1 = (real *)malloc(n * sizeof(real));
        real *prod1_p2 = (real *)malloc(n * sizeof(real));
        real *prod2_p1 = (real *)malloc(n * sizeof(real));
        real *prod2_p2 = (real *)malloc(n * sizeof(real));

        for (size_t i1 = 0; i1 < afsi->volshape[1]; i1++)
        {
            for (size_t i2 = 0; i2 < afsi->volshape[2]; i2++)
            {
                size_t spatial_index =
                    i0 * afsi->volshape[1] * afsi->volshape[2] +
                    i1 * afsi->volshape[2] + i2;
                real r = afsi->r[spatial_index];
                real z = afsi->z[spatial_index];
                real phi = afsi->phi[spatial_index];
                real vol = afsi->vol[spatial_index];

                real psi, rho[2];
                if (B_field_eval_psi(&psi, r, phi, z, time, &sim->bfield) ||
                    B_field_eval_rho(rho, psi, &sim->bfield))
                {
                    continue;
                }

                real density1, density2;
                afsi_sample_reactant_momenta_2d(
                    sim, afsi, m1, m2, vol, n, i0, i1, i2, r, phi, z, time,
                    rho[0], &density1, ppara1, pperp1, &density2, ppara2,
                    pperp2);
                if (density1 == 0 || density2 == 0)
                {
                    continue;
                }
                for (size_t i = 0; i < n; i++)
                {
                    real vcom2;
                    afsi_compute_product_momenta_2d(
                        i, m1, m2, mprod1, mprod2, Q, prod_mom_space, ppara1,
                        pperp1, ppara2, pperp2, &vcom2, prod1_p1, prod1_p2,
                        prod2_p1, prod2_p2);
                    real E = 0.5 * (m1 * m2) / (m1 + m2) * vcom2;

                    real weight = density1 * density2 * sqrt(vcom2) *
                                  boschhale_sigma(afsi->reaction, E) / n * vol;

                    size_t ip1 = math_bin_index(
                        prod1_p1[i], prod1->axes[p1coord].n,
                        prod1->axes[p1coord].min, prod1->axes[p1coord].max);
                    size_t ip2 = math_bin_index(
                        prod1_p2[i], prod1->axes[p2coord].n,
                        prod1->axes[p2coord].min, prod1->axes[p2coord].max);
                    if (ip1 < prod1->axes[p1coord].n &&
                        ip2 < prod1->axes[p2coord].n)
                    {
                        size_t index = i0 * prod1->strides[i0coord] +
                                       i1 * prod1->strides[i1coord] +
                                       i2 * prod1->strides[i2coord] +
                                       ip1 * prod1->strides[p1coord] +
                                       ip2 * prod1->strides[p2coord];
                        prod1->bins[index] += weight * afsi->mult;
                    }

                    ip1 = math_bin_index(
                        prod2_p1[i], prod2->axes[p1coord].n,
                        prod2->axes[p1coord].min, prod2->axes[p1coord].max);
                    ip2 = math_bin_index(
                        prod2_p2[i], prod2->axes[p2coord].n,
                        prod2->axes[p2coord].min, prod2->axes[p2coord].max);
                    if (ip1 < prod2->axes[p1coord].n &&
                        ip2 < prod2->axes[p2coord].n)
                    {
                        size_t index = i0 * prod2->strides[i0coord] +
                                       i1 * prod2->strides[i1coord] +
                                       i2 * prod2->strides[i2coord] +
                                       ip1 * prod2->strides[p1coord] +
                                       ip2 * prod2->strides[p2coord];
                        prod2->bins[index] += weight * afsi->mult;
                    }
                }
            }
        }
        free(ppara1);
        free(ppara2);
        free(pperp1);
        free(pperp2);
        free(prod1_p1);
        free(prod1_p2);
        free(prod2_p1);
        free(prod2_p2);
    }

    /*
    m1     = m1 / CONST_U;
    m2     = m2 / CONST_U;
    mprod1 = mprod1 / CONST_U;
    mprod2 = mprod2 / CONST_U;
    int c1     = (int)rint(q1 / CONST_E);
    int c2     = (int)rint(q2 / CONST_E);
    int cprod1 = (int)rint(qprod1 / CONST_E);
    int cprod2 = (int)rint(qprod2 / CONST_E);
    */
}

void ascot_solve_nbi(
    sim_data *sim, int nprt, real t1, real t2, particle_state **p)
{

    /* Initialize input data */
    // simulate_init(sim);
    random_init(&sim->random_data, time(NULL));

    /* Calculate total NBI power so that we can distribute markers along
     * the injectors according to their power */
    real total_power = 0;
    for (int i = 0; i < sim->nbi.ninj; i++)
    {
        total_power += sim->nbi.inj[i].power;
    }

    /* Initialize particle struct */
    *p = (particle_state *)malloc(nprt * sizeof(particle_state));

    /* Generate markers at the injectors */
    int nprt_generated = 0;
    for (int i = 0; i < sim->nbi.ninj; i++)
    {

        /* Number of markers generated is proportional to NBI power */
        int nprt_inj = (sim->nbi.inj[i].power / total_power) * nprt;
        if (i == sim->nbi.ninj - 1)
        {
            /* All "remaining" markers goes to the last injector to avoid any
             * rounding issues */
            nprt_inj = nprt - nprt_generated;
        }

        /* Generates markers at the injector location and traces them until
         * they enter the region with magnetic field data */
        bbnbi_inject_markers(
            &((*p)[nprt_generated]), nprt_inj, nprt_generated, t1, t2,
            &(sim->nbi.inj[i]), sim);

        nprt_generated += nprt_inj;
    }

    /* Place markers in a queue */
    particle_queue pq;
    pq.n = 0;
    for (int i = 0; i < nprt; i++)
    {
        pq.n++;
    }
    pq.p = (particle_state **)malloc(pq.n * sizeof(particle_state *));
    pq.finished = 0;

    pq.next = 0;
    for (int i = 0; i < nprt; i++)
    {
        pq.p[pq.next++] = &((*p)[i]);
    }
    pq.next = 0;

/* Trace neutrals until they are ionized or lost to the wall */
#pragma omp parallel
    bbnbi_trace_markers(&pq, sim);
}
