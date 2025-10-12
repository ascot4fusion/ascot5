/**
 * Main simulation routines (see ascot.h).
 */
#include "ascot.h"
#include "data/atomic.h"
#include "data/bfield.h"
#include "data/boozer.h"
#include "data/diag.h"
#include "data/efield.h"
#include "data/marker.h"
#include "data/mhd.h"
#include "data/neutral.h"
#include "data/plasma.h"
#include "data/rfof.h"
#include "data/wall.h"
#include "endcond.h"
#include "simulate/fusion_source.h"
#include "simulate/nbi_source.h"
#include "simulate/simulate.h"
#include "utils/gctransform.h"
#include "utils/mathlib.h"
#include "utils/random.h"
#include <string.h>
#include <unistd.h>

void ascot_solve_distribution(Simulation *sim, size_t nmrk, State mrk[nmrk])
{

#ifdef GPU
    if (sim->options->simulation_mode != 1)
    {
        print_err("Only GO mode ported to GPU. Please set SIM_MODE=1.");
        exit(1);
    }
    if (sim->options->record_mode)
    {
        print_err("RECORD_MODE=1 not ported to GPU. Please disable it.");
        exit(1);
    }
    if (sim->options->enable_atomic)
    {
        print_err("Atomic not yet ported to GPU. Please set ENABLE_ATOMIC=0.");
        exit(1);
    }
    if (sim->options->enable_mhd)
    {
        print_err("MHD not yet ported to GPU. Please set ENABLE_MHD=0.");
        exit(1);
    }
    if (sim->options->collect_orbit)
    {
        print_err("ENABLE_ORBITWRITE=1 not ported to GPU. Please disable it.");
        exit(1);
    }
    if (sim->options->collect_transport_coefficient)
    {
        print_err("ENABLE_TRANSCOEF=1 not ported to GPU. Please disable it.");
        exit(1);
    }
#endif

    GPU_MAP_TO_DEVICE(sim [0:1])
    Bfield_offload(&sim->bfield);
    Efield_offload(&sim->efield);
    plasma_offload(&sim->plasma);
    neutral_offload(&sim->neutral);
    wall_offload(&sim->wall);
    boozer_offload(sim->boozer);
    mhd_offload(&sim->mhd);
    asigma_offload(sim->atomic);
    random_init(&sim->random_data, 0);

    if (sim->options->disable_first_order_gctransformation)
    {
        gctransform_setorder(0);
    }
    asigma_extrapolate(sim->options->enable_atomic == 2);

#ifdef GPU
    size_t vector_size = nmrk;
#else
    size_t vector_size = NSIMD;
#endif

    MarkerQueue queue;
    queue.n = 0;
    for (size_t i = 0; i < nmrk; i++)
    {
        queue.n++;
    }

    queue.p = (State **)malloc(queue.n * sizeof(State *));
    queue.finished = 0;

    queue.next = 0;
    for (size_t i = 0; i < nmrk; i++)
    {
        queue.p[queue.next++] = &mrk[i];
    }
    queue.next = 0;
    if (queue.n > 0 && (sim->options->simulation_mode == simulate_mode_gc ||
                        sim->options->simulation_mode == simulate_mode_hybrid))
    {
        if (sim->options->enable_adaptive)
        {
            OMP_PARALLEL_CPU_ONLY
            simulate_gc_adaptive(sim, &queue, vector_size);
        }
        else
        {
            OMP_PARALLEL_CPU_ONLY
            simulate_gc_fixed(sim, &queue, vector_size);
        }
    }
    else if (queue.n > 0 && sim->options->simulation_mode == simulate_mode_fo)
    {
        OMP_PARALLEL_CPU_ONLY
        simulate_go_fixed(sim, &queue, vector_size);
    }
    else if (queue.n > 0 && sim->options->simulation_mode == simulate_mode_ml)
    {
        simulate_fl_adaptive(sim, &queue, vector_size);
    }

    size_t n_new = 0;
    if (sim->options->simulation_mode == simulate_mode_hybrid)
    {

        /* Determine the number markers that should be run
         * in go after previous gc simulation */
        for (size_t i = 0; i < queue.n; i++)
        {
            if (queue.p[i]->endcond == ENDCOND_HYBRID)
            {
                /* Check that there was no wall between when moving from
                   gc to fo */
                real w_coll;
                size_t tile = Wall_hit_wall(
                    queue.p[i]->r, queue.p[i]->phi, queue.p[i]->z,
                    queue.p[i]->rprt, queue.p[i]->phiprt, queue.p[i]->zprt,
                    &w_coll, &sim->wall);
                if (tile > 0)
                {
                    queue.p[i]->walltile = tile;
                    queue.p[i]->endcond |= ENDCOND_WALL;
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
        for (size_t i = 0; i < queue.n; i++)
        {
            if (queue.p[i]->endcond & ENDCOND_HYBRID)
            {
                queue.p[i]->endcond ^= ENDCOND_HYBRID;
            }
        }
        queue.next = 0;
        queue.finished = 0;
        OMP_PARALLEL_CPU_ONLY
        simulate_go_fixed(sim, &queue, vector_size);
    }
    free(queue.p);
}

void ascot_solve_fusion(
    Simulation *sim, FusionSource *source, size_t nsample, histogram *product1,
    histogram *product2)
{

    random_init(&rdata, time((NULL)));

    real m1, q1, m2, q2, mprod1, qprod1, mprod2, qprod2, Q;
    boschhale_reaction(
        source->reaction, &m1, &q1, &m2, &q2, &mprod1, &qprod1, &mprod2,
        &qprod2, &Q);

    int prod_mom_space;
    int i0coord, i1coord, i2coord, p1coord, p2coord;
    if (product1->axes[0].n)
    {
        i0coord = 0;
        i1coord = 1;
        i2coord = 2;
    }
    else if (product1->axes[3].n)
    {
        i0coord = 3;
        i1coord = 4;
        i2coord = 1;
    }
    else
    {
        return;
    }
    if (product1->axes[5].n)
    {
        p1coord = 5;
        p2coord = 6;
        prod_mom_space = PPARPPERP;
    }
    else if (product1->axes[10].n)
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
    OMP_PARALLEL_CPU_ONLY
    for (size_t i0 = 0; i0 < source->volshape[0]; i0++)
    {
        real *ppara1 = (real *)malloc(nsample * sizeof(real));
        real *pperp1 = (real *)malloc(nsample * sizeof(real));
        real *ppara2 = (real *)malloc(nsample * sizeof(real));
        real *pperp2 = (real *)malloc(nsample * sizeof(real));
        real *prod1_p1 = (real *)malloc(nsample * sizeof(real));
        real *prod1_p2 = (real *)malloc(nsample * sizeof(real));
        real *prod2_p1 = (real *)malloc(nsample * sizeof(real));
        real *prod2_p2 = (real *)malloc(nsample * sizeof(real));

        for (size_t i1 = 0; i1 < source->volshape[1]; i1++)
        {
            for (size_t i2 = 0; i2 < source->volshape[2]; i2++)
            {
                size_t spatial_index =
                    i0 * source->volshape[1] * source->volshape[2] +
                    i1 * source->volshape[2] + i2;
                real r = source->r[spatial_index];
                real z = source->z[spatial_index];
                real phi = source->phi[spatial_index];
                real vol = source->vol[spatial_index];

                real psi, rho[2];
                if (Bfield_eval_psi(&psi, r, phi, z, time, &sim->bfield) ||
                    Bfield_eval_rho(rho, psi, &sim->bfield))
                {
                    continue;
                }

                real density1, density2;
                afsi_sample_reactant_momenta_2d(
                    sim, source, m1, m2, vol, nsample, i0, i1, i2, r, phi, z,
                    time, rho[0], &density1, ppara1, pperp1, &density2, ppara2,
                    pperp2);
                if (density1 == 0 || density2 == 0)
                {
                    continue;
                }
                for (size_t i = 0; i < nsample; i++)
                {
                    real vcom2;
                    afsi_compute_product_momenta_2d(
                        i, m1, m2, mprod1, mprod2, Q, prod_mom_space, ppara1,
                        pperp1, ppara2, pperp2, &vcom2, prod1_p1, prod1_p2,
                        prod2_p1, prod2_p2);
                    real E = 0.5 * (m1 * m2) / (m1 + m2) * vcom2;

                    real weight = density1 * density2 * sqrt(vcom2) *
                                  boschhale_sigma(source->reaction, E) /
                                  nsample * vol;

                    size_t ip1 = math_bin_index(
                        prod1_p1[i], product1->axes[p1coord].n,
                        product1->axes[p1coord].min,
                        product1->axes[p1coord].max);
                    size_t ip2 = math_bin_index(
                        prod1_p2[i], product1->axes[p2coord].n,
                        product1->axes[p2coord].min,
                        product1->axes[p2coord].max);
                    if (ip1 < product1->axes[p1coord].n &&
                        ip2 < product1->axes[p2coord].n)
                    {
                        size_t index = i0 * product1->strides[i0coord] +
                                       i1 * product1->strides[i1coord] +
                                       i2 * product1->strides[i2coord] +
                                       ip1 * product1->strides[p1coord] +
                                       ip2 * product1->strides[p2coord];
                        product1->bins[index] += weight * source->mult;
                    }

                    ip1 = math_bin_index(
                        prod2_p1[i], product2->axes[p1coord].n,
                        product2->axes[p1coord].min,
                        product2->axes[p1coord].max);
                    ip2 = math_bin_index(
                        prod2_p2[i], product2->axes[p2coord].n,
                        product2->axes[p2coord].min,
                        product2->axes[p2coord].max);
                    if (ip1 < product2->axes[p1coord].n &&
                        ip2 < product2->axes[p2coord].n)
                    {
                        size_t index = i0 * product2->strides[i0coord] +
                                       i1 * product2->strides[i1coord] +
                                       i2 * product2->strides[i2coord] +
                                       ip1 * product2->strides[p1coord] +
                                       ip2 * product2->strides[p2coord];
                        product2->bins[index] += weight * source->mult;
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
    Simulation *sim, size_t ninj, Nbi injectors[ninj], real tlim[2],
    size_t nmrk, State mrk[nmrk])
{

    /* Initialize input data */
    random_init(&sim->random_data, time(NULL));

    /* Calculate total NBI power so that we can distribute markers along
     * the injectors according to their power */
    real total_power = 0;
    for (size_t i = 0; i < ninj; i++)
    {
        total_power += injectors[i].power;
    }

    /* Generate markers at the injectors */
    size_t nprt_generated = 0;
    for (size_t i = 0; i < ninj; i++)
    {

        /* Number of markers generated is proportional to NBI power */
        size_t nprt_inj = (injectors[i].power / total_power) * nmrk;
        if (i == ninj - 1)
        {
            /* All "remaining" markers goes to the last injector to avoid any
             * rounding issues */
            nprt_inj = nmrk - nprt_generated;
        }

        nbi_source_inject_markers(
            sim, &(injectors[i]), tlim, nprt_inj, &(mrk[nprt_generated]));

        nprt_generated += nprt_inj;
    }

    MarkerQueue queue;
    queue.n = 0;
    for (size_t i = 0; i < nmrk; i++)
    {
        mrk[i].id = i + 1;
        queue.n++;
    }
    queue.p = (State **)malloc(queue.n * sizeof(State *));
    queue.finished = 0;

    queue.next = 0;
    for (size_t i = 0; i < nmrk; i++)
    {
        queue.p[queue.next++] = &(mrk[i]);
    }
    queue.next = 0;

    /* Trace neutrals until they are ionized or lost to the wall */
    OMP_PARALLEL_CPU_ONLY
    nbi_source_trace_markers(&queue, sim);
}

void ascot_solve_field(
    size_t npnt, size_t ncoil, real coilxyz[3][ncoil], real xyz[3][npnt],
    real bxyz[3][npnt])
{
    real *x = xyz[0];
    real *y = xyz[1];
    real *z = xyz[2];

    OMP_PARALLEL_CPU_ONLY
    for (size_t ix = 0; ix < npnt; ix++)
    {
        bxyz[0][ix] = 0;
        bxyz[1][ix] = 0;
        bxyz[2][ix] = 0;

        real p1[3], p2[3];
        p2[0] = coilxyz[0][0];
        p2[1] = coilxyz[1][0];
        p2[2] = coilxyz[2][0];

        for (size_t i = 1; i < ncoil; i++)
        {
            math_copy(p1, p2);

            p2[0] = coilxyz[0][i * 3];
            p2[1] = coilxyz[1][i * 3];
            p2[2] = coilxyz[2][i * 3];

            real p1p2[3];
            p1p2[0] = p2[0] - p1[0];
            p1p2[1] = p2[1] - p1[1];
            p1p2[2] = p2[2] - p1[2];

            real p1x[3];
            p1x[0] = x[ix] - p1[0];
            p1x[1] = y[ix] - p1[1];
            p1x[2] = z[ix] - p1[2];

            real p2x[3];
            p2x[0] = x[ix] - p2[0];
            p2x[1] = y[ix] - p2[1];
            p2x[2] = z[ix] - p2[2];

            real d1 = math_norm(p1x);
            real d2 = math_norm(p2x);
            real l = math_norm(p1p2);
            real s = math_dot(p1p2, p1x) / math_dot(p1p2, p1p2);
            real h = s * l;

            real xs[3];
            xs[0] = p1[0] + s * p1p2[0] - x[ix];
            xs[1] = p1[1] + s * p1p2[1] - y[ix];
            xs[2] = p1[2] + s * p1p2[2] - z[ix];

            real d = math_norm(xs);

            real bnorm =
                CONST_MU0 / (4 * CONST_PI) * ((l - h) / d2 + h / d1) / d;

            real bdir[3];
            math_cross(p1x, p2x, bdir);

            bxyz[0][ix] += bnorm * bdir[0] / math_norm(bdir);
            bxyz[1][ix] += bnorm * bdir[1] / math_norm(bdir);
            bxyz[2][ix] += bnorm * bdir[2] / math_norm(bdir);
        }
    }
}
