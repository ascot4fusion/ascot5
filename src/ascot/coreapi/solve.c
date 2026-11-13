/**
 * Main simulation routines (see ascot.h).
 */
#include "ascot.h"
#include "data/atomic.h"
#include "data/bfield.h"
#include "data/boozer.h"
#include "data/diag.h"
#include "data/diag_hist.h"
#include "data/efield.h"
#include "data/marker.h"
#include "data/mhd.h"
#include "data/neutral.h"
#include "data/plasma.h"
#include "data/rfof.h"
#include "data/wall.h"
#include "simulate/fusion_source.h"
#include "simulate/nbi_source.h"
#include "simulate/simulate.h"
#include "utils/gctransform.h"
#include "utils/mathlib.h"
#include <signal.h>
#include <stdatomic.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

volatile sig_atomic_t stop_flag = 0;

static void handle_signal(int sig)
{
    (void)sig;
    stop_flag = 1;
}

/**
 * Helper function that allocates marker buffers for fusion source calculations.
 */
static int allocate_fusion_source_buffers(
    size_t nsample, MarkerGyroOrbit *react1, MarkerGyroOrbit *react2,
    MarkerGyroOrbit *prod1, MarkerGyroOrbit *prod2, MarkerGyroOrbit *mrk_copy)
{
    if (MarkerGyroOrbit_allocate(react1, nsample))
        return 1;
    if (MarkerGyroOrbit_allocate(react2, nsample))
        return 1;
    if (MarkerGyroOrbit_allocate(prod1, nsample))
        return 1;
    if (MarkerGyroOrbit_allocate(prod2, nsample))
        return 1;
    if (MarkerGyroOrbit_allocate(mrk_copy, nsample))
        return 1;
    return 0;
}

void ascot_setup_signal_handlers(void)
{
    struct sigaction sa = {0};
    sa.sa_handler = handle_signal;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = 0;

    sigaction(SIGINT, &sa, NULL);
    sigaction(SIGTERM, &sa, NULL);
}

int ascot_solve_distribution(Simulation *sim, size_t nmrk, State mrk[nmrk])
{
    GPU_MAP_TO_DEVICE(sim [0:1])
    Bfield_offload(&sim->bfield);
    Efield_offload(&sim->efield);
    Plasma_offload(&sim->plasma);
    Neutral_offload(&sim->neutral);
    Wall_offload(&sim->wall);
    Boozer_offload(sim->boozer);
    Mhd_offload(&sim->mhd);
    Atomic_offload(sim->atomic);

#ifdef GPU
    size_t vector_size = nmrk;
#else
    size_t vector_size = NSIMD;
#endif
    MarkerQueue queue;
    queue.n = 0;
    for (size_t i = 0; i < nmrk; i++)
        queue.n++;

    int err = 0;
    queue.p = (State **)xmalloc(&err, queue.n * sizeof(State *));
    queue.finished = 0;
    queue.next = 0;
    for (size_t i = 0; i < nmrk; i++)
        queue.p[queue.next++] = &mrk[i];

    queue.next = 0;
    if (queue.n > 0 && (sim->options->mode == simulate_mode_gc ||
                        sim->options->mode == simulate_mode_hybrid))
    {
        if (sim->options->enable_adaptive)
        {
            OMP_PARALLEL_CPU_ONLY
            err = simulate_gc_adaptive(sim, &queue, vector_size);
        }
        else
        {
            OMP_PARALLEL_CPU_ONLY
            err = simulate_gc_fixed(sim, &queue, vector_size);
        }
    }
    else if (queue.n > 0 && sim->options->mode == simulate_mode_fo)
    {
        OMP_PARALLEL_CPU_ONLY
        err = simulate_go_fixed(sim, &queue, vector_size);
    }
    else if (queue.n > 0 && sim->options->mode == simulate_mode_ml)
    {
        err = simulate_fl_adaptive(sim, &queue, vector_size);
    }

    size_t n_new = 0;
    if (sim->options->mode == simulate_mode_hybrid)
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
                size_t tile = Wall_eval_intersection(
                    &w_coll, queue.p[i]->r, queue.p[i]->phi, queue.p[i]->z,
                    queue.p[i]->rprt, queue.p[i]->phiprt, queue.p[i]->zprt,
                    &sim->wall);
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
        err = simulate_go_fixed(sim, &queue, vector_size);
    }
    free(queue.p);
    Diag_onload(&sim->diagnostics);
    return 0;
}

int ascot_solve_fusion(
    Simulation *sim, FusionSource *source, size_t nsample, real time,
    DiagHist *product1, DiagHist *product2)
{

    int failed = 0;
    OMP_PARALLEL_CPU_ONLY
    {
        MarkerGyroOrbit react1, react2, prod1, prod2, mrk_copy;
        int fail = allocate_fusion_source_buffers(
            nsample, &react1, &react2, &prod1, &prod2, &mrk_copy);
        if (fail)
        {
            GPU_ATOMIC
            failed++;
        }

        OMP_BARRIER
        if (failed)
            goto CLEANUP;

        react1.mass = source->mass[0];
        react2.mass = source->mass[1];
        react1.charge[0] = source->charge[0];
        react2.charge[0] = source->charge[1];
        prod1.mass = source->mass[2];
        prod2.mass = source->mass[3];
        prod1.charge[0] = source->charge[2];
        prod2.charge[0] = source->charge[3];
        for (size_t i = 0; i < nsample; i++)
        {
            react1.id[i] = 1 + i;
            react1.time[i] = time;

            prod1.id[i] = 1 + i;
            prod1.time[i] = time;
        }

        CPU_PARALLELIZE_LOOP
        for (size_t j = 0; j < source->nspatial; j++)
        {
            real r = source->r[j];
            real z = source->z[j];
            real phi = source->phi[j];
            real vol = source->vol[j];

            real psi, rho[2];
            if (Bfield_eval_psi(&psi, r, phi, z, time, &sim->bfield) ||
                Bfield_eval_rho(rho, psi, &sim->bfield))
                continue;

            for (size_t i = 0; i < nsample; i++)
            {
                react1.r[i] = r;
                react1.phi[i] = phi;
                react1.z[i] = z;
                react1.rho[i] = rho[0];
                react1.id[i] = 1 + i + j * nsample;
                react1.rng_count[i] = 0;
            }

            real density1, density2;
            if (source->beam1 == NULL)
            {
                size_t species = 0;
                fusion_source_sample_thermal(sim, species, &react1, &density1);
            }
            else
                fusion_source_sample_beam(
                    source->beam1, j, vol, &react1, &density1);
            if (source->beam2 == NULL)
            {
                size_t species = 0;
                fusion_source_sample_thermal(sim, species, &react2, &density2);
            }
            else
                fusion_source_sample_beam(
                    source->beam2, j, vol, &react2, &density2);
            if (density1 == 0 || density2 == 0)
                continue;

            fusion_source_compute_product(
                source->reaction, vol, &react1, &react2, &prod1, &prod2);

            for (size_t i = 0; i < nsample; i++)
            {
                MarkerGyroOrbit_copy(&mrk_copy, &prod1, i);
                mrk_copy.time[i] = time + 1.;
            }
            DiagHist_update_go(product1, &sim->bfield, &prod1, &mrk_copy);

            for (size_t i = 0; i < nsample; i++)
            {
                MarkerGyroOrbit_copy(&mrk_copy, &prod2, i);
                mrk_copy.time[i] = time + 1.;
            }
            DiagHist_update_go(product2, &sim->bfield, &prod2, &mrk_copy);
        }

    CLEANUP:
        MarkerGyroOrbit_deallocate(&react1);
        MarkerGyroOrbit_deallocate(&react2);
        MarkerGyroOrbit_deallocate(&prod1);
        MarkerGyroOrbit_deallocate(&prod2);
        MarkerGyroOrbit_deallocate(&mrk_copy);
    }
    return 0;
}

int ascot_solve_nbi(
    Simulation *sim, size_t ninj, Nbi injectors[ninj], real time_step,
    real max_dist, real tlim[2], size_t nmrk, State mrk[nmrk])
{

    /* Calculate total NBI power so that we can distribute markers along
     * the injectors according to their power */
    real total_power = 0;
    for (size_t i = 0; i < ninj; i++)
        total_power += injectors[i].power;

    /* Generate markers at the injectors */
    size_t nprt_generated = 0;
    for (size_t i = 0; i < ninj; i++)
    {

        /* Number of markers generated is proportional to NBI power. All
         * "remaining" markers goes to the last injector to avoid any rounding
         * issues. */
        size_t nprt_inj = (injectors[i].power / total_power) * nmrk;
        if (i == ninj - 1)
            nprt_inj = nmrk - nprt_generated;

        nbi_source_inject_markers(
            sim, &(injectors[i]), time_step, max_dist, tlim, nprt_inj,
            &(mrk[nprt_generated]));

        nprt_generated += nprt_inj;
    }
    return 0;
}

void ascot_solve_field(
    size_t npnt, size_t ncoil, real coilxyz[3][ncoil], real xyz[3][npnt],
    real bxyz[3][npnt])
{
    real *x = xyz[0];
    real *y = xyz[1];
    real *z = xyz[2];

#pragma omp parallel for
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
            p1[0] = p2[0];
            p1[1] = p2[1];
            p1[2] = p2[2];

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
