/**
 * Implements nbi_source.h.
 */
#include "nbi_source.h"
#include "ascot.h"
#include "data/atomic.h"
#include "data/bfield.h"
#include "consts.h"
#include "defines.h"
#include "data/diag.h"
#include "endcond.h"
#include "data/marker.h"
#include "utils/mathlib.h"
#include "data/nbi.h"
#include "data/neutral.h"
#include "parallel.h"
#include "utils/physlib.h"
#include "data/plasma.h"
#include "utils/random.h"
#include "datatypes.h"
#include "utils/suzuki.h"
#include "data/wall.h"
#include <getopt.h>
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

void nbi_source_inject_markers(
    Simulation *sim, Nbi *inj, real tlim[2], size_t nmrk, State *mrk)
{

    /* Set marker weights assuming a large number is created so that the energy
     * fractions of generated markers are close to the injector values */
    real f = 1.0 * inj->efrac[0] + (1.0 / 2) * inj->efrac[1] +
             (1.0 / 3) * inj->efrac[2];
    real weight = (inj->power / inj->energy) / (f * nmrk);

    OMP_PARALLEL_CPU_ONLY
    for (size_t i = 0; i < nmrk; i++)
    {
        real time =
            tlim[0] + random_uniform(sim->random_data) * (tlim[1] - tlim[0]);

        real xyz[3], vxyz[3], rpz[3], vhat[3];
        nbi_inject(xyz, vxyz, inj, sim->random_data);
        math_xyz2rpz(xyz, rpz);
        math_unit(vxyz, vhat);

        real psi;
        real ds = 1e-3;
        err_t err =
            Bfield_eval_psi(&psi, rpz[0], rpz[1], rpz[2], time, &sim->bfield);
        while (err)
        {
            xyz[0] += ds * vhat[0];
            xyz[1] += ds * vhat[1];
            xyz[2] += ds * vhat[2];
            math_xyz2rpz(xyz, rpz);
            err = Bfield_eval_psi(
                &psi, rpz[0], rpz[1], rpz[2], time, &sim->bfield);
        }

        real vrpz[3];
        math_vec_xyz2rpz(vxyz, vrpz, rpz[1]);
        real gamma = physlib_gamma_vnorm(math_norm(vrpz));

        mrk[i].err = 0;
        mrk[i].rprt = rpz[0];
        mrk[i].zprt = rpz[2];
        mrk[i].phiprt = rpz[1];
        mrk[i].pr = vrpz[0] * gamma * inj->mass;
        mrk[i].pz = vrpz[2] * gamma * inj->mass;
        mrk[i].pphi = vrpz[1] * gamma * inj->mass;
        mrk[i].mass = inj->mass;
        mrk[i].anum = inj->anum;
        mrk[i].znum = inj->znum;
        mrk[i].time = time;
        mrk[i].weight = weight;
        mrk[i].charge = 0.0;
        mrk[i].mileage = 0.0;
        mrk[i].cputime = 0.0;
        mrk[i].endcond = 0;
        mrk[i].walltile = 0;
    }
}

void nbi_source_trace_markers(MarkerQueue *pq, Simulation *sim)
{
    size_t cycle[NSIMD] __memalign__;
    real hin[NSIMD] __memalign__;
    int shinethrough[NSIMD] __memalign__;
    real remaining[NSIMD] __memalign__;
    real threshold[NSIMD] __memalign__;
    MarkerGyroOrbit p, p0, pdiag;

    int n_species = Plasma_get_n_species(&sim->plasma);
    const int *pls_anum = Plasma_get_species_anum(&sim->plasma);
    const int *pls_znum = Plasma_get_species_znum(&sim->plasma);

    /* Init dummy markers */
    marker_allocate_go(&p, NSIMD);
    marker_allocate_go(&p0, NSIMD);
    marker_allocate_go(&pdiag, NSIMD);
    for (size_t i = 0; i < NSIMD; i++)
    {
        p.id[i] = -1;
        p.running[i] = 0;
        hin[i] = 1e-10;
        threshold[i] = random_uniform(&sim->random_data);
        remaining[i] = 1.0;
        shinethrough[i] = 0;
    }

    /* Initialize running particles */
    size_t n_running = marker_cycle_go(pq, &p, &sim->bfield, cycle);
    while (n_running > 0)
    {

        OMP_PARALLEL_CPU_ONLY
        for (size_t i = 0; i < NSIMD; i++)
        {
            /* Store marker states */
            marker_copy_go(&p, i, &p0, i);

            if (p.running[i])
            {
                err_t err = 0;

                /* These are needed later */
                real pnorm = math_normc(p.p_r[i], p.p_phi[i], p.p_z[i]);
                real gamma = physlib_gamma_pnorm(p.mass[i], pnorm);
                real ekin = physlib_Ekin_pnorm(p.mass[i], pnorm);

                /* Advance ballistic trajectory by converting momentum to
                 * cartesian coordinates */
                real prpz[3] = {p.p_r[i], p.p_phi[i], p.p_z[i]};
                real pxyz[3];
                math_vec_rpz2xyz(prpz, pxyz, p.phi[i]);

                real posrpz[3] = {p.r[i], p.phi[i], p.z[i]};
                real posxyz[3], fposxyz[3];
                math_rpz2xyz(posrpz, posxyz);
                fposxyz[0] = posxyz[0] + pxyz[0] * hin[i] / (gamma * p.mass[i]);
                fposxyz[1] = posxyz[1] + pxyz[1] * hin[i] / (gamma * p.mass[i]);
                fposxyz[2] = posxyz[2] + pxyz[2] * hin[i] / (gamma * p.mass[i]);

                /* Back to cylindrical coordinates (note phi is cumulative) */
                p.r[i] =
                    sqrt(fposxyz[0] * fposxyz[0] + fposxyz[1] * fposxyz[1]);
                p.phi[i] += atan2(
                    posxyz[0] * fposxyz[1] - posxyz[1] * fposxyz[0],
                    posxyz[0] * fposxyz[0] + posxyz[1] * fposxyz[1]);
                p.z[i] = fposxyz[2];

                real cosp = cos(p.phi[i]);
                real sinp = sin(p.phi[i]);
                p.p_r[i] = pxyz[0] * cosp + pxyz[1] * sinp;
                p.p_phi[i] = -pxyz[0] * sinp + pxyz[1] * cosp;
                p.p_z[i] = pxyz[2];

                real ds = hin[i];
                p.mileage[i] += hin[i];

                /* Update background values at the new position */
                real psi, rho[2], pls_dens[MAX_SPECIES], pls_temp[MAX_SPECIES];
                err = Bfield_eval_psi(
                    &psi, p.r[i], p.phi[i], p.z[i], p.time[i], &sim->bfield);
                if (!err)
                {
                    err = Bfield_eval_rho(rho, psi, &sim->bfield);
                }
                if (!err && p.rho[i] <= 1.0 && rho[0] > 1.0)
                {
                    shinethrough[i] = 1;
                }
                p.rho[i] = rho[0];
                /* Update theta value */
                real axisrz[2];
                if (!err)
                {
                    Bfield_eval_axis_rz(axisrz, &sim->bfield, p.phi[i]);
                }
                p.theta[i] = atan2(p.z[i] - axisrz[1], p.r[i] - axisrz[0]);

                if (!err)
                {
                    err = Plasma_eval_densandtemp(
                        pls_dens, pls_temp, rho[0], p.r[i], p.phi[i], p.z[i],
                        p.time[i], &sim->plasma);
                }

                /* Calculate ionization rate */
                real rate = 0.0;
                if (!err)
                {
                    real sigmav;
                    if (asigma_eval_bms(
                            &sigmav, p.znum[i], p.anum[i], ekin, p.mass[i],
                            n_species - 1, pls_znum, pls_anum, pls_temp[0],
                            &(pls_dens[1]), sim->atomic))
                    {
                        err = 1;
                    }
                    rate = pls_dens[0] * sigmav;
                }
                remaining[i] *= exp(-rate * ds);

                /* Check for end conditions */
                if (!err)
                {
                    real w_coll = 0;
                    int tile = 0;
                    if (shinethrough[i])
                    {
                        tile = Wall_hit_wall(
                            p0.r[i], p0.phi[i], p0.z[i], p.r[i], p.phi[i],
                            p.z[i], &w_coll, &sim->wall);
                    }
                    if (tile > 0)
                    {
                        real w = w_coll;
                        p.time[i] = p0.time[i] + w * (p.time[i] - p0.time[i]);
                        p.r[i] = p0.r[i] + w * (p.r[i] - p0.r[i]);
                        p.phi[i] = p0.phi[i] + w * (p.phi[i] - p0.phi[i]);
                        p.z[i] = p0.z[i] + w * (p.z[i] - p0.z[i]);

                        p.walltile[i] = tile;
                        p.endcond[i] |= ENDCOND_WALL;
                        p.running[i] = 0;
                    }
                    if (p.mileage[i] > NBI_MAX_DISTANCE)
                    {
                        p.endcond[i] |= ENDCOND_TLIM;
                        p.running[i] = 0;
                    }
                    if (remaining[i] < threshold[i])
                    {
                        p.charge[i] = 1 * CONST_E;
                        p.endcond[i] |= ENDCOND_NEUTR;
                        p.running[i] = 0;
                    }
                }
                else
                {
                    p.err[i] = err;
                    p.running[i] = 0;
                }
            }
        }

        /* Update markers that just finished */
        OMP_PARALLEL_CPU_ONLY
        for (size_t i = 0; i < NSIMD; i++)
        {
            /* Use this as a flag for which markers to update in diagnostics */
            pdiag.running[i] = 0;
            if (!p.running[i] && p.id[i] > 0)
            {
                p.time[i] += p.mileage[i];

                /* Reset these for the next marker */
                threshold[i] = random_uniform(&sim->random_data);
                remaining[i] = 1.0;
                shinethrough[i] = 0;

                /* Update the magnetic field at the marker position */
                if (!p.err[i])
                {
                    real B_dB[15];
                    Bfield_eval_b_db(
                        B_dB, p.r[i], p.phi[i], p.z[i], p.time[i],
                        &sim->bfield);
                    p.B_r[i] = B_dB[0];
                    p.B_r_dr[i] = B_dB[1];
                    p.B_r_dphi[i] = B_dB[2];
                    p.B_r_dz[i] = B_dB[3];

                    p.B_phi[i] = B_dB[4];
                    p.B_phi_dr[i] = B_dB[5];
                    p.B_phi_dphi[i] = B_dB[6];
                    p.B_phi_dz[i] = B_dB[7];

                    p.B_z[i] = B_dB[8];
                    p.B_z_dr[i] = B_dB[9];
                    p.B_z_dphi[i] = B_dB[10];
                    p.B_z_dz[i] = B_dB[11];
                }
            }
            marker_copy_go(&p, i, &pdiag, i);
            if (!p.running[i] && p.id[i] > 0)
            {
                pdiag.running[i] = 1;
            }

            /* Normalize weight with time and add hin so that we don't divide
             * with zero when updating distributions */
            pdiag.time[i] += hin[i];
            pdiag.weight[i] /= hin[i];
        }

        /* Update distributions for markers that finished */
        // diag_update_go(&sim->diag_data, &sim->bfield, &pdiag, &p);

        /* Update running particles */
        n_running = marker_cycle_go(pq, &p, &sim->bfield, cycle);
    }
}
