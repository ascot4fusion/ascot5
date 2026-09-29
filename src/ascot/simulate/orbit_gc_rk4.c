/**
 * Guiding center integrator implemented with RK4 with fixed time-step
 * (see "orbit_following.h").
 **/
#include "consts.h"
#include "data/bfield.h"
#include "data/boozer.h"
#include "data/efield.h"
#include "data/marker.h"
#include "data/mhd.h"
#include "defines.h"
#include "orbit_following.h"
#include "utils/mathlib.h"
#include "utils/physlib.h"
#include <math.h>
#include <stdio.h>

void step_gc_rk4(
    MarkerGuidingCenter *mrk, const real *h, Bfield *bfield, Efield *efield,
    int aldforce)
{
    GPU_DATA_IS_MAPPED(h [0:mrk->size])
    GPU_PARALLEL_LOOP_ALL_LEVELS
    for (size_t i = 0; i < mrk->size; i++)
    {
        if (mrk->running[i])
        {
            err_t errflag = 0;

            real k1[6], k2[6], k3[6], k4[6];
            real tempy[6];
            real yprev[6];
            real y[6];

            real mass = mrk->mass;
            real charge = mrk->charge[i] * CONST_E;

            real b_db[15], E[3];
            real alpha[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
            real Phi[5] = {0.0, 0.0, 0.0, 0.0, 0.0};

            real R0 = mrk->r[i];
            real z0 = mrk->z[i];
            real t0 = mrk->time[i];

            /* Coordinates are copied from the struct into an array to make
             * passing parameters easier */
            yprev[0] = mrk->r[i];
            yprev[1] = mrk->phi[i];
            yprev[2] = mrk->z[i];
            yprev[3] = mrk->ppar[i];
            yprev[4] = mrk->mu[i];
            yprev[5] = mrk->zeta[i];

            /* Magnetic field at initial position already known */
            b_db[0] = mrk->br[i];
            b_db[3] = mrk->dbrdr[i];
            b_db[4] = mrk->dbrdphi[i];
            b_db[5] = mrk->dbrdz[i];

            b_db[1] = mrk->bphi[i];
            b_db[6] = mrk->dbphidr[i];
            b_db[7] = mrk->dbphidphi[i];
            b_db[8] = mrk->dbphidz[i];

            b_db[2] = mrk->bz[i];
            b_db[9] = mrk->dbzdr[i];
            b_db[10] = mrk->dbzdphi[i];
            b_db[11] = mrk->dbzdz[i];

            if (!errflag)
            {
                errflag = Efield_eval_e(
                    E, yprev[0], yprev[1], yprev[2], t0, efield, bfield);
            }
            if (!errflag)
            {
                step_gceom(
                    k1, yprev, mass, charge, b_db, E, alpha, Phi, aldforce);
            }

            /* particle coordinates for the subsequent ydot evaluations are
             * stored in tempy */
            for (int j = 0; j < 6; j++)
            {
                tempy[j] = yprev[j] + h[i] * k1[j] / 2.0;
            }

            if (!errflag)
            {
                errflag = Bfield_eval_b_db(
                    b_db, tempy[0], tempy[1], tempy[2], t0 + h[i] / 2.0,
                    bfield);
            }
            if (!errflag)
            {
                errflag = Efield_eval_e(
                    E, tempy[0], tempy[1], tempy[2], t0 + h[i] / 2.0, efield,
                    bfield);
            }
            if (!errflag)
            {
                step_gceom(
                    k2, tempy, mass, charge, b_db, E, alpha, Phi, aldforce);
            }
            for (int j = 0; j < 6; j++)
            {
                tempy[j] = yprev[j] + h[i] * k2[j] / 2.0;
            }

            if (!errflag)
            {
                errflag = Bfield_eval_b_db(
                    b_db, tempy[0], tempy[1], tempy[2], t0 + h[i] / 2.0,
                    bfield);
            }
            if (!errflag)
            {
                errflag = Efield_eval_e(
                    E, tempy[0], tempy[1], tempy[2], t0 + h[i] / 2.0, efield,
                    bfield);
            }
            if (!errflag)
            {
                step_gceom(
                    k3, tempy, mass, charge, b_db, E, alpha, Phi, aldforce);
            }
            for (int j = 0; j < 6; j++)
            {
                tempy[j] = yprev[j] + h[i] * k3[j];
            }

            if (!errflag)
            {
                errflag = Bfield_eval_b_db(
                    b_db, tempy[0], tempy[1], tempy[2], t0 + h[i], bfield);
            }
            if (!errflag)
            {
                errflag = Efield_eval_e(
                    E, tempy[0], tempy[1], tempy[2], t0 + h[i], efield, bfield);
            }
            if (!errflag)
            {
                step_gceom(
                    k4, tempy, mass, charge, b_db, E, alpha, Phi, aldforce);
            }
            for (int j = 0; j < 6; j++)
            {
                y[j] = yprev[j] +
                       h[i] / 6.0 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
            }

            /* Test that results are physical */
            errflag = ERROR_CHECK(
                errflag, y[0] <= 0, ERR_UNPHYSICAL_RESULT,
                SIMULATE_ORBIT_GC_RK4_C);
            errflag = ERROR_CHECK(
                errflag, y[4] < 0, ERR_UNPHYSICAL_RESULT,
                SIMULATE_ORBIT_GC_RK4_C);

            /* Update gc phase space position */
            if (!errflag)
            {
                mrk->r[i] = y[0];
                mrk->phi[i] = y[1];
                mrk->z[i] = y[2];
                mrk->ppar[i] = y[3];
                mrk->mu[i] = y[4];
                mrk->zeta[i] = fmod(y[5], CONST_2PI);
                if (mrk->zeta[i] < 0)
                {
                    mrk->zeta[i] = CONST_2PI + mrk->zeta[i];
                }
            }

            /* Evaluate magnetic field (and gradient) and rho at new position */
            real psi[1];
            real rho[2];
            if (!errflag)
            {
                errflag = Bfield_eval_b_db(
                    b_db, mrk->r[i], mrk->phi[i], mrk->z[i], t0 + h[i], bfield);
            }
            if (!errflag)
            {
                errflag = Bfield_eval_psi(
                    psi, mrk->r[i], mrk->phi[i], mrk->z[i], t0 + h[i], bfield);
            }
            if (!errflag)
            {
                errflag = Bfield_eval_rho(rho, psi[0], bfield);
            }

            if (!errflag)
            {
                mrk->br[i] = b_db[0];
                mrk->dbrdr[i] = b_db[3];
                mrk->dbrdphi[i] = b_db[4];
                mrk->dbrdz[i] = b_db[5];

                mrk->bphi[i] = b_db[1];
                mrk->dbphidr[i] = b_db[6];
                mrk->dbphidphi[i] = b_db[7];
                mrk->dbphidz[i] = b_db[8];

                mrk->bz[i] = b_db[2];
                mrk->dbzdr[i] = b_db[9];
                mrk->dbzdphi[i] = b_db[10];
                mrk->dbzdz[i] = b_db[11];

                mrk->rho[i] = rho[0];

                /* Evaluate theta angle so that it is cumulative */
                real axisrz[2];
                errflag = Bfield_eval_axis_rz(axisrz, bfield, mrk->phi[i]);
                mrk->theta[i] += atan2(
                    (R0 - axisrz[0]) * (mrk->z[i] - axisrz[1]) -
                        (z0 - axisrz[1]) * (mrk->r[i] - axisrz[0]),
                    (R0 - axisrz[0]) * (mrk->r[i] - axisrz[0]) +
                        (z0 - axisrz[1]) * (mrk->z[i] - axisrz[1]));
            }

            /* Error handling */
            if (errflag)
            {
                mrk->err[i] = errflag;
                mrk->running[i] = 0;
            }
        }
    }
}

void step_gc_rk4_mhd(
    MarkerGuidingCenter *mrk, const real *h, Bfield *bfield, Efield *efield,
    Boozer *boozer, Mhd *mhd, int aldforce)
{
    GPU_DATA_IS_MAPPED(h [0:mrk->size])
    GPU_PARALLEL_LOOP_ALL_LEVELS
    for (size_t i = 0; i < mrk->size; i++)
    {
        if (mrk->running[i])
        {
            err_t errflag = 0;

            real k1[6], k2[6], k3[6], k4[6];
            real tempy[6];
            real yprev[6];
            real y[6];

            real mass = mrk->mass;
            real charge = mrk->charge[i] * CONST_E;
            real b_db[15], E[3], alpha[5], Phi[5];

            real R0 = mrk->r[i];
            real z0 = mrk->z[i];
            real t0 = mrk->time[i];

            /* Coordinates are copied from the struct into an array to make
             * passing parameters easier */
            yprev[0] = mrk->r[i];
            yprev[1] = mrk->phi[i];
            yprev[2] = mrk->z[i];
            yprev[3] = mrk->ppar[i];
            yprev[4] = mrk->mu[i];
            yprev[5] = mrk->zeta[i];

            /* Magnetic field at initial position already known */
            b_db[0] = mrk->br[i];
            b_db[3] = mrk->dbrdr[i];
            b_db[4] = mrk->dbrdphi[i];
            b_db[5] = mrk->dbrdz[i];

            b_db[1] = mrk->bphi[i];
            b_db[6] = mrk->dbphidr[i];
            b_db[7] = mrk->dbphidphi[i];
            b_db[8] = mrk->dbphidz[i];

            b_db[2] = mrk->bz[i];
            b_db[9] = mrk->dbzdr[i];
            b_db[10] = mrk->dbzdphi[i];
            b_db[11] = mrk->dbzdz[i];

            if (!errflag)
            {
                errflag = Efield_eval_e(
                    E, yprev[0], yprev[1], yprev[2], t0, efield, bfield);
            }
            if (!errflag)
            {
                errflag = Mhd_eval_alpha_Phi(
                    alpha, Phi, yprev[0], yprev[1], yprev[2], t0,
                    MHD_INCLUDE_ALL, mhd, bfield, boozer);
            }
            if (!errflag)
            {
                step_gceom(
                    k1, yprev, mass, charge, b_db, E, alpha, Phi, aldforce);
            }

            /* particle coordinates for the subsequent ydot evaluations are
             * stored in tempy */
            for (int j = 0; j < 6; j++)
            {
                tempy[j] = yprev[j] + h[i] / 2.0 * k1[j];
            }

            if (!errflag)
            {
                errflag = Bfield_eval_b_db(
                    b_db, tempy[0], tempy[1], tempy[2], t0 + h[i] / 2.0,
                    bfield);
            }
            if (!errflag)
            {
                errflag = Efield_eval_e(
                    E, tempy[0], tempy[1], tempy[2], t0 + h[i] / 2.0, efield,
                    bfield);
            }
            if (!errflag)
            {
                errflag = Mhd_eval_alpha_Phi(
                    alpha, Phi, yprev[0], yprev[1], yprev[2], t0 + h[i] / 2.0,
                    MHD_INCLUDE_ALL, mhd, bfield, boozer);
            }
            if (!errflag)
            {
                step_gceom(
                    k2, tempy, mass, charge, b_db, E, alpha, Phi, aldforce);
            }
            for (int j = 0; j < 6; j++)
            {
                tempy[j] = yprev[j] + h[i] / 2.0 * k2[j];
            }

            if (!errflag)
            {
                errflag = Bfield_eval_b_db(
                    b_db, tempy[0], tempy[1], tempy[2], t0 + h[i] / 2.0,
                    bfield);
            }
            if (!errflag)
            {
                errflag = Efield_eval_e(
                    E, tempy[0], tempy[1], tempy[2], t0 + h[i] / 2.0, efield,
                    bfield);
            }
            if (!errflag)
            {
                errflag = Mhd_eval_alpha_Phi(
                    alpha, Phi, yprev[0], yprev[1], yprev[2], t0 + h[i] / 2.0,
                    MHD_INCLUDE_ALL, mhd, bfield, boozer);
            }
            if (!errflag)
            {
                step_gceom(
                    k3, tempy, mass, charge, b_db, E, alpha, Phi, aldforce);
            }
            for (int j = 0; j < 6; j++)
            {
                tempy[j] = yprev[j] + h[i] * k3[j];
            }

            if (!errflag)
            {
                errflag = Bfield_eval_b_db(
                    b_db, tempy[0], tempy[1], tempy[2], t0 + h[i], bfield);
            }
            if (!errflag)
            {
                errflag = Efield_eval_e(
                    E, tempy[0], tempy[1], tempy[2], t0 + h[i], efield, bfield);
            }
            if (!errflag)
            {
                errflag = Mhd_eval_alpha_Phi(
                    alpha, Phi, yprev[0], yprev[1], yprev[2], t0 + h[i] / 2.0,
                    MHD_INCLUDE_ALL, mhd, bfield, boozer);
            }
            if (!errflag)
            {
                step_gceom(
                    k4, tempy, mass, charge, b_db, E, alpha, Phi, aldforce);
            }
            for (int j = 0; j < 6; j++)
            {
                y[j] = yprev[j] +
                       h[i] / 6.0 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
            }

            errflag = ERROR_CHECK(
                errflag, y[0] <= 0, ERR_UNPHYSICAL_RESULT,
                SIMULATE_ORBIT_GC_RK4_C);
            errflag = ERROR_CHECK(
                errflag, y[4] < 0, ERR_UNPHYSICAL_RESULT,
                SIMULATE_ORBIT_GC_RK4_C);

            /* Update gc phase space position */
            if (!errflag)
            {
                mrk->r[i] = y[0];
                mrk->phi[i] = y[1];
                mrk->z[i] = y[2];
                mrk->ppar[i] = y[3];
                mrk->mu[i] = y[4];
                mrk->zeta[i] = fmod(y[5], CONST_2PI);
                if (mrk->zeta[i] < 0)
                {
                    mrk->zeta[i] = CONST_2PI + mrk->zeta[i];
                }
            }

            /* Evaluate magnetic field (and gradient) and rho at new position */
            real psi[1];
            real rho[2];
            if (!errflag)
            {
                errflag = Bfield_eval_b_db(
                    b_db, mrk->r[i], mrk->phi[i], mrk->z[i], t0 + h[i], bfield);
            }
            if (!errflag)
            {
                errflag = Bfield_eval_psi(
                    psi, mrk->r[i], mrk->phi[i], mrk->z[i], t0 + h[i], bfield);
            }
            if (!errflag)
            {
                errflag = Bfield_eval_rho(rho, psi[0], bfield);
            }

            if (!errflag)
            {
                mrk->br[i] = b_db[0];
                mrk->dbrdr[i] = b_db[3];
                mrk->dbrdphi[i] = b_db[4];
                mrk->dbrdz[i] = b_db[5];

                mrk->bphi[i] = b_db[1];
                mrk->dbphidr[i] = b_db[6];
                mrk->dbphidphi[i] = b_db[7];
                mrk->dbphidz[i] = b_db[8];

                mrk->bz[i] = b_db[2];
                mrk->dbzdr[i] = b_db[9];
                mrk->dbzdphi[i] = b_db[10];
                mrk->dbzdz[i] = b_db[11];
                mrk->rho[i] = rho[0];

                /* Evaluate pol angle so that it is cumulative */
                real axisrz[2];
                errflag = Bfield_eval_axis_rz(axisrz, bfield, mrk->phi[i]);
                mrk->theta[i] += atan2(
                    (R0 - axisrz[0]) * (mrk->z[i] - axisrz[1]) -
                        (z0 - axisrz[1]) * (mrk->r[i] - axisrz[0]),
                    (R0 - axisrz[0]) * (mrk->r[i] - axisrz[0]) +
                        (z0 - axisrz[1]) * (mrk->z[i] - axisrz[1]));
            }

            /* Error handling */
            if (errflag)
            {
                mrk->err[i] = errflag;
                mrk->running[i] = 0;
            }
        }
    }
}
