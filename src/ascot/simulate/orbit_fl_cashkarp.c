/**
 * Field line integrator implemented with adaptive Cash-Karp method
 * (see orbit_following.h).
 **/
#include "data/bfield.h"
#include "data/boozer.h"
#include "data/marker.h"
#include "data/mhd.h"
#include "defines.h"
#include "orbit_following.h"
#include "utils/mathlib.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>

void step_fl_cashkarp(
    MarkerFieldLine *mrk, const real *h, real *hnext, real tol, Bfield *bfield)
{
    GPU_DATA_IS_MAPPED(h [0:mrk->size], hnext [0:mrk->size])
    GPU_PARALLEL_LOOP_ALL_LEVELS
    for (size_t i = 0; i < mrk->size; i++)
    {
        if (mrk->running[i])
        {
            err_t errflag = 0;

            real k2[3], k3[3], k4[3], k5[3], k6[3];
            real tempy[3];

            real normB;

            real R0 = mrk->r[i];
            real z0 = mrk->z[i];
            real t0 = mrk->time[i];
            int direction = 1 - 2 * (mrk->pitch[i] < 0);

            real yprev[3] = {mrk->r[i], mrk->phi[i], mrk->z[i]};
            real k1[3] = {mrk->br[i], mrk->bphi[i], mrk->bz[i]};
            normB = (math_normc(k1[0], k1[1], k1[2])) * direction;
            k1[0] /= normB;
            k1[1] /= normB * yprev[0];
            k1[2] /= normB;

            for (int j = 0; j < 3; j++)
            {
                tempy[j] = yprev[j] + h[i] * ((1.0 / 5) * k1[j]);
            }
            if (!errflag)
            {
                errflag = Bfield_eval_b(
                    k2, tempy[0], tempy[1], tempy[2], t0 + (1.0 / 5) * h[i],
                    bfield);
            }
            normB = (math_normc(k2[0], k2[1], k2[2])) * direction;
            k2[0] /= normB;
            k2[1] /= normB;
            k2[2] /= normB;
            k2[1] /= tempy[0];

            for (int j = 0; j < 3; j++)
            {
                tempy[j] =
                    yprev[j] + h[i] * ((3.0 / 40) * k1[j] + (9.0 / 40) * k2[j]);
            }
            if (!errflag)
            {
                errflag = Bfield_eval_b(
                    k3, tempy[0], tempy[1], tempy[2], t0 + (3.0 / 10) * h[i],
                    bfield);
            }
            normB = (math_normc(k3[0], k3[1], k3[2])) * direction;
            k3[0] /= normB;
            k3[1] /= normB;
            k3[2] /= normB;
            k3[1] /= tempy[0];

            for (int j = 0; j < 3; j++)
            {
                tempy[j] =
                    yprev[j] + h[i] * ((3.0 / 10) * k1[j] +
                                       (-9.0 / 10) * k2[j] + (6.0 / 5) * k3[j]);
            }
            if (!errflag)
            {
                errflag = Bfield_eval_b(
                    k4, tempy[0], tempy[1], tempy[2], t0 + (3.0 / 5) * h[i],
                    bfield);
            }
            normB = (math_normc(k4[0], k4[1], k4[2])) * direction;
            k4[0] /= normB;
            k4[1] /= normB;
            k4[2] /= normB;
            k4[1] /= tempy[0];

            for (int j = 0; j < 3; j++)
            {
                tempy[j] = yprev[j] +
                           h[i] * ((-11.0 / 54) * k1[j] + (5.0 / 2) * k2[j] +
                                   (-70.0 / 27) * k3[j] + (35.0 / 27) * k4[j]);
            }
            if (!errflag)
            {
                errflag = Bfield_eval_b(
                    k5, tempy[0], tempy[1], tempy[2], t0 + h[i], bfield);
            }
            normB = (math_normc(k5[0], k5[1], k5[2])) * direction;
            k5[0] /= normB;
            k5[1] /= normB;
            k5[2] /= normB;
            k5[1] /= tempy[0];

            for (int j = 0; j < 3; j++)
            {
                tempy[j] = yprev[j] + h[i] * ((1631.0 / 55296) * k1[j] +
                                              (175.0 / 512) * k2[j] +
                                              (575.0 / 13824) * k3[j] +
                                              (44275.0 / 110592) * k4[j] +
                                              (253.0 / 4096) * k5[j]);
            }
            if (!errflag)
            {
                errflag = Bfield_eval_b(
                    k6, tempy[0], tempy[1], tempy[2], t0 + (7.0 / 8) * h[i],
                    bfield);
            }
            normB = (math_normc(k6[0], k6[1], k6[2])) * direction;
            k6[0] /= normB;
            k6[1] /= normB;
            k6[2] /= normB;
            k6[1] /= tempy[0];

            real rk5[3];
            real err = 0.0;
            for (int j = 0; j < 3; j++)
            {
                rk5[j] =
                    yprev[j] +
                    h[i] * ((37.0 / 378) * k1[j] + (250.0 / 621) * k3[j] +
                            (125.0 / 594) * k4[j] + (512.0 / 1771) * k6[j]);

                real rk4 = yprev[j] +
                           h[i] * ((2825.0 / 27648) * k1[j] +
                                   (18575.0 / 48384) * k3[j] +
                                   (13525.0 / 55296) * k4[j] +
                                   (277.0 / 14336) * k5[j] + (1.0 / 4) * k6[j]);

                real yerr = fabs(rk5[j] - rk4);
                real ytol = fabs(yprev[j]) + fabs(k1[j] * h[i]) + DBL_EPSILON;
                err = fmax(err, yerr / ytol);
            }

            err = err / tol;
            if (err <= 1)
            {
                /* Time step accepted */
                hnext[i] = 0.85 * h[i] * pow(err, -0.2);
            }
            else
            {
                /* Time step rejected */
                hnext[i] = -0.85 * h[i] * pow(err, -0.25);
            }

            mrk->r[i] = rk5[0];
            mrk->phi[i] = rk5[1];
            mrk->z[i] = rk5[2];

            /* Evaluate magnetic field (and gradient) and rho at new position */
            real b_db[15];
            if (!errflag)
            {
                errflag = Bfield_eval_b_db(
                    b_db, mrk->r[i], mrk->phi[i], mrk->z[i],
                    mrk->time[i] + h[i], bfield);
            }
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

            real psi[1];
            real rho[2];
            if (!errflag)
            {
                errflag = Bfield_eval_psi(
                    psi, mrk->r[i], mrk->phi[i], mrk->z[i], mrk->time[i] + h[i],
                    bfield);
            }
            if (!errflag)
            {
                errflag = Bfield_eval_rho(rho, psi[0], bfield);
            }
            mrk->rho[i] = rho[0];

            /* Evaluate theta angle so that it is cumulative */
            if (!errflag)
            {
                real axisrz[2];
                errflag = Bfield_eval_axis_rz(axisrz, bfield, mrk->phi[i]);
                mrk->theta[i] += atan2(
                    (R0 - axisrz[0]) * (mrk->z[i] - axisrz[1]) -
                        (z0 - axisrz[1]) * (mrk->r[i] - axisrz[0]),
                    (R0 - axisrz[0]) * (mrk->r[i] - axisrz[0]) +
                        (z0 - axisrz[1]) * (mrk->z[i] - axisrz[1]));
            }

            if (errflag)
            {
                mrk->err[i] = errflag;
            }
        }
    }
}

void step_fl_cashkarp_mhd(
    MarkerFieldLine *mrk, const real *h, real *hnext, real tol, Bfield *bfield,
    Boozer *boozerdata, Mhd *mhddata)
{
    GPU_DATA_IS_MAPPED(h [0:mrk->size], h_next [0:mrk->size])
    GPU_PARALLEL_LOOP_ALL_LEVELS
    for (size_t i = 0; i < mrk->size; i++)
    {
        if (mrk->running[i])
        {
            err_t errflag = 0;

            real k1[3], k2[3], k3[3], k4[3], k5[3], k6[3];
            real tempy[3];
            real yprev[3];

            real normB;
            real bpert[3], epert[3], Phipert[1];

            real R0 = mrk->r[i];
            real z0 = mrk->z[i];
            real t0 = mrk->time[i];

            /* Direction */
            int direction = 1;
            if (mrk->pitch[i] < 0)
            {
                direction = -1;
            }

            int pertonly = 0;

            /* Coordinates are copied from the struct into an array to make
             * passing parameters easier */
            yprev[0] = mrk->r[i];
            yprev[1] = mrk->phi[i];
            yprev[2] = mrk->z[i];

            if (!errflag)
            {
                errflag = Mhd_eval_perturbation(
                    bpert, epert, Phipert, tempy[0], tempy[1], tempy[2], t0,
                    pertonly, MHD_INCLUDE_ALL, mhddata, bfield, boozerdata);
            }
            k1[0] = bpert[0];
            k1[1] = bpert[1];
            k1[2] = bpert[2];

            normB = (math_normc(k1[0], k1[1], k1[2])) * direction;
            k1[0] /= normB;
            k1[1] /= normB;
            k1[2] /= normB;
            k1[1] /= yprev[0];

            for (int j = 0; j < 3; j++)
            {
                tempy[j] = yprev[j] + h[i] * ((1.0 / 5) * k1[j]);
            }
            if (!errflag)
            {
                errflag = Mhd_eval_perturbation(
                    bpert, epert, Phipert, tempy[0], tempy[1], tempy[2], t0,
                    pertonly, MHD_INCLUDE_ALL, mhddata, bfield, boozerdata);
            }
            k2[0] = bpert[0];
            k2[1] = bpert[1];
            k2[2] = bpert[2];

            normB = (math_normc(k2[0], k2[1], k2[2])) * direction;
            k2[0] /= normB;
            k2[1] /= normB;
            k2[2] /= normB;
            k2[1] /= tempy[0];

            for (int j = 0; j < 3; j++)
            {
                tempy[j] =
                    yprev[j] + h[i] * ((3.0 / 40) * k1[j] + (9.0 / 40) * k2[j]);
            }
            if (!errflag)
            {
                errflag = Mhd_eval_perturbation(
                    bpert, epert, Phipert, tempy[0], tempy[1], tempy[2], t0,
                    pertonly, MHD_INCLUDE_ALL, mhddata, bfield, boozerdata);
            }
            k3[0] = bpert[0];
            k3[1] = bpert[1];
            k3[2] = bpert[2];

            normB = (math_normc(k3[0], k3[1], k3[2])) * direction;
            k3[0] /= normB;
            k3[1] /= normB;
            k3[2] /= normB;
            k3[1] /= tempy[0];

            for (int j = 0; j < 3; j++)
            {
                tempy[j] =
                    yprev[j] + h[i] * ((3.0 / 10) * k1[j] +
                                       (-9.0 / 10) * k2[j] + (6.0 / 5) * k3[j]);
            }
            if (!errflag)
            {
                errflag = Mhd_eval_perturbation(
                    bpert, epert, Phipert, tempy[0], tempy[1], tempy[2], t0,
                    pertonly, MHD_INCLUDE_ALL, mhddata, bfield, boozerdata);
            }
            k4[0] = bpert[0];
            k4[1] = bpert[1];
            k4[2] = bpert[2];

            normB = (math_normc(k4[0], k4[1], k4[2])) * direction;
            k4[0] /= normB;
            k4[1] /= normB;
            k4[2] /= normB;
            k4[1] /= tempy[0];

            for (int j = 0; j < 3; j++)
            {
                tempy[j] = yprev[j] +
                           h[i] * ((-11.0 / 54) * k1[j] + (5.0 / 2) * k2[j] +
                                   (-70.0 / 27) * k3[j] + (35.0 / 27) * k4[j]);
            }
            if (!errflag)
            {
                errflag = Mhd_eval_perturbation(
                    bpert, epert, Phipert, tempy[0], tempy[1], tempy[2], t0,
                    pertonly, MHD_INCLUDE_ALL, mhddata, bfield, boozerdata);
            }
            k5[0] = bpert[0];
            k5[1] = bpert[1];
            k5[2] = bpert[2];

            normB = (math_normc(k5[0], k5[1], k5[2])) * direction;
            k5[0] /= normB;
            k5[1] /= normB;
            k5[2] /= normB;
            k5[1] /= tempy[0];

            for (int j = 0; j < 3; j++)
            {
                tempy[j] = yprev[j] + h[i] * ((1631.0 / 55296) * k1[j] +
                                              (175.0 / 512) * k2[j] +
                                              (575.0 / 13824) * k3[j] +
                                              (44275.0 / 110592) * k4[j] +
                                              (253.0 / 4096) * k5[j]);
            }
            if (!errflag)
            {
                errflag = Mhd_eval_perturbation(
                    bpert, epert, Phipert, tempy[0], tempy[1], tempy[2], t0,
                    pertonly, MHD_INCLUDE_ALL, mhddata, bfield, boozerdata);
            }
            k6[0] = bpert[0];
            k6[1] = bpert[1];
            k6[2] = bpert[2];

            normB = (math_normc(k6[0], k6[1], k6[2])) * direction;
            k6[0] /= normB;
            k6[1] /= normB;
            k6[2] /= normB;
            k6[1] /= tempy[0];

            /* Error estimate is a difference between RK4 and RK5 solutions. If
             * time-step is accepted, the RK5 solution will be used to advance
             * marker. */
            real rk5[3];
            real err = 0.0;
            for (int j = 0; j < 3; j++)
            {
                rk5[j] =
                    yprev[j] +
                    h[i] * ((37.0 / 378) * k1[j] + (250.0 / 621) * k3[j] +
                            (125.0 / 594) * k4[j] + (512.0 / 1771) * k6[j]);

                real rk4 = yprev[j] +
                           h[i] * ((2825.0 / 27648) * k1[j] +
                                   (18575.0 / 48384) * k3[j] +
                                   (13525.0 / 55296) * k4[j] +
                                   (277.0 / 14336) * k5[j] + (1.0 / 4) * k6[j]);

                real yerr = fabs(rk5[j] - rk4);
                real ytol = fabs(yprev[j]) + fabs(k1[j] * h[i]) + DBL_EPSILON;
                err = fmax(err, yerr / ytol);
            }

            err = err / tol;
            if (err <= 1)
            {
                /* Time step accepted */
                hnext[i] = 0.85 * h[i] * pow(err, -0.2);
            }
            else
            {
                /* Time step rejected */
                hnext[i] = -0.85 * h[i] * pow(err, -0.25);
            }

            mrk->r[i] = rk5[0];
            mrk->phi[i] = rk5[1];
            mrk->z[i] = rk5[2];

            /* Evaluate magnetic field (and gradient) and rho at new position */
            real b_db[15];
            if (!errflag)
            {
                errflag = Bfield_eval_b_db(
                    b_db, mrk->r[i], mrk->phi[i], mrk->z[i],
                    mrk->time[i] + h[i], bfield);
            }
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

            real psi[1];
            real rho[2];
            if (!errflag)
            {
                errflag = Bfield_eval_psi(
                    psi, mrk->r[i], mrk->phi[i], mrk->z[i], mrk->time[i] + h[i],
                    bfield);
            }
            if (!errflag)
            {
                errflag = Bfield_eval_rho(rho, psi[0], bfield);
            }
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
    }
}
