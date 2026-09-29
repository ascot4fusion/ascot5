/**
 * Gyro-orbit integrator implemented with VPA with fixed time-step
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
#include "parallel.h"
#include "utils/mathlib.h"
#include "utils/physlib.h"
#include <math.h>
#include <stdio.h>

void step_go_vpa(
    MarkerGyroOrbit *mrk, const real *h, Bfield *bfield, Efield *efield,
    int aldforce)
{
    GPU_DATA_IS_MAPPED(h [0:mrk->size])
    GPU_PARALLEL_LOOP_ALL_LEVELS
    for (size_t i = 0; i < mrk->size; i++)
    {
        if (mrk->running[i])
        {
            err_t errflag = 0;

            real R0 = mrk->r[i];
            real z0 = mrk->z[i];
            real t0 = mrk->time[i];
            real mass = mrk->mass;

            /* Convert velocity to cartesian coordinates */
            real prpz[3] = {mrk->p_r[i], mrk->p_phi[i], mrk->p_z[i]};
            real pxyz[3];
            math_vec_rpz2xyz(prpz, pxyz, mrk->phi[i]);

            real posrpz[3] = {mrk->r[i], mrk->phi[i], mrk->z[i]};
            real posxyz0[3], posxyz[3];
            math_rpz2xyz(posrpz, posxyz0);

            /* Take a half step and evaluate fields at that position */
            real gamma = physlib_gamma_pnorm(mass, math_norm(pxyz));
            posxyz[0] = posxyz0[0] + pxyz[0] * h[i] / (2.0 * gamma * mass);
            posxyz[1] = posxyz0[1] + pxyz[1] * h[i] / (2.0 * gamma * mass);
            posxyz[2] = posxyz0[2] + pxyz[2] * h[i] / (2.0 * gamma * mass);

            math_xyz2rpz(posxyz, posrpz);

            real Brpz[3];
            real Erpz[3];
            if (!errflag)
            {
                errflag = Bfield_eval_b(
                    Brpz, posrpz[0], posrpz[1], posrpz[2], t0 + h[i] / 2.0,
                    bfield);
            }
            if (!errflag)
            {
                errflag = Efield_eval_e(
                    Erpz, posrpz[0], posrpz[1], posrpz[2], t0 + h[i] / 2.0,
                    efield, bfield);
            }

            real fposxyz[3]; // final position in cartesian coordinates

            if (!errflag)
            {
                /* Electromagnetic fields to cartesian coordinates */
                real Bxyz[3];
                real Exyz[3];

                math_vec_rpz2xyz(Brpz, Bxyz, posrpz[1]);
                math_vec_rpz2xyz(Erpz, Exyz, posrpz[1]);

                /* Evaluate helper variable pminus */
                real pminus[3];
                real sigma =
                    mrk->charge[i] * CONST_E * h[i] / (2 * mrk->mass * CONST_C);
                pminus[0] = pxyz[0] / (mass * CONST_C) + sigma * Exyz[0];
                pminus[1] = pxyz[1] / (mass * CONST_C) + sigma * Exyz[1];
                pminus[2] = pxyz[2] / (mass * CONST_C) + sigma * Exyz[2];

                /* Second helper variable pplus*/
                real d = (mrk->charge[i] * CONST_E * h[i] / (2 * mrk->mass)) /
                         sqrt(1 + math_dot(pminus, pminus));
                real d2 = d * d;

                real Bhat[9] = {0,       Bxyz[2], -Bxyz[1], -Bxyz[2], 0,
                                Bxyz[0], Bxyz[1], -Bxyz[0], 0};
                real Bhat2[9];
                math_matmul(Bhat, Bhat, Bhat2);

                real B2 =
                    Bxyz[0] * Bxyz[0] + Bxyz[1] * Bxyz[1] + Bxyz[2] * Bxyz[2];

                real A[9];
                for (int j = 0; j < 9; j++)
                {
                    A[j] =
                        (Bhat[j] + d * Bhat2[j]) * (2.0 * d / (1.0 + d2 * B2));
                }

                real pplus[3];
                math_matvecmul(pminus, A, pplus);

                /* Take the step */
                real pfinal[3];
                pfinal[0] = pminus[0] + pplus[0] + sigma * Exyz[0];
                pfinal[1] = pminus[1] + pplus[1] + sigma * Exyz[1];
                pfinal[2] = pminus[2] + pplus[2] + sigma * Exyz[2];

                pxyz[0] = pfinal[0] * mass * CONST_C;
                pxyz[1] = pfinal[1] * mass * CONST_C;
                pxyz[2] = pfinal[2] * mass * CONST_C;
            }

            gamma = physlib_gamma_pnorm(mass, math_norm(pxyz));
            fposxyz[0] = posxyz[0] + h[i] * pxyz[0] / (2.0 * gamma * mass);
            fposxyz[1] = posxyz[1] + h[i] * pxyz[1] / (2.0 * gamma * mass);
            fposxyz[2] = posxyz[2] + h[i] * pxyz[2] / (2.0 * gamma * mass);

            if (!errflag)
            {
                /* Back to cylindrical coordinates */
                mrk->r[i] =
                    sqrt(fposxyz[0] * fposxyz[0] + fposxyz[1] * fposxyz[1]);

                /* phi is evaluated like this to make sure it is cumulative */
                mrk->phi[i] += atan2(
                    posxyz0[0] * fposxyz[1] - posxyz0[1] * fposxyz[0],
                    posxyz0[0] * fposxyz[0] + posxyz0[1] * fposxyz[1]);
                mrk->z[i] = fposxyz[2];

                real cosp = cos(mrk->phi[i]);
                real sinp = sin(mrk->phi[i]);
                mrk->p_r[i] = pxyz[0] * cosp + pxyz[1] * sinp;
                mrk->p_phi[i] = -pxyz[0] * sinp + pxyz[1] * cosp;
                mrk->p_z[i] = pxyz[2];
            }

            /* Evaluate magnetic field (and gradient) and rho at new position */
            real b_db[15];
            real psi[1];
            real rho[2];
            if (!errflag)
            {
                errflag = Bfield_eval_b_db(
                    b_db, mrk->r[i], mrk->phi[i], mrk->z[i], t0 + h[i],
                    bfield);
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

                /* Evaluate phi and theta angles so that they are cumulative */
                real axisrz[2];
                errflag = Bfield_eval_axis_rz(axisrz, bfield, mrk->phi[i]);
                mrk->theta[i] += atan2(
                    (R0 - axisrz[0]) * (mrk->z[i] - axisrz[1]) -
                        (z0 - axisrz[1]) * (mrk->r[i] - axisrz[0]),
                    (R0 - axisrz[0]) * (mrk->r[i] - axisrz[0]) +
                        (z0 - axisrz[1]) * (mrk->z[i] - axisrz[1]));
            }

            /* Evaluate Abraham-Lorentz-Dirac force (if enabled) is evaluated
             * separately using the Euler method */
            real Bnorm = math_normc(mrk->br[i], mrk->bphi[i], mrk->bz[i]);
            real pnorm = math_normc(mrk->p_r[i], mrk->p_phi[i], mrk->p_z[i]);
            real t_ald = phys_ald_force_chartime(
                             mrk->charge[i] * CONST_E, mrk->mass, Bnorm, gamma) *
                         aldforce;
            real pparbhatperB =
                (mrk->p_r[i] * mrk->br[i] + mrk->p_phi[i] * mrk->bphi[i] +
                 mrk->p_z[i] * mrk->bz[i]) /
                (Bnorm * Bnorm * pnorm);
            real pperpvec[3] = {
                mrk->p_r[i] - pparbhatperB * mrk->br[i],
                mrk->p_phi[i] - pparbhatperB * mrk->bphi[i],
                mrk->p_z[i] - pparbhatperB * mrk->bz[i]};
            real C = (pperpvec[0] * pperpvec[0] + pperpvec[1] * pperpvec[1] +
                      pperpvec[2] * pperpvec[2]) /
                     (mrk->mass * mrk->mass * CONST_C2);
            mrk->p_r[i] -= t_ald * (pperpvec[0] + C * mrk->p_r[i]);
            mrk->p_phi[i] -= t_ald * (pperpvec[1] + C * mrk->p_phi[i]);
            mrk->p_z[i] -= t_ald * (pperpvec[2] + C * mrk->p_z[i]);

            /* Error handling */
            if (errflag)
            {
                mrk->err[i] = errflag;
                mrk->running[i] = 0;
            }
        }
    }
}

void step_go_vpa_mhd(
    MarkerGyroOrbit *mrk, const real *h, Bfield *bfield, Efield *efield,
    Boozer *boozer, Mhd *mhd, int aldforce)
{
    (void)aldforce; // TODO
    GPU_DATA_IS_MAPPED(h [0:mrk->size])
    GPU_PARALLEL_LOOP_ALL_LEVELS
    for (size_t i = 0; i < mrk->size; i++)
    {
        if (mrk->running[i])
        {
            err_t errflag = 0;

            real R0 = mrk->r[i];
            real z0 = mrk->z[i];
            real t0 = mrk->time[i];
            real mass = mrk->mass;

            /* Convert velocity to cartesian coordinates */
            real prpz[3] = {mrk->p_r[i], mrk->p_phi[i], mrk->p_z[i]};
            real pxyz[3];
            math_vec_rpz2xyz(prpz, pxyz, mrk->phi[i]);

            real posrpz[3] = {mrk->r[i], mrk->phi[i], mrk->z[i]};
            real posxyz0[3], posxyz[3];
            math_rpz2xyz(posrpz, posxyz0);

            /* Take a half step and evaluate fields at that position */
            real gamma = physlib_gamma_pnorm(mass, math_norm(pxyz));
            posxyz[0] = posxyz0[0] + pxyz[0] * h[i] / (2 * gamma * mass);
            posxyz[1] = posxyz0[1] + pxyz[1] * h[i] / (2 * gamma * mass);
            posxyz[2] = posxyz0[2] + pxyz[2] * h[i] / (2 * gamma * mass);

            math_xyz2rpz(posxyz, posrpz);

            real Brpz[3], Erpz[3], Epert[3], Psi[1];
            if (!errflag)
            {
                errflag = Efield_eval_e(
                    Erpz, posrpz[0], posrpz[1], posrpz[2], t0 + h[i] / 2,
                    efield, bfield);
            }

            int pertonly = 0;
            if (!errflag)
            {
                errflag = Mhd_eval_perturbation(
                    Brpz, Epert, Psi, posrpz[0], posrpz[1], posrpz[2],
                    t0 + h[i] / 2, pertonly, MHD_INCLUDE_ALL, mhd, bfield,
                    boozer);
            }
            Erpz[0] += Epert[0];
            Erpz[1] += Epert[1];
            Erpz[2] += Epert[2];

            real fposxyz[3]; // final position in cartesian coordinates

            if (!errflag)
            {
                /* Electromagnetic fields to cartesian coordinates */
                real Bxyz[3];
                real Exyz[3];

                math_vec_rpz2xyz(Brpz, Bxyz, posrpz[1]);
                math_vec_rpz2xyz(Erpz, Exyz, posrpz[1]);

                /* Evaluate helper variable pminus */
                real pminus[3];
                real sigma =
                    mrk->charge[i] * CONST_E * h[i] / (2 * mrk->mass * CONST_C);
                pminus[0] = pxyz[0] / (mass * CONST_C) + sigma * Exyz[0];
                pminus[1] = pxyz[1] / (mass * CONST_C) + sigma * Exyz[1];
                pminus[2] = pxyz[2] / (mass * CONST_C) + sigma * Exyz[2];

                /* Second helper variable pplus*/
                real d = (mrk->charge[i] * CONST_E * h[i] / (2 * mrk->mass)) /
                         sqrt(1 + math_dot(pminus, pminus));
                real d2 = d * d;

                real Bhat[9] = {0,       Bxyz[2], -Bxyz[1], -Bxyz[2], 0,
                                Bxyz[0], Bxyz[1], -Bxyz[0], 0};
                real Bhat2[9];
                math_matmul(Bhat, Bhat, Bhat2);

                real B2 =
                    Bxyz[0] * Bxyz[0] + Bxyz[1] * Bxyz[1] + Bxyz[2] * Bxyz[2];

                real A[9];
                for (int j = 0; j < 9; j++)
                {
                    A[j] = (Bhat[j] + d * Bhat2[j]) * (2.0 * d / (1 + d2 * B2));
                }

                real pplus[3];
                math_matvecmul(pminus, A, pplus);

                /* Take the step */
                real pfinal[3];
                pfinal[0] = pminus[0] + pplus[0] + sigma * Exyz[0];
                pfinal[1] = pminus[1] + pplus[1] + sigma * Exyz[1];
                pfinal[2] = pminus[2] + pplus[2] + sigma * Exyz[2];

                pxyz[0] = pfinal[0] * mass * CONST_C;
                pxyz[1] = pfinal[1] * mass * CONST_C;
                pxyz[2] = pfinal[2] * mass * CONST_C;
            }

            gamma = physlib_gamma_pnorm(mass, math_norm(pxyz));
            fposxyz[0] = posxyz[0] + h[i] * pxyz[0] / (2 * gamma * mass);
            fposxyz[1] = posxyz[1] + h[i] * pxyz[1] / (2 * gamma * mass);
            fposxyz[2] = posxyz[2] + h[i] * pxyz[2] / (2 * gamma * mass);

            if (!errflag)
            {
                /* Back to cylindrical coordinates */
                mrk->r[i] =
                    sqrt(fposxyz[0] * fposxyz[0] + fposxyz[1] * fposxyz[1]);

                /* phi is evaluated like this to make sure it is cumulative */
                mrk->phi[i] += atan2(
                    posxyz0[0] * fposxyz[1] - posxyz0[1] * fposxyz[0],
                    posxyz0[0] * fposxyz[0] + posxyz0[1] * fposxyz[1]);
                mrk->z[i] = fposxyz[2];

                real cosp = cos(mrk->phi[i]);
                real sinp = sin(mrk->phi[i]);
                mrk->p_r[i] = pxyz[0] * cosp + pxyz[1] * sinp;
                mrk->p_phi[i] = -pxyz[0] * sinp + pxyz[1] * cosp;
                mrk->p_z[i] = pxyz[2];
            }

            /* Evaluate magnetic field (and gradient) and rho at new position */
            real b_db[15];
            real psi[1];
            real rho[2];
            if (!errflag)
            {
                errflag = Bfield_eval_b_db(
                    b_db, mrk->r[i], mrk->phi[i], mrk->z[i], t0 + h[i],
                    bfield);
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

                /* Evaluate phi and theta angles so that they are cumulative */
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
