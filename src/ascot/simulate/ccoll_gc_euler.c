/**
 * @file mccc_gc_euler.c
 * @brief Euler-Maruyama integrator for collision operator in GC picture.
 */
#include "consts.h"
#include "coulomb_collisions.h"
#include "data/bfield.h"
#include "data/marker.h"
#include "data/plasma.h"
#include "defines.h"
#include "utils/mathlib.h"
#include "utils/physlib.h"
#include "utils/random.h"
#include <math.h>

void mccc_gc_euler(
    MarkerGuidingCenter *p, real *h, Bfield *bfield, Plasma *plasma,
    mccc_data *mdata, real *rnd)
{

    /* Get plasma information before going to the  SIMD loop */
    int n_species = Plasma_get_n_species(plasma);
    const real *qb = Plasma_get_species_charge(plasma);
    const real *mb = Plasma_get_species_mass(plasma);

#pragma omp simd
    for (int i = 0; i < NSIMD; i++)
    {
        if (p->running[i])
        {
            err_t errflag = 0;

            /* Initial (R,z) position and magnetic field are needed for later */
            real Brpz[3] = {p->B_r[i], p->B_phi[i], p->B_z[i]};
            real Bnorm = math_norm(Brpz);
            real Bxyz[3];
            math_vec_rpz2xyz(Brpz, Bxyz, p->phi[i]);
            real R0 = p->r[i];
            real z0 = p->z[i];

            /* Move guiding center to (x, y, z, vnorm, xi) coordinates */
            real vin, pin, vflow, gamma, ppar_flow, xiin, Xin_xyz[3];
            Xin_xyz[0] = p->r[i] * cos(p->phi[i]);
            Xin_xyz[1] = p->r[i] * sin(p->phi[i]);
            Xin_xyz[2] = p->z[i];
            if (!errflag)
            {
                errflag = Plasma_eval_flow(
                    &vflow, p->rho[i], p->r[i], p->phi[i], p->z[i], p->time[i],
                    plasma);
            }
            gamma = physlib_gamma_ppar(p->mass[i], p->mu[i], p->ppar[i], Bnorm);
            ppar_flow = p->ppar[i] - gamma * vflow * p->mass[i];
            pin = physlib_gc_p(p->mass[i], p->mu[i], ppar_flow, Bnorm);
            xiin = physlib_gc_xi(p->mass[i], p->mu[i], ppar_flow, Bnorm);
            vin = physlib_vnorm_pnorm(p->mass[i], pin);

            /* Evaluate plasma density and temperature */
            real nb[MAX_SPECIES], Tb[MAX_SPECIES];
            if (!errflag)
            {
                errflag = Plasma_eval_nT(
                    nb, Tb, p->rho[i], p->r[i], p->phi[i], p->z[i], p->time[i],
                    plasma);
            }

            /* Coulomb logarithm */
            real clogab[MAX_SPECIES];
            mccc_coefs_clog(
                clogab, p->mass[i], p->charge[i], vin, n_species, mb, qb, nb,
                Tb);

            /* Evaluate collision coefficients and sum them for each *
             * species                                               */
            real gyrofreq =
                phys_gyrofreq_pnorm(p->mass[i], p->charge[i], pin, Bnorm);
            real K = 0, Dpara = 0, nu = 0, DX = 0;
            for (int j = 0; j < n_species; j++)
            {
                real vb = sqrt(2 * Tb[j] / mb[j]);
                real x = vin / vb;
                real mufun[3];
                mccc_coefs_mufun(mufun, x, mdata); // eq. 2.83 PhD Hirvijoki

                real Qb = mccc_coefs_Q(
                    p->mass[i], p->charge[i], mb[j], qb[j], nb[j], vb,
                    clogab[j], mufun[0]);
                real Dparab = mccc_coefs_Dpara(
                    p->mass[i], p->charge[i], vin, qb[j], nb[j], vb, clogab[j],
                    mufun[0]);
                real Dperpb = mccc_coefs_Dperp(
                    p->mass[i], p->charge[i], vin, qb[j], nb[j], vb, clogab[j],
                    mufun[1]);
                real dDparab = mccc_coefs_dDpara(
                    p->mass[i], p->charge[i], vin, qb[j], nb[j], vb, clogab[j],
                    mufun[0], mufun[2]);

                K += mccc_coefs_K(vin, Dparab, dDparab, Qb);
                Dpara += Dparab;
                nu += mccc_coefs_nu(vin, Dperpb); // eq.41
                DX += mccc_coefs_DX(xiin, Dparab, Dperpb, gyrofreq);
            }

            /* Evaluate collisions */
            real sdt = sqrt(h[i]);
            real dW[5];
            dW[0] = sdt * rnd[0 * NSIMD + i]; // For X_1
            dW[1] = sdt * rnd[1 * NSIMD + i]; // For X_2
            dW[2] = sdt * rnd[2 * NSIMD + i]; // For X_3
            dW[3] = sdt * rnd[3 * NSIMD + i]; // For v
            dW[4] = sdt * rnd[4 * NSIMD + i]; // For xi

            real bhat[3];
            math_unit(Bxyz, bhat);

            real k1 = sqrt(2 * DX);
            real k2 = math_dot(bhat, dW);

            real vout, xiout, Xout_xyz[3];
            Xout_xyz[0] = Xin_xyz[0] + k1 * (dW[0] - k2 * bhat[0]);
            Xout_xyz[1] = Xin_xyz[1] + k1 * (dW[1] - k2 * bhat[1]);
            Xout_xyz[2] = Xin_xyz[2] + k1 * (dW[2] - k2 * bhat[2]);
            vout = vin + K * h[i] + sqrt(2 * Dpara) * dW[3];
            xiout =
                xiin - xiin * nu * h[i] + sqrt((1 - xiin * xiin) * nu) * dW[4];

            /* Enforce boundary conditions */
            real cutoff = MCCC_CUTOFF * sqrt(Tb[0] / p->mass[i]);
            if (vout < cutoff)
            {
                vout = 2 * cutoff - vout;
            }

            if (fabs(xiout) > 1)
            {
                xiout = ((xiout > 0) - (xiout < 0)) * (2 - fabs(xiout));
            }

            /* Remove energy or pitch change or spatial diffusion from the    *
             * results if that is requested                                   */
            if (!mdata->include_energy)
            {
                vout = vin;
            }
            if (!mdata->include_pitch)
            {
                xiout = xiin;
            }
            if (!mdata->include_gcdiff)
            {
                Xout_xyz[0] = Xin_xyz[0];
                Xout_xyz[1] = Xin_xyz[1];
                Xout_xyz[2] = Xin_xyz[2];
            }
            real pout = physlib_pnorm_vnorm(p->mass[i], vout);

            /* Back to cylindrical coordinates */
            real Xout_rpz[3];
            math_xyz2rpz(Xout_xyz, Xout_rpz);

            /* Evaluate magnetic field (and gradient) and rho at new position */
            real B_dB[15], psi[1], rho[2];
            if (!errflag)
            {
                errflag = Bfield_eval_b_db(
                    B_dB, Xout_rpz[0], Xout_rpz[1], Xout_rpz[2],
                    p->time[i] + h[i], bfield);
            }
            if (!errflag)
            {
                errflag = Bfield_eval_psi(
                    psi, Xout_rpz[0], Xout_rpz[1], Xout_rpz[2],
                    p->time[i] + h[i], bfield);
            }
            if (!errflag)
            {
                errflag = Bfield_eval_rho(rho, psi[0], bfield);
            }

            if (!errflag)
            {
                /* Update marker coordinates at the new position */
                p->B_r[i] = B_dB[0];
                p->B_r_dr[i] = B_dB[1];
                p->B_r_dphi[i] = B_dB[2];
                p->B_r_dz[i] = B_dB[3];

                p->B_phi[i] = B_dB[4];
                p->B_phi_dr[i] = B_dB[5];
                p->B_phi_dphi[i] = B_dB[6];
                p->B_phi_dz[i] = B_dB[7];

                p->B_z[i] = B_dB[8];
                p->B_z_dr[i] = B_dB[9];
                p->B_z_dphi[i] = B_dB[10];
                p->B_z_dz[i] = B_dB[11];

                p->rho[i] = rho[0];

                Bnorm = math_normc(B_dB[0], B_dB[4], B_dB[8]);

                p->r[i] = Xout_rpz[0];
                p->z[i] = Xout_rpz[2];

                /*Since we use xiout (in "flow frame") here, we will get the mu
                in the "flow frame". If we assume that the flow frame has the
                same perpendicular velocity, mu remains unchanged under the
                co-ordinate transformation and thus this mu is then the same as
                the mu in the lab frame.*/
                p->mu[i] = physlib_gc_mu(p->mass[i], pout, xiout, Bnorm);
                gamma = physlib_gamma_pnorm(p->mass[i], pout);

                /* difference in ppar between the lab frame and "flow frame"  */
                real dppar = gamma * p->mass[i] * vflow;

                /* p->ppar is in the lab frame; xiout and pout are in the
                "flow frame" */
                p->ppar[i] = physlib_gc_ppar(pout, xiout) + dppar;

                /* Evaluate phi and theta angles so that they are cumulative */
                real axisrz[2];
                errflag = Bfield_eval_axis_rz(axisrz, bfield, p->phi[i]);
                p->theta[i] += atan2(
                    (R0 - axisrz[0]) * (p->z[i] - axisrz[1]) -
                        (z0 - axisrz[1]) * (p->r[i] - axisrz[0]),
                    (R0 - axisrz[0]) * (p->r[i] - axisrz[0]) +
                        (z0 - axisrz[1]) * (p->z[i] - axisrz[1]));
                p->phi[i] += atan2(
                    Xin_xyz[0] * Xout_xyz[1] - Xin_xyz[1] * Xout_xyz[0],
                    Xin_xyz[0] * Xout_xyz[0] + Xin_xyz[1] * Xout_xyz[1]);
            }

            /* Error handling */
            if (errflag)
            {
                p->err[i] = errflag;
                p->running[i] = 0;
            }
        }
    }
}
