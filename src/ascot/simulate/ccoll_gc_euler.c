/**
 * Implements Euler-Maruyama integrator for collision operator in GC picture
 * (see coulomb_collisions.h).
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
    (void)mdata;
    /* Get plasma information before going to the  SIMD loop */
    size_t n_species = Plasma_get_n_species(plasma);
    const real *qb = Plasma_get_species_charge(plasma);
    const real *mb = Plasma_get_species_mass(plasma);

    GPU_DATA_IS_MAPPED(h[0:p->size], rnd[0:3*p->size])
    GPU_PARALLEL_LOOP_ALL_LEVELS
    for (size_t i = 0; i < p->size; i++)
    {
        if (p->running[i])
        {
            err_t errflag = 0;

            /* Initial (R,z) position and magnetic field are needed for later */
            real Brpz[3] = {p->br[i], p->bphi[i], p->bz[i]};
            real Bnorm = math_norm(Brpz);
            real Bxyz[3];
            math_vec_rpz2xyz(Brpz, Bxyz, p->phi[i]);
            real R0 = p->r[i];
            real z0 = p->z[i];

            /* Move guiding center to (x, y, z, vnorm, xi) coordinates */
            real vin, pin, vflow, vpar, vperp2, xiin, Xin_xyz[3];
            Xin_xyz[0] = p->r[i] * cos(p->phi[i]);
            Xin_xyz[1] = p->r[i] * sin(p->phi[i]);
            Xin_xyz[2] = p->z[i];
            if (!errflag)
            {
                errflag = Plasma_eval_flow(
                    &vflow, p->rho[i], p->r[i], p->phi[i], p->z[i], p->time[i],
                    plasma);
            }
            pin  = physlib_gc_p(p->mass, p->mu[i], p->ppar[i], Bnorm);
            xiin = physlib_gc_xi(p->mass, p->mu[i], p->ppar[i], Bnorm);
            vin = physlib_vnorm_pnorm(p->mass, pin);
            vpar = xiin * vin;
            vperp2 = (1 - xiin * xiin) * vin * vin;
            vin = sqrt((vpar - vflow) * (vpar - vflow) + vperp2);
            xiin =  (vpar - vflow) / vin;

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
                clogab, p->mass, p->charge[i] * CONST_E, vin, n_species, mb, qb, nb,
                Tb);

            /* Evaluate collision coefficients and sum them for each *
             * species                                               */
            real gyrofreq =
                phys_gyrofreq_pnorm(p->mass, p->charge[i] * CONST_E, pin, Bnorm);
            real K = 0, Dpara = 0, nu = 0, DX = 0;
            GPU_SEQUENTIAL_LOOP
            for (size_t j = 0; j < n_species; j++)
            {
                real vb = sqrt(2 * Tb[j] / mb[j]);
                real x = vin / vb;
                real mufun[3];
                mccc_coefs_mufun(mufun, x); // eq. 2.83 PhD Hirvijoki

                real Qb = mccc_coefs_Q(
                    p->mass, p->charge[i] * CONST_E, mb[j], qb[j], nb[j], vb,
                    clogab[j], mufun[0]);
                real Dparab = mccc_coefs_Dpara(
                    p->mass, p->charge[i] * CONST_E, vin, qb[j], nb[j], vb, clogab[j],
                    mufun[0]);
                real Dperpb = mccc_coefs_Dperp(
                    p->mass, p->charge[i] * CONST_E, vin, qb[j], nb[j], vb, clogab[j],
                    mufun[1]);
                real dDparab = mccc_coefs_dDpara(
                    p->mass, p->charge[i] * CONST_E, vin, qb[j], nb[j], vb, clogab[j],
                    mufun[0], mufun[2]);

                K += mccc_coefs_K(vin, Dparab, dDparab, Qb);
                Dpara += Dparab;
                nu += mccc_coefs_nu(vin, Dperpb); // eq.41
                DX += mccc_coefs_DX(xiin, Dparab, Dperpb, gyrofreq);
            }

            /* Evaluate collisions */
            real sdt = sqrt(h[i]);
            real dW[5];
            dW[0]=sdt*rnd[0*p->size + i]; // For X_1
            dW[1]=sdt*rnd[1*p->size + i]; // For X_2
            dW[2]=sdt*rnd[2*p->size + i]; // For X_3
            dW[3]=sdt*rnd[3*p->size + i]; // For v
            dW[4]=sdt*rnd[4*p->size + i]; // For xi

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
            real cutoff = MCCC_CUTOFF * sqrt(Tb[0] / p->mass);
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

            vpar = xiout * vout;
            vperp2 = (1 - xiout * xiout) * vout * vout;
            vout = sqrt((vpar + vflow) * (vpar + vflow) + vperp2);
            xiout =  (vpar + vflow) / vout;
            real pout = physlib_pnorm_vnorm(p->mass, vout);

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
                p->br[i] = B_dB[0];
                p->dbrdr[i] = B_dB[3];
                p->dbrdphi[i] = B_dB[4];
                p->dbrdz[i] = B_dB[5];

                p->bphi[i] = B_dB[1];
                p->dbphidr[i] = B_dB[6];
                p->dbphidphi[i] = B_dB[7];
                p->dbphidz[i] = B_dB[8];

                p->bz[i] = B_dB[2];
                p->dbzdr[i] = B_dB[9];
                p->dbzdphi[i] = B_dB[10];
                p->dbzdz[i] = B_dB[11];

                p->rho[i] = rho[0];

                Bnorm = math_normc(B_dB[0], B_dB[1], B_dB[2]);

                p->r[i] = Xout_rpz[0];
                p->z[i] = Xout_rpz[2];

                p->ppar[i] = physlib_gc_ppar(pout, xiout);

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
