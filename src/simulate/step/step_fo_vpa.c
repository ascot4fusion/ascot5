/**
 * @file step_fo_vpa.c
 * @brief Calculate a full orbit step for a struct of particles with VPA
 **/
#include <math.h>
#include <stdio.h>
#include "../../offload.h"
#include "../../ascot5.h"
#include "../../math.h"
#include "../../consts.h"
#include "../../physlib.h"
#include "../../error.h"
#include "../../B_field.h"
#include "../../E_field.h"
#include "../../boozer.h"
#include "../../mhd.h"
#include "../../particle.h"
#include "step_fo_vpa.h"

/**
 * @brief Integrate a full orbit step for a struct of particles with VPA
 *
 * The integration is performed for a struct of NSIMD particles using the
 * volume preserving algorithm (Boris method for relativistic particles) see
 * Zhang 2015.
 *
 * This algorithm is valid for neutral particles as well, in which case the
 * motion reduces to ballistic motion where momentum remains constant.
 *
 * @param p particle_simd_fo struct that will be updated
 * @param h pointer to array containing time steps
 * @param Bdata pointer to magnetic field data
 * @param Edata pointer to electric field data
 */
void step_fo_vpa(particle_simd_fo* p, real* h, B_field_data* Bdata,
                 E_field_data* Edata, int n_running_ref) {
    GPU_DATA_IS_MAPPED(h[0:p->n_mrk])
    GPU_PARALLEL_LOOP_ALL_LEVELS
      for(int i = 0; i < n_running_ref; i++) {
        if(p->running[i]) {
            a5err errflag = 0;

            real R0   = p->r[i];
            real z0   = p->z[i];
            real t0   = p->time[i];
            real mass = p->mass[i];

            /* Convert velocity to cartesian coordinates */
            real prpz[3] = {p->p_r[i], p->p_phi[i], p->p_z[i]};
            real pxyz[3];
            math_vec_rpz2xyz(prpz, pxyz, p->phi[i]);

            real posrpz[3] = {p->r[i], p->phi[i], p->z[i]};
            real posxyz0[3],posxyz[3];
            math_rpz2xyz(posrpz, posxyz0);

            /* Take a half step and evaluate fields at that position */
            real gamma = physlib_gamma_pnorm(mass, math_norm(pxyz));
            posxyz[0] = posxyz0[0] + pxyz[0] * h[i] / (2.0 * gamma * mass);
            posxyz[1] = posxyz0[1] + pxyz[1] * h[i] / (2.0 * gamma * mass);
            posxyz[2] = posxyz0[2] + pxyz[2] * h[i] / (2.0 * gamma * mass);

            math_xyz2rpz(posxyz, posrpz);

            real Brpz[3];
            real Erpz[3];
            if(!errflag) {
                errflag = B_field_eval_B(Brpz, posrpz[0], posrpz[1], posrpz[2],
                                         t0 + h[i]/2.0, Bdata);
            }
            if(!errflag) {
                errflag = E_field_eval_E(Erpz, posrpz[0], posrpz[1], posrpz[2],
                                         t0 + h[i]/2.0, Edata, Bdata);
            }

            real fposxyz[3]; // final position in cartesian coordinates

            if(!errflag) {
                /* Electromagnetic fields to cartesian coordinates */
                real Bxyz[3];
                real Exyz[3];

                math_vec_rpz2xyz(Brpz, Bxyz, posrpz[1]);
                math_vec_rpz2xyz(Erpz, Exyz, posrpz[1]);

                /* Evaluate helper variable pminus */
                real pminus[3];
                real sigma = p->charge[i]*h[i]/(2*p->mass[i]*CONST_C);
                pminus[0] = pxyz[0] / (mass * CONST_C) + sigma * Exyz[0];
                pminus[1] = pxyz[1] / (mass * CONST_C) + sigma * Exyz[1];
                pminus[2] = pxyz[2] / (mass * CONST_C) + sigma * Exyz[2];

                /* Second helper variable pplus*/
                real d = (p->charge[i]*h[i]/(2*p->mass[i])) /
                    sqrt( 1 + math_dot(pminus,pminus) );
                real d2 = d*d;

                real Bhat[9] = {       0,  Bxyz[2], -Bxyz[1],
                                -Bxyz[2],        0,  Bxyz[0],
                                 Bxyz[1], -Bxyz[0],        0};
                real Bhat2[9];
                math_matmul(Bhat, Bhat, 3, 3, 3, Bhat2);

                real B2 = Bxyz[0]*Bxyz[0] + Bxyz[1]*Bxyz[1] + Bxyz[2]*Bxyz[2];

                real A[9];
                for(int j=0; j<9; j++) {
                    A[j] = (Bhat[j] + d*Bhat2[j]) * (2.0*d/(1.0+d2*B2));
                }

                real pplus[3];
                math_matmul(pminus, A, 1, 3, 3, pplus);

                /* Take the step */
                real pfinal[3];
                pfinal[0] = pminus[0] + pplus[0] + sigma*Exyz[0];
                pfinal[1] = pminus[1] + pplus[1] + sigma*Exyz[1];
                pfinal[2] = pminus[2] + pplus[2] + sigma*Exyz[2];

                pxyz[0] = pfinal[0] * mass * CONST_C;
                pxyz[1] = pfinal[1] * mass * CONST_C;
                pxyz[2] = pfinal[2] * mass * CONST_C;
            }

            gamma = physlib_gamma_pnorm(mass, math_norm(pxyz));
            fposxyz[0] = posxyz[0] + h[i] * pxyz[0] / (2.0 * gamma * mass);
            fposxyz[1] = posxyz[1] + h[i] * pxyz[1] / (2.0 * gamma * mass);
            fposxyz[2] = posxyz[2] + h[i] * pxyz[2] / (2.0 * gamma * mass);

            if(!errflag) {
                /* Back to cylindrical coordinates */
                p->r[i] = sqrt(fposxyz[0]*fposxyz[0]+fposxyz[1]*fposxyz[1]);

                /* phi is evaluated like this to make sure it is cumulative */
                p->phi[i] += atan2(
                    posxyz0[0] * fposxyz[1] - posxyz0[1] * fposxyz[0],
                    posxyz0[0] * fposxyz[0] + posxyz0[1] * fposxyz[1] );
                p->z[i] = fposxyz[2];

                real cosp = cos(p->phi[i]);
                real sinp = sin(p->phi[i]);
                p->p_r[i]   =  pxyz[0] * cosp + pxyz[1] * sinp;
                p->p_phi[i] = -pxyz[0] * sinp + pxyz[1] * cosp;
                p->p_z[i]   =  pxyz[2];
            }

            /* Evaluate magnetic field (and gradient) and rho at new position */
            real BdBrpz[15];
            real psi[1];
            real rho[2];
            if(!errflag) {
                errflag = B_field_eval_B_dB(BdBrpz, p->r[i], p->phi[i], p->z[i],
                                            t0 + h[i], Bdata);
            }
            if(!errflag) {
                errflag = B_field_eval_psi(psi, p->r[i], p->phi[i], p->z[i],
                                           t0 + h[i], Bdata);
            }
            if(!errflag) {
                errflag = B_field_eval_rho(rho, psi[0], Bdata);
            }

            if(!errflag) {
                p->B_r[i]        = BdBrpz[0];
                p->B_r_dr[i]     = BdBrpz[1];
                p->B_r_dphi[i]   = BdBrpz[2];
                p->B_r_dz[i]     = BdBrpz[3];

                p->B_phi[i]      = BdBrpz[4];
                p->B_phi_dr[i]   = BdBrpz[5];
                p->B_phi_dphi[i] = BdBrpz[6];
                p->B_phi_dz[i]   = BdBrpz[7];

                p->B_z[i]        = BdBrpz[8];
                p->B_z_dr[i]     = BdBrpz[9];
                p->B_z_dphi[i]   = BdBrpz[10];
                p->B_z_dz[i]     = BdBrpz[11];

                p->rho[i] = rho[0];

                /* Evaluate phi and theta angles so that they are cumulative */
                real axisrz[2];
                errflag = B_field_get_axis_rz(axisrz, Bdata, p->phi[i]);
                p->theta[i] += atan2(   (R0-axisrz[0]) * (p->z[i]-axisrz[1])
                                      - (z0-axisrz[1]) * (p->r[i]-axisrz[0]),
                                        (R0-axisrz[0]) * (p->r[i]-axisrz[0])
                                      + (z0-axisrz[1]) * (p->z[i]-axisrz[1]) );
            }

            /* Error handling */
            if(errflag) {
                p->err[i] = errflag;
                p->running[i] = 0;
            }
        }
    }
}


/**
 * @brief Integrate a full orbit step with VPA and MHd modes present.
 *
 * Same as previous method but with MHD present.
 *
 * @param p particle_simd_fo struct that will be updated
 * @param h pointer to array containing time steps
 * @param Bdata pointer to magnetic field data
 * @param Edata pointer to electric field data
 * @param boozer pointer to boozer data
 * @param mhd pointer to MHD data
 */
void step_fo_vpa_mhd(particle_simd_fo* p, real* h, B_field_data* Bdata,
                     E_field_data* Edata, boozer_data* boozer, mhd_data* mhd) {

    int i;
    /* Following loop will be executed simultaneously for all i */
    #pragma omp simd  aligned(h : 64)
    for(i = 0; i < NSIMD; i++) {
        if(p->running[i]) {
            a5err errflag = 0;

            real R0   = p->r[i];
            real z0   = p->z[i];
            real t0   = p->time[i];
            real mass = p->mass[i];

            /* Convert velocity to cartesian coordinates */
            real prpz[3] = {p->p_r[i], p->p_phi[i], p->p_z[i]};
            real pxyz[3];
            math_vec_rpz2xyz(prpz, pxyz, p->phi[i]);

            real posrpz[3] = {p->r[i], p->phi[i], p->z[i]};
            real posxyz0[3],posxyz[3];
            math_rpz2xyz(posrpz,posxyz0);

            /* Take a half step and evaluate fields at that position */
            real gamma = physlib_gamma_pnorm(mass, math_norm(pxyz));
            posxyz[0] = posxyz0[0] + pxyz[0] * h[i] / (2 * gamma * mass);
            posxyz[1] = posxyz0[1] + pxyz[1] * h[i] / (2 * gamma * mass);
            posxyz[2] = posxyz0[2] + pxyz[2] * h[i] / (2 * gamma * mass);

            math_xyz2rpz(posxyz,posrpz);

            real Brpz[3];
            real Erpz[3];
            if(!errflag) {
                errflag = B_field_eval_B(Brpz, posrpz[0], posrpz[1], posrpz[2],
                                         t0 + h[i]/2, Bdata);
            }
            if(!errflag) {
                errflag = E_field_eval_E(Erpz, posrpz[0], posrpz[1], posrpz[2],
                                         t0 + h[i]/2, Edata, Bdata);
            }

            real pert[7];
            int pertonly = 0;
            if(!errflag) {
                errflag = mhd_perturbations(
                    pert, posrpz[0], posrpz[1], posrpz[2], t0 + h[i]/2,
                    pertonly, MHD_INCLUDE_ALL, boozer, mhd, Bdata);
            }
            Brpz[0] = pert[0];
            Brpz[1] = pert[1];
            Brpz[2] = pert[2];
            Erpz[0] += pert[3];
            Erpz[1] += pert[4];
            Erpz[2] += pert[5];

            real fposxyz[3]; // final position in cartesian coordinates

            if(!errflag) {
                /* Electromagnetic fields to cartesian coordinates */
                real Bxyz[3];
                real Exyz[3];

                math_vec_rpz2xyz(Brpz, Bxyz, posrpz[1]);
                math_vec_rpz2xyz(Erpz, Exyz, posrpz[1]);

                /* Evaluate helper variable pminus */
                real pminus[3];
                real sigma = p->charge[i]*h[i]/(2*p->mass[i]*CONST_C);
                pminus[0] = pxyz[0] / (mass * CONST_C) + sigma * Exyz[0];
                pminus[1] = pxyz[1] / (mass * CONST_C) + sigma * Exyz[1];
                pminus[2] = pxyz[2] / (mass * CONST_C) + sigma * Exyz[2];

                /* Second helper variable pplus*/
                real d = (p->charge[i]*h[i]/(2*p->mass[i])) /
                    sqrt( 1 + math_dot(pminus,pminus) );
                real d2 = d*d;

                real Bhat[9] = {       0,  Bxyz[2], -Bxyz[1],
                                -Bxyz[2],        0,  Bxyz[0],
                        Bxyz[1], -Bxyz[0],        0};
                real Bhat2[9];
                math_matmul(Bhat, Bhat, 3, 3, 3, Bhat2);

                real B2 = Bxyz[0]*Bxyz[0] + Bxyz[1]*Bxyz[1] + Bxyz[2]*Bxyz[2];

                real A[9];
                for(int j=0; j<9; j++) {
                    A[j] = (Bhat[j] + d*Bhat2[j]) * (2.0*d/(1+d2*B2));
                }

                real pplus[3];
                math_matmul(pminus, A, 1, 3, 3, pplus);

                /* Take the step */
                real pfinal[3];
                pfinal[0] = pminus[0] + pplus[0] + sigma*Exyz[0];
                pfinal[1] = pminus[1] + pplus[1] + sigma*Exyz[1];
                pfinal[2] = pminus[2] + pplus[2] + sigma*Exyz[2];

                pxyz[0] = pfinal[0] * mass * CONST_C;
                pxyz[1] = pfinal[1] * mass * CONST_C;
                pxyz[2] = pfinal[2] * mass * CONST_C;
            }

            gamma = physlib_gamma_pnorm(mass, math_norm(pxyz));
            fposxyz[0] = posxyz[0] + h[i] * pxyz[0] / (2 * gamma * mass);
            fposxyz[1] = posxyz[1] + h[i] * pxyz[1] / (2 * gamma * mass);
            fposxyz[2] = posxyz[2] + h[i] * pxyz[2] / (2 * gamma * mass);

            if(!errflag) {
                /* Back to cylindrical coordinates */
                p->r[i] = sqrt(fposxyz[0]*fposxyz[0]+fposxyz[1]*fposxyz[1]);

                /* phi is evaluated like this to make sure it is cumulative */
                p->phi[i] += atan2(
                    posxyz0[0] * fposxyz[1] - posxyz0[1] * fposxyz[0],
                    posxyz0[0] * fposxyz[0] + posxyz0[1] * fposxyz[1] );
                p->z[i] = fposxyz[2];

                real cosp = cos(p->phi[i]);
                real sinp = sin(p->phi[i]);
                p->p_r[i]   =  pxyz[0] * cosp + pxyz[1] * sinp;
                p->p_phi[i] = -pxyz[0] * sinp + pxyz[1] * cosp;
                p->p_z[i]   =  pxyz[2];
            }

            /* Evaluate magnetic field (and gradient) and rho at new position */
            real BdBrpz[15];
            real psi[1];
            real rho[2];
            if(!errflag) {
                errflag = B_field_eval_B_dB(BdBrpz, p->r[i], p->phi[i], p->z[i],
                                            t0 + h[i], Bdata);
            }
            if(!errflag) {
                errflag = B_field_eval_psi(psi, p->r[i], p->phi[i], p->z[i],
                                           t0 + h[i], Bdata);
            }
            if(!errflag) {
                errflag = B_field_eval_rho(rho, psi[0], Bdata);
            }

            if(!errflag) {
                p->B_r[i]        = BdBrpz[0];
                p->B_r_dr[i]     = BdBrpz[1];
                p->B_r_dphi[i]   = BdBrpz[2];
                p->B_r_dz[i]     = BdBrpz[3];

                p->B_phi[i]      = BdBrpz[4];
                p->B_phi_dr[i]   = BdBrpz[5];
                p->B_phi_dphi[i] = BdBrpz[6];
                p->B_phi_dz[i]   = BdBrpz[7];

                p->B_z[i]        = BdBrpz[8];
                p->B_z_dr[i]     = BdBrpz[9];
                p->B_z_dphi[i]   = BdBrpz[10];
                p->B_z_dz[i]     = BdBrpz[11];

                p->rho[i] = rho[0];

                /* Evaluate phi and theta angles so that they are cumulative */
                real axisrz[2];
                errflag = B_field_get_axis_rz(axisrz, Bdata, p->phi[i]);
                p->theta[i] += atan2(   (R0-axisrz[0]) * (p->z[i]-axisrz[1])
                                      - (z0-axisrz[1]) * (p->r[i]-axisrz[0]),
                                        (R0-axisrz[0]) * (p->r[i]-axisrz[0])
                                      + (z0-axisrz[1]) * (p->z[i]-axisrz[1]) );
            }

            /* Error handling */
            if(errflag) {
                p->err[i] = errflag;
                p->running[i] = 0;
            }
        }
    }
}
