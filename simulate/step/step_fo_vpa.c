/**
 * @file step_fo_lf.c
 * @brief Calculate a full orbit step for a struct of particles with leap-frog
 **/
#include <math.h>
#include <stdio.h>
#include "../../ascot5.h"
#include "../../math.h"
#include "../../consts.h"
#include "../../error.h"
#include "../../B_field.h"
#include "../../E_field.h"
#include "../../particle.h"
#include "step_fo_vpa.h"

/**
 * @brief Integrate a full orbit step for a struct of particles with VPA
 *
 * The integration is performed for a struct of NSIMD particles using the
 * volume preserving algorithm (Boris method for relativistic particles) see Zhang 2015.
 *
 * @param p particle_simd_fo struct that will be updated
 * @param h pointer to array containing time steps
 * @param Bdata pointer to magnetic field data
 * @param Edata pointer to electric field data
 */
void step_fo_vpa(particle_simd_fo* p, real* h, B_field_data* Bdata, E_field_data* Edata) {

    int i;
    /* Following loop will be executed simultaneously for all i */
#pragma omp simd  aligned(h : 64)
    for(i = 0; i < NSIMD; i++) {
        if(p->running[i]) {
            a5err errflag = 0;

            real R0   = p->r[i];
            real z0   = p->z[i];

            /* Convert velocity to cartesian coordinates */
            real vrpz[3] = {p->rdot[i], p->phidot[i]*p->r[i], p->zdot[i]};
            real vxyz[3];
            math_vec_rpz2xyz(vrpz, vxyz, p->phi[i]);

            real posrpz[3] = {p->r[i], p->phi[i], p->z[i]};
            real posxyz0[3],posxyz[3];
            math_rpz2xyz(posrpz,posxyz0);

            /* Take a half step and evaluate fields at that position */
            posxyz[0] = posxyz0[0] + vxyz[0]*h[i]/2;
            posxyz[1] = posxyz0[1] + vxyz[1]*h[i]/2;
            posxyz[2] = posxyz0[2] + vxyz[2]*h[i]/2;

            math_xyz2rpz(posxyz,posrpz);

            real Brpz[3];
            real Erpz[3];
            if(!errflag) {errflag = B_field_eval_B(Brpz, posrpz[0], posrpz[1], posrpz[2], Bdata);}
            if(!errflag) {errflag = E_field_eval_E(Erpz, posrpz[0], posrpz[1], posrpz[2], Edata, Bdata);}


            real fposxyz[3]; // final position in cartesian coordinates

            if(!errflag) {
                /* Electromagnetic fields to cartesian coordinates */
                real Bxyz[3];
                real Exyz[3];

                math_vec_rpz2xyz(Brpz, Bxyz, posrpz[1]);
                math_vec_rpz2xyz(Erpz, Exyz, posrpz[1]);

                /* Evaluate helper variable pminus */
                real pminus[3];
                real vnorm = math_norm(vxyz);
                real gamma = sqrt(1 / ( (1 - vnorm/CONST_C)*(1 + vnorm/CONST_C) ));
                real sigma = p->charge[i]*h[i]/(2*p->mass[i]*CONST_C);
                pminus[0] = gamma*vxyz[0]/(CONST_C) + sigma*Exyz[0];
                pminus[1] = gamma*vxyz[1]/(CONST_C) + sigma*Exyz[1];
                pminus[2] = gamma*vxyz[2]/(CONST_C) + sigma*Exyz[2];

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
                gamma = sqrt(1/(1 + math_dot(pfinal,pfinal)));

                vxyz[0] = pfinal[0]*CONST_C*gamma;
                vxyz[1] = pfinal[1]*CONST_C*gamma;
                vxyz[2] = pfinal[2]*CONST_C*gamma;

                fposxyz[0] = posxyz[0] + h[i]*vxyz[0]/2;
                fposxyz[1] = posxyz[1] + h[i]*vxyz[1]/2;
                fposxyz[2] = posxyz[2] + h[i]*vxyz[2]/2;
            }

            /* Test that the results are reasonable */
            if(!errflag && ( posxyz[0] == 0 && posxyz[1] == 0 && posxyz[2] == 0 ))       {errflag = error_raise(ERR_UNPHYSICAL_FO, __LINE__);}
            if(!errflag && (vxyz[0]*vxyz[0]+vxyz[1]*vxyz[1]+vxyz[2]*vxyz[2]) > CONST_C2) {errflag = error_raise(ERR_UNPHYSICAL_FO, __LINE__);}

            if(!errflag) {
                /* Back to cylindrical coordinates */
                p->r[i] = sqrt(fposxyz[0]*fposxyz[0]+fposxyz[1]*fposxyz[1]);

                // We need to evaluate phi like this to make sure it is cumulative
                p->phi[i] += atan2( posxyz0[0] * fposxyz[1] - posxyz0[1] * fposxyz[0],
                                    posxyz0[0] * fposxyz[0] + posxyz0[1] * fposxyz[1] );
                p->z[i] = fposxyz[2];

                p->rdot[i] = vxyz[0] * cos(p->phi[i]) + vxyz[1] * sin(p->phi[i]);
                p->phidot[i] = ( -vxyz[0] * sin(p->phi[i]) + vxyz[1] * cos(p->phi[i]) ) / p->r[i];
                p->zdot[i] = vxyz[2];
            }

            /* Evaluate magnetic field (and gradient) and rho at new position */
            real BdBrpz[12];
            real psi[1];
            real rho[1];
            if(!errflag) {errflag = B_field_eval_B_dB(BdBrpz, p->r[i], p->phi[i], p->z[i], Bdata);}
            if(!errflag) {errflag = B_field_eval_psi(psi, p->r[i], p->phi[i], p->z[i], Bdata);}
            if(!errflag) {errflag = B_field_eval_rho(rho, psi[0], Bdata);}

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

                /* Evaluate phi and pol angles so that they are cumulative */
                real axis_r = B_field_get_axis_r(Bdata, p->phi[i]);
                real axis_z = B_field_get_axis_z(Bdata, p->phi[i]);
                p->pol[i] += atan2( (R0-axis_r) * (p->z[i]-axis_z) - (z0-axis_z) * (p->r[i]-axis_r),
                             (R0-axis_r) * (p->r[i]-axis_r) + (z0-axis_z) * (p->z[i]-axis_z) );
            }

            /* Error handling */
            if(errflag) {
                p->err[i] = error_module(errflag, ERRMOD_ORBSTEP);
                p->running[i] = 0;
            }
        }
    }
}
