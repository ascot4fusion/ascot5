/**
 * @file mccc_gc_milstein.c
 * @brief Milstein integrator for collision operator in GC picture.
 */
#include <math.h>
#include <float.h>
#include "../../ascot5.h"
#include "../../consts.h"
#include "../../math.h"
#include "../../physlib.h"
#include "../../error.h"
#include "../../particle.h"
#include "../../B_field.h"
#include "../../plasma.h"
#include "../../random.h"
#include "mccc_wiener.h"
#include "mccc_coefs.h"
#include "mccc.h"

/**
 * @brief Integrate collisions for one time-step
 *
 * @param p pointer to gc simd struct
 * @param hin time-steps for NSIMD markers
 * @param hout suggestions for the next timesteps for NSIMD markers
 * @param tol relative error tolerance
 * @param w array holding wiener structs for NSIMD markers
 * @param Bdata pointer to magnetic field data
 * @param pdata pointer to plasma data
 * @param mdata pointer to collision data struct
 * @param rnd array of normally distributed random numbers used to resolve
 *        collisions. Values for marker i are rnd[i*NSIMD + j]
 */
void mccc_gc_milstein(particle_simd_gc* p, real* hin, real* hout, real tol,
                      mccc_wienarr* w, B_field_data* Bdata, plasma_data* pdata,
                      mccc_data* mdata, real* rnd) {

    /* Get plasma information before going to the  SIMD loop */
    int n_species  = plasma_get_n_species(pdata);
    const real* qb = plasma_get_species_charge(pdata);
    const real* mb = plasma_get_species_mass(pdata);

    #pragma omp simd
    for(int i = 0; i < NSIMD; i++) {
        if(p->running[i]) {
            a5err errflag = 0;

            /* Initial (R,z) position and magnetic field are needed for later */
            real Brpz[3] = {p->B_r[i], p->B_phi[i], p->B_z[i]};
            real Bnorm   = math_norm(Brpz);
            real Bxyz[3];
            math_vec_rpz2xyz(Brpz, Bxyz, p->phi[i]);
            real R0   = p->r[i];
            real z0   = p->z[i];

            /* Move guiding center to (x, y, z, vnorm, xi) coordinates */
            real vin, pin, xiin, Xin_xyz[3];
            pin  = physlib_gc_p( p->mass[i], p->mu[i], p->ppar[i], Bnorm);
            xiin = physlib_gc_xi(p->mass[i], p->mu[i], p->ppar[i], Bnorm);
            vin  = physlib_vnorm_pnorm(p->mass[i], pin);
            Xin_xyz[0] = p->r[i] * cos(p->phi[i]);
            Xin_xyz[1] = p->r[i] * sin(p->phi[i]);
            Xin_xyz[2] = p->z[i];

            /* Evaluate plasma density and temperature */
            real nb[MAX_SPECIES], Tb[MAX_SPECIES];
            if(!errflag) {
                errflag = plasma_eval_densandtemp(nb, Tb, p->rho[i],
                                                  p->r[i], p->phi[i], p->z[i],
                                                  p->time[i], pdata);
            }

            /* Coulomb logarithm */
            real clogab[MAX_SPECIES];
            mccc_coefs_clog(clogab, p->mass[i], p->charge[i], vin,
                            n_species, mb, qb, nb, Tb);

            /* Evaluate collision coefficients and sum them for each *
             * species                                               */
            real gyrofreq = phys_gyrofreq_pnorm(p->mass[i], p->charge[i], pin,
                                                Bnorm);
            real K = 0, Dpara = 0, dDpara = 0, dQ = 0, nu = 0, DX = 0;
            for(int j = 0; j < n_species; j++) {
                real vb = sqrt( 2 * Tb[j] / mb[j] );
                real x  = vin / vb;
                real mufun[3];
                mccc_coefs_mufun(mufun, x, mdata);

                real Qb      = mccc_coefs_Q(p->mass[i], p->charge[i], mb[j],
                                            qb[j], nb[j], vb, clogab[j],
                                            mufun[0]);
                real Dparab  = mccc_coefs_Dpara(p->mass[i], p->charge[i], vin,
                                               qb[j], nb[j], vb, clogab[j],
                                               mufun[0]);
                real Dperpb  = mccc_coefs_Dperp(p->mass[i], p->charge[i], vin,
                                               qb[j], nb[j], vb, clogab[j],
                                               mufun[1]);
                real dDparab = mccc_coefs_dDpara(p->mass[i], p->charge[i], vin,
                                                 qb[j], nb[j], vb, clogab[j],
                                                 mufun[0], mufun[2]);

                K      += mccc_coefs_K(vin, Dparab, dDparab, Qb);
                dQ     += mccc_coefs_dQ(p->mass[i], p->charge[i], mb[j],
                                        qb[j], nb[j], vb, clogab[j], mufun[2]);
                Dpara  += Dparab;
                dDpara += dDparab;
                nu     += mccc_coefs_nu(vin, Dperpb);
                DX     += mccc_coefs_DX(xiin, Dparab, Dperpb, gyrofreq);
            }

            /* Generate Wiener process for this step */
            int tindex;
            if(!errflag) {
                errflag = mccc_wiener_generate(&w[i], w[i].time[0]+hin[i],
                                               &tindex, &rnd[i*5]);
            }
            real dW[5] = {0, 0, 0, 0, 0};
            if(!errflag) {
                dW[0] = w[i].wiener[tindex*5 + 0] - w[i].wiener[0]; // For X_1
                dW[1] = w[i].wiener[tindex*5 + 1] - w[i].wiener[1]; // For X_2
                dW[2] = w[i].wiener[tindex*5 + 2] - w[i].wiener[2]; // For X_3
                dW[3] = w[i].wiener[tindex*5 + 3] - w[i].wiener[3]; // For v
                dW[4] = w[i].wiener[tindex*5 + 4] - w[i].wiener[4]; // For xi
            }

            /* Evaluate collisions */

            real bhat[3];
            math_unit(Bxyz,bhat);

            real k1 = sqrt(2*DX);
            real k2 = math_dot(bhat, dW);

            real vout, xiout, Xout_xyz[3];
            Xout_xyz[0] = Xin_xyz[0] + k1 * ( dW[0] - k2 * bhat[0] );
            Xout_xyz[1] = Xin_xyz[1] + k1 * ( dW[1] - k2 * bhat[1] );
            Xout_xyz[2] = Xin_xyz[2] + k1 * ( dW[2] - k2 * bhat[2] );
            vout  = vin + K*hin[i] + sqrt( 2 * Dpara ) * dW[3]
                  + 0.5 * dDpara * ( dW[3]*dW[3] - hin[i] );
            xiout = xiin - xiin*nu*hin[i] + sqrt( ( 1 - xiin*xiin ) * nu )*dW[4]
                  - 0.5 * xiin * nu * ( dW[4]*dW[4] - hin[i] );

            /* Enforce boundary conditions */
            real cutoff = MCCC_CUTOFF * sqrt( Tb[0] / p->mass[i] );
            if(vout < cutoff){
                vout = 2 * cutoff - vout;
            }

            if(fabs(xiout) > 1){
                xiout = ( (xiout > 0) - (xiout < 0) )
                        * ( 2 - fabs( xiout ) );
            }

            /* Compute error estimates for drift and diffusion limits */

            // xi is limited to interval [-1, 1] but for v we need some value
            // to translate relative error to absolute error.
            real v0    = ( vin + fabs(K) * hin[i] + sqrt( 2*Dpara*hin[i] ) )
                         + DBL_EPSILON;
            real verr  = fabs( K*dQ ) / (2*tol*v0);
            real xierr = fabs( xiin*nu*nu ) / (2*tol);

            // kappa_k is error due to drift
            real kappa_k;
            if(verr > xierr){
                kappa_k = verr*hin[i]*hin[i];
            }
            else{
                kappa_k = xierr*hin[i]*hin[i];
            }

            // kappa_d is error due to diffusion (v and xi are both needed)
            real kappa_d0 = fabs(  dW[3]*dW[3]*dW[3]
                                 * dDpara*dDpara / sqrt( Dpara ) ) / (6*tol*v0);
            real kappa_d1 = sqrt( 1 - xiin*xiin ) * nu * sqrt( nu )
                          * fabs( dW[4] + sqrt( hin[i]/3 ) ) * hin[i] / (2*tol);

            /* Remove energy or pitch change or spatial diffusion from the    *
             * results if that is requested                                   */
            if(!mdata->include_energy) {
                vout = vin;
            }
            if(!mdata->include_pitch) {
                xiout = xiin;
            }
            if(!mdata->include_gcdiff) {
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
            if(!errflag) {
                errflag = B_field_eval_B_dB(B_dB, Xout_rpz[0], Xout_rpz[1],
                                            Xout_rpz[2], p->time[i] + hin[i],
                                            Bdata);
            }
            if(!errflag) {
                errflag = B_field_eval_psi(psi, Xout_rpz[0], Xout_rpz[1],
                                           Xout_rpz[2], p->time[i] + hin[i],
                                           Bdata);
            }
            if(!errflag) {
                errflag = B_field_eval_rho(rho, psi[0], Bdata);
            }

            if(!errflag) {
                /* Update marker coordinates at the new position */
                p->B_r[i]        = B_dB[0];
                p->B_r_dr[i]     = B_dB[1];
                p->B_r_dphi[i]   = B_dB[2];
                p->B_r_dz[i]     = B_dB[3];

                p->B_phi[i]      = B_dB[4];
                p->B_phi_dr[i]   = B_dB[5];
                p->B_phi_dphi[i] = B_dB[6];
                p->B_phi_dz[i]   = B_dB[7];

                p->B_z[i]        = B_dB[8];
                p->B_z_dr[i]     = B_dB[9];
                p->B_z_dphi[i]   = B_dB[10];
                p->B_z_dz[i]     = B_dB[11];

                p->rho[i] = rho[0];

                Bnorm = math_normc(B_dB[0], B_dB[4], B_dB[8]);

                p->r[i]    = Xout_rpz[0];
                p->z[i]    = Xout_rpz[2];
                p->ppar[i] = physlib_gc_ppar(pout, xiout);
                p->mu[i]   = physlib_gc_mu(p->mass[i], pout, xiout, Bnorm);

                /* Evaluate phi and theta angles so that they are cumulative */
                real axisrz[2];
                errflag  = B_field_get_axis_rz(axisrz, Bdata, p->phi[i]);
                p->theta[i] += atan2(   (R0-axisrz[0]) * (p->z[i]-axisrz[1])
                                      - (z0-axisrz[1]) * (p->r[i]-axisrz[0]),
                                        (R0-axisrz[0]) * (p->r[i]-axisrz[0])
                                      + (z0-axisrz[1]) * (p->z[i]-axisrz[1]) );
                p->phi[i] += atan2(   Xin_xyz[0] * Xout_xyz[1]
                                    - Xin_xyz[1] * Xout_xyz[0],
                                      Xin_xyz[0] * Xout_xyz[0]
                                    + Xin_xyz[1] * Xout_xyz[1] );
            }

            /* Check whether timestep was rejected and suggest next time step */

            if( kappa_k >= kappa_d0 && kappa_k >= kappa_d1 ) {
                /* Drift error dominates */
                hout[i] = 0.8 * hin[i] / sqrt( kappa_k );
            }
            else if( kappa_d0 >= kappa_k && kappa_d0 >= kappa_d1 ) {
                /* Velocity diffusion error dominates */
                hout[i] = 0.9 * hin[i] * pow( kappa_d0, -2.0/3.0 );
            }
            else {
                /* Pitch diffusion error dominates */
                hout[i] = 0.9 * hin[i] * pow( kappa_d1, -2.0/3.0 );
            }

            /* Negative value indicates time step was rejected*/
            if( kappa_k > 1 || kappa_d0 > 1 || kappa_d1 > 1 ){
                hout[i]   = -hout[i];
            }
            else if(hout[i] > 1.5*hin[i]) {
                /* Make sure we don't increase time step too much */
                hout[i] = 1.5*hin[i];
            }

            /* Error handling */
            if(errflag) {
                p->err[i]     = errflag;
                p->running[i] = 0;
            }
        }
    }
}
