/**
 * @file mccc_gc_euler.c
 * @brief Euler-Maruyama integrator for collision operator in GC picture.
 */
#include <math.h>
#include "../../ascot5.h"
#include "../../consts.h"
#include "../../math.h"
#include "../../physlib.h"
#include "../../error.h"
#include "../../particle.h"
#include "../../B_field.h"
#include "../../plasma.h"
#include "../../random.h"
#include "mccc_coefs.h"
#include "mccc.h"

/**
 * @brief Integrate collisions for one time-step
 *
 * @param p gc struct
 * @param h time-steps for NSIMD markers
 * @param Bdata pointer to magnetic field
 * @param pdata pointer to plasma data
 * @param rdata pointer to random-generator data
 * @param mdata pointer to collision data struct
 */
void mccc_gc_euler(particle_simd_gc* p, real* h, B_field_data* Bdata,
                   plasma_data* pdata, random_data* rdata, mccc_data* mdata) {

    /* Generate random numbers and get plasma information before going to the *
     * SIMD loop                                                              */
    real rnd[5*NSIMD];
    random_normal_simd(rdata, 5*NSIMD, rnd);

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
            real vin, xiin, Xin_xyz[3];
            vin  = physlib_gc_v(p->mass[i], p->mu[i], p->vpar[i], Bnorm);
            xiin = physlib_gc_xi(p->mass[i], p->mu[i], p->vpar[i], Bnorm);
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
            real gyrofreq = phys_gyrofreq_vnorm(p->mass[i], p->charge[i], vin,
                                                Bnorm);
            real K = 0, Dpara = 0, nu = 0, DX = 0;
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

                K     += mccc_coefs_K(vin, Dparab, dDparab, Qb);
                Dpara += Dparab;
                nu    += mccc_coefs_nu(vin, Dperpb);
                DX    += mccc_coefs_DX(xiin, Dparab, Dperpb, gyrofreq);
            }

            /* Evaluate collisions */
            real sdt = sqrt(h[i]);
            real dW[5];
            dW[0]=sdt*rnd[0*NSIMD + i]; // For X_1
            dW[1]=sdt*rnd[1*NSIMD + i]; // For X_2
            dW[2]=sdt*rnd[2*NSIMD + i]; // For X_3
            dW[3]=sdt*rnd[3*NSIMD + i]; // For v
            dW[4]=sdt*rnd[4*NSIMD + i]; // For xi

            real bhat[3];
            math_unit(Bxyz, bhat);

            real k1 = sqrt(2*DX);
            real k2 = math_dot(bhat, dW);

            real vout, xiout, Xout_xyz[3];
            Xout_xyz[0] = Xin_xyz[0] + k1 * ( dW[0] - k2 * bhat[0] );
            Xout_xyz[1] = Xin_xyz[1] + k1 * ( dW[1] - k2 * bhat[1] );
            Xout_xyz[2] = Xin_xyz[2] + k1 * ( dW[2] - k2 * bhat[2] );
            vout  = vin + K*h[i] + sqrt( 2 * Dpara ) * dW[3];
            xiout = xiin - xiin*nu*h[i] + sqrt(( 1 - xiin*xiin ) * nu) * dW[4];

            /* Enforce boundary conditions */
            real cutoff = MCCC_CUTOFF * sqrt( Tb[0] / p->mass[i] );
            if(vout < cutoff){
                vout = 2 * cutoff - vout;
            }

            if(fabs(xiout) > 1){
                xiout = ( (xiout > 0) - (xiout < 0) )
                        * ( 2 - fabs( xiout ) );
            }

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

            /* Back to cylindrical coordinates */
            real Xout_rpz[3];
            math_xyz2rpz(Xout_xyz, Xout_rpz);

            /* Evaluate magnetic field (and gradient) and rho at new position */
            real B_dB[15], psi[1], rho[1];
            if(!errflag) {
                errflag = B_field_eval_B_dB(B_dB, Xout_rpz[0], Xout_rpz[1],
                                            Xout_rpz[2], p->time[i] + h[i],
                                            Bdata);
            }
            if(!errflag) {
                errflag = B_field_eval_psi(psi, Xout_rpz[0], Xout_rpz[1],
                                           Xout_rpz[2], p->time[i] + h[i],
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
                p->vpar[i] = physlib_gc_vpar(vout, xiout);
                p->mu[i]   = physlib_gc_mu(p->mass[i], vout, xiout, Bnorm);

                /* Evaluate phi and pol angles so that they are cumulative */
                real axis_r = B_field_get_axis_r(Bdata, p->phi[i]);
                real axis_z = B_field_get_axis_z(Bdata, p->phi[i]);
                p->pol[i] += atan2(   (R0-axis_r) * (p->z[i]-axis_z)
                                    - (z0-axis_z) * (p->r[i]-axis_r),
                                      (R0-axis_r) * (p->r[i]-axis_r)
                                    + (z0-axis_z) * (p->z[i]-axis_z) );
                p->phi[i] += atan2(   Xin_xyz[0] * Xout_xyz[1]
                                    - Xin_xyz[1] * Xout_xyz[0],
                                      Xin_xyz[0] * Xout_xyz[0]
                                    + Xin_xyz[1] * Xout_xyz[1] );
            }

            /* Error handling */
            if(errflag) {
                p->err[i]     = errflag;
                p->running[i] = 0;
            }
        }
    }
}
