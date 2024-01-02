/**
 * @file mccc_fo_euler.c
 * @brief Euler-Maruyama integrator for collision operator in FO picture.
 */
#include <math.h>
#include "../../ascot5.h"
#include "../../consts.h"
#include "../../math.h"
#include "../../error.h"
#include "../../physlib.h"
#include "../../particle.h"
#include "../../plasma.h"
#include "../../random.h"
#include "mccc_coefs.h"
#include "mccc.h"

/**
 * @brief Integrate collisions for one time-step
 *
 * @param p fo struct
 * @param h time-steps for NSIMD markers
 * @param pdata pointer to plasma data
 * @param rdata pointer to random-generator data
 * @param mdata pointer collision data struct
 */
void mccc_fo_euler(particle_simd_fo* p, real* h, plasma_data* pdata,
                   mccc_data* mdata, real rnd[3*NSIMD]) {

    /* Get plasma information before going to the  SIMD loop */
    int n_species  = plasma_get_n_species(pdata);
    const real* qb = plasma_get_species_charge(pdata);
    const real* mb = plasma_get_species_mass(pdata);

    #pragma omp simd
    for(int i = 0; i < NSIMD; i++) {
        if(p->running[i]) {
            a5err errflag = 0;

            /* These are needed twice to transform velocity to cartesian and *
             * back to cylindrical coordinates. Position does not change     */
            real sinphi = sin(p->phi[i]);
            real cosphi = cos(p->phi[i]);

            real pnorm = sqrt( p->p_r[i] * p->p_r[i] + p->p_phi[i] * p->p_phi[i]
                               + p->p_z[i] * p->p_z[i] );
            real gamma = physlib_gamma_pnorm(p->mass[i], pnorm);

            real vin_xyz[3];
            vin_xyz[0] = ( p->p_r[i] * cosphi - p->p_phi[i] * sinphi )
                / ( gamma * p->mass[i] );
            vin_xyz[1] = ( p->p_r[i] * sinphi + p->p_phi[i] * cosphi )
                / ( gamma * p->mass[i] );
            vin_xyz[2] = p->p_z[i] / ( gamma * p->mass[i] );
            real vin   = math_norm(vin_xyz);

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
            real F = 0, Dpara = 0, Dperp = 0;
            for(int j = 0; j < n_species; j++) {
                real vb = sqrt( 2 * Tb[j] / mb[j] );
                real x  = vin / vb;
                real mufun[3];
                mccc_coefs_mufun(mufun, x, mdata);

                F     += mccc_coefs_F(p->mass[i], p->charge[i], mb[j], qb[j],
                                      nb[j], vb, clogab[j], mufun[0]);
                Dpara += mccc_coefs_Dpara(p->mass[i], p->charge[i], vin, qb[j],
                                          nb[j], vb, clogab[j], mufun[0]);
                Dperp += mccc_coefs_Dperp(p->mass[i], p->charge[i], vin, qb[j],
                                          nb[j], vb, clogab[j], mufun[1]);
            }

            /* Evaluate collisions */
            real sdt = sqrt(h[i]);
            real dW[3];
            dW[0] = sdt * rnd[0*NSIMD + i];
            dW[1] = sdt * rnd[1*NSIMD + i];
            dW[2] = sdt * rnd[2*NSIMD + i];

            real vhat[3];
            math_unit(vin_xyz, vhat);

            real t1 = math_dot(vhat, dW);
            real k1 = F*h[i];
            real k2 = sqrt(2*Dpara)*t1;
            real k3 = sqrt(2*Dperp);

            real vout_xyz[3];
            vout_xyz[0] = vin_xyz[0] + k1*vhat[0] + k2*vhat[0]
                        + k3*(dW[0]  - t1*vhat[0]);
            vout_xyz[1] = vin_xyz[1] + k1*vhat[1] + k2*vhat[1]
                        + k3*(dW[1]  - t1*vhat[1]);
            vout_xyz[2] = vin_xyz[2] + k1*vhat[2] + k2*vhat[2]
                        + k3*(dW[2]  - t1*vhat[2]);

            /* Transform back to cylindrical coordinates.  */

            real vnorm = math_norm(vout_xyz);
            gamma = physlib_gamma_vnorm(vnorm);
            if(!errflag) {
                p->p_r[i]   = (  vout_xyz[0] * cosphi + vout_xyz[1] * sinphi )
                    * gamma * p->mass[i];
                p->p_phi[i] = ( -vout_xyz[0] * sinphi + vout_xyz[1] * cosphi )
                               * gamma * p->mass[i];
                p->p_z[i]   =    vout_xyz[2] * gamma * p->mass[i];
            }

            /* Error handling */
            if(errflag) {
                p->err[i]     = errflag;
                p->running[i] = 0;
            }
        }
    }
}
