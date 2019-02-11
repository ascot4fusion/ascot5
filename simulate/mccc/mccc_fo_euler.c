/**
 * @file mccc_fo_euler.c
 * @brief Euler-Maruyama integrator for collision operator in FO picture.
 */
#include <math.h>
#include "../../ascot5.h"
#include "../../consts.h"
#include "../../math.h"
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
 * @param p fo struct
 * @param h time-steps for NSIMD markers
 * @param Bdata pointer to magnetic field
 * @param pdata pointer to plasma data
 * @param rdata pointer to random-generator data
 * @param coldata pointer collision coefficient data for interpolation or NULL
 *        if coefficients are evaluated exactly
 */
void mccc_fo_euler(particle_simd_fo* p, real* h, B_field_data* Bdata,
                   plasma_data* pdata, random_data* rdata, mccc_data* mdata) {

    /* Generate random numbers and get plasma information before going to the *
     * SIMD loop                                                              */
    real rnd[3*NSIMD];
    random_normal_simd(rdata, 3*NSIMD, rnd);

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

            real vin_xyz[3];
            vin_xyz[0] = p->rdot[i] * cosphi - (p->phidot[i]*p->r[i]) * sinphi;
            vin_xyz[1] = p->rdot[i] * sinphi + (p->phidot[i]*p->r[i]) * cosphi;
            vin_xyz[2] = p->zdot[i];
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

            /* Remove energy or pitch change from the results if that is *
             * requested and transform back to cylindrical coordinates.  */
            if(!mdata->include_energy) {
                real vinpervout = vin / math_norm(vout_xyz);
                vout_xyz[0] *= vinpervout;
                vout_xyz[1] *= vinpervout;
                vout_xyz[2] *= vinpervout;
            }
            if(!mdata->include_pitch) {
                real vinpervout = vin / math_norm(vout_xyz);
                vout_xyz[0] = vin_xyz[0] / vinpervout;
                vout_xyz[1] = vin_xyz[1] / vinpervout;
                vout_xyz[2] = vin_xyz[2] / vinpervout;
            }

            if(!errflag) {
                p->rdot[i]   =    vout_xyz[0] * cosphi + vout_xyz[1] * sinphi;
                p->phidot[i] = ( -vout_xyz[0] * sinphi + vout_xyz[1] * cosphi )
                               / p->r[i];
                p->zdot[i]   =    vout_xyz[2];
            }

            /* Error handling */
            if(errflag) {
                p->err[i]     = errflag;
                p->running[i] = 0;
            }
        }
    }
}
