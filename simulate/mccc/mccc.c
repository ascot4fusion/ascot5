/**
 * @file mccc.c
 * @brief Interface for using mccc package within ascot5
 */
#include <stdlib.h>
#include <math.h>
#include "../../ascot5.h"
#include "../../B_field.h"
#include "../../plasma.h"
#include "mccc.h"
#include "mccc_coefs.h"

/**
 * @brief Set collision operator data.
 *
 * @param mdata pointer to collision operator data struct
 * @param include_energy can collisions change marker energy, either 0 or 1
 * @param include_pitch  can collisions change marker pitch, either 0 or 1
 * @param include_gcdiff can collisions change GC position, either 0 or 1
 */
void mccc_init(mccc_data* mdata, int include_energy, int include_pitch,
               int include_gcdiff) {
    mdata->include_energy = include_energy;
    mdata->include_pitch  = include_pitch;
    mdata->include_gcdiff = include_gcdiff;

    mdata->usetabulated = 0;
}

/**
 * @brief Evaluate collision coefficients
 *
 * This function is not called during the simulation but is used as a way
 * to get easy access to collision coefficients. Coefficients are evaluated
 * for plasma parameters found on given coordinates, for an array of given
 * velocities.
 *
 * Evaluated coefficients are stored in the given arrays as
 * [nv*i_species + i_v]. If a NULL pointer is given, then that coefficient is
 * not evaluated.
 *
 * @param ma particle mass [kg]
 * @param qa particle charge [C]
 * @param r R-coordinate where plasma is evaluated [m]
 * @param phi phi-coordinate where plasma is evaluated [rad]
 * @param z z-coordinate where plasma is evaluated [m]
 * @param t time-coordinate where plasma is evaluated [s]
 * @param va array of velocities [m/s]
 * @param nv number of velocities
 * @param pdata pointer to plasma data
 * @param Bdata pointer to magnetic field data
 * @param F
 * @param Dpara
 * @param Dperp
 * @param K
 * @param nu
 *
 * @return zero if plasma parameters at (r, phi, z, t) evaluated succesfully
 */
int mccc_eval_coefs(real ma, real qa, real r, real phi, real z, real t,
                    real* va, int nv, plasma_data* pdata,  B_field_data* Bdata,
                    real* F, real* Dpara, real* Dperp, real* K, real* nu) {

    mccc_data mdata;
    mdata.usetabulated = 0;

    /* Evaluate rho as it is needed to evaluate plasma parameters */
    real psi, rho;
    if( B_field_eval_psi(&psi, r, phi, z, t, Bdata) ) {
        return 1;
    }
    if( B_field_eval_rho(&rho, psi, Bdata) ) {
        return 1;
    }

    /* Evaluate plasma parameters */
    int n_species  = plasma_get_n_species(pdata);
    const real* qb = plasma_get_species_charge(pdata);
    const real* mb = plasma_get_species_mass(pdata);

    real nb[MAX_SPECIES], Tb[MAX_SPECIES];
    if( plasma_eval_densandtemp(nb, Tb, rho, r, phi, z, t, pdata) ) {
        return 1;
    }

    /* Loop through all plasma species */
    for(int ib=0; ib<n_species; ib++) {

        /* Evaluate coefficients for different velocities */
        for(int iv=0; iv<nv; iv++) {

            /* Coulomb logarithm */
            real clogab[MAX_SPECIES];
            mccc_coefs_clog(clogab, ma, qa, va[iv], n_species, mb, qb, nb, Tb);

            /* Special functions */
            real vb = sqrt( 2 * Tb[ib] / mb[ib] );
            real x  = va[iv] / vb;
            real mufun[3];
            mccc_coefs_mufun(mufun, x, &mdata);

            /* Coefficients */
            real Fb      = mccc_coefs_F(ma, qa, mb[ib], qb[ib], nb[ib], vb,
                                        clogab[ib], mufun[0]);
            real Qb      = mccc_coefs_Q(ma, qa, mb[ib], qb[ib], nb[ib], vb,
                                        clogab[ib], mufun[0]);
            real Dparab  = mccc_coefs_Dpara(ma, qa, va[iv], qb[ib], nb[ib], vb,
                                            clogab[ib], mufun[0]);
            real Dperpb  = mccc_coefs_Dperp(ma, qa, va[iv], qb[ib], nb[ib], vb,
                                            clogab[ib], mufun[1]);
            real dDparab = mccc_coefs_dDpara(ma, qa, va[iv], qb[ib], nb[ib], vb,
                                             clogab[ib], mufun[0], mufun[2]);
            real Kb      = mccc_coefs_K(va[iv], Dparab, dDparab, Qb);
            real nub     = mccc_coefs_nu(va[iv], Dperpb);

            /* Store requested quantities */
            if(F != NULL) {
                F[nv*ib + iv]     = Fb;
            }
            if(Dpara != NULL) {
                Dpara[nv*ib + iv] = Dparab;
            }
            if(Dperp != NULL) {
                Dperp[nv*ib + iv] = Dperpb;
            }
            if(K != NULL) {
                K[nv*ib + iv]     = Kb;
            }
            if(nu != NULL) {
                nu[nv*ib + iv]    = nub;
            }
        }
    }

    return 0;
}
