/**
 * Implements boozer.h.
 */
#include "boozer.h"
#include "bfield.h"
#include "consts.h"
#include "defines.h"
#include "utils/interp.h"
#include "utils/mathlib.h"
#include <math.h>
#include <stdlib.h>

int boozer_init(
    Boozer *boozer, size_t npsi, real psi_min, real psi_max, size_t ntheta,
    size_t nthetag, int npadding, real *nu, real *theta, size_t nrzs, real *rs,
    real *zs)
{

    int err = 0;
    real THETAMIN = 0;
    real THETAMAX = CONST_2PI;
    real padding = (npadding * CONST_2PI) / (nthetag - 2 * npadding - 1);
    boozer->psilim[0] = psi_min;
    boozer->psilim[1] = psi_max;

    err = interp2Dcomp_setup(
        &boozer->nu, nu, npsi, ntheta, NATURALBC, PERIODICBC, psi_min, psi_max,
        THETAMIN, THETAMAX);
    err = interp2Dcomp_setup(
        &boozer->theta, theta, npsi, nthetag, NATURALBC, NATURALBC, psi_min,
        psi_max, THETAMIN - padding, THETAMAX + padding);

    boozer->nrz = nrzs;
    boozer->rlim = (real *)malloc(nrzs * sizeof(real));
    boozer->zlim = (real *)malloc(nrzs * sizeof(real));
    for (size_t i = 0; i < nrzs; i++)
    {
        boozer->rlim[i] = rs[i];
        boozer->zlim[i] = zs[i];
    }
    return err;
}

void boozer_free(Boozer *boozer)
{
    free(boozer->rlim);
    free(boozer->zlim);
    free(boozer->nu.c);
    free(boozer->theta.c);
}

void boozer_offload(Boozer *boozer) { SUPPRESS_UNUSED_WARNING(boozer); }

err_t boozer_eval_psithetazeta(
    real psithetazeta[12], int *isinside, real r, real phi, real z,
    Bfield *bfield, Boozer *boozer)
{
    err_t err = 0;
    int interperr = 0;

    /* Winding number to test whether we are inside the plasma (and not in the
       private plasma region) */
    isinside[0] = 0;
    if (math_point_in_polygon(r, z, boozer->rlim, boozer->zlim, boozer->nrz))
    {
        /* Get the psi value and check that it is within the psi grid (the grid
           does not extend all the way to the axis) Use t = 0.0 s */
        real psi[4], rho[2], r0, z0;
        err = Bfield_eval_psi_dpsi(psi, r, phi, z, 0.0, bfield);
        if (!err)
        {
            err = Bfield_eval_rho(rho, psi[0], bfield);
        }
        if (!err)
        {
            real rz[2];
            err = Bfield_eval_axis_rz(rz, bfield, phi);
            r0 = rz[0];
            z0 = rz[1];
        }

        if (!err && psi[0] < boozer->psilim[1] && psi[0] > boozer->psilim[0])
        {

            /* Update the flag, and we are good to go */
            isinside[0] = 1;

            /* Geometrical theta */
            real thgeo;
            thgeo = fmod(atan2(z - z0, r - r0) + CONST_2PI, CONST_2PI);

            /* Boozer theta and derivatives */
            real theta[6];
            interperr +=
                interp2Dcomp_eval_df(theta, &boozer->theta, psi[0], thgeo);

            /* Boozer nu function and derivatives */
            real nu[6];
            interperr +=
                interp2Dcomp_eval_df(nu, &boozer->nu, psi[0], theta[0]);

            /* Set up data for returning the requested values */

            /* Psi and derivatives */
            psithetazeta[0] = psi[0]; /* psi       */
            psithetazeta[1] = psi[1]; /* dpsi_dr   */
            psithetazeta[2] = 0;      /* dpsi_dphi */
            psithetazeta[3] = psi[3]; /* dpsi_dz   */

            /* Helpers */
            real asq;
            asq = (r - r0) * (r - r0) + (z - z0) * (z - z0);
            real dthgeo_dr;
            dthgeo_dr = -(z - z0) / asq;
            real dthgeo_dz;
            dthgeo_dz = (r - r0) / asq;

            /* Theta and derivatives */
            psithetazeta[4] = theta[0]; /* theta       */
            psithetazeta[5] =
                theta[1] * psi[1] + theta[2] * dthgeo_dr; /* dtheta_dr   */
            psithetazeta[6] = 0;                          /* dtheta_dphi */
            psithetazeta[7] =
                theta[1] * psi[3] + theta[2] * dthgeo_dz; /* dtheta_dz   */

            /* Zeta and derivatives */
            psithetazeta[8] = fmod(phi + nu[0], CONST_2PI); /* zeta      */
            psithetazeta[9] =
                nu[1] * psi[1] + nu[2] * psithetazeta[5]; /* dzeta_dR  */
            psithetazeta[10] = 1.0;                       /* dzeta_dphi*/
            psithetazeta[11] =
                nu[1] * psi[3] + nu[2] * psithetazeta[7]; /* dzeta_dz  */

            /* Make sure zeta is between [0, 2pi]*/
            psithetazeta[8] = fmod(psithetazeta[8] + CONST_2PI, CONST_2PI);
        }
    }

    err = ERROR_CHECK(
        err, interperr, ERR_INTERPOLATED_OUTSIDE_RANGE, DATA_BOOZER_C);

    return err;
}
