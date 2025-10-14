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

int Boozer_init(
    Boozer *boozer, size_t npsi, size_t ntheta, size_t nthetag, size_t nrz,
    int npadding, real psilim[2], real nu[npsi * ntheta],
    real theta[npsi * nthetag], real rs[nrz], real zs[nrz])
{
    int err = 0;
    real THETAMIN = 0;
    real THETAMAX = CONST_2PI;
    real padding = (npadding * CONST_2PI) / (nthetag - 2 * npadding - 1);

    err = interp2Dcomp_setup(
        &boozer->nu, nu, npsi, ntheta, NATURALBC, PERIODICBC, psilim[0],
        psilim[1], THETAMIN, THETAMAX);
    err = interp2Dcomp_setup(
        &boozer->theta, theta, npsi, nthetag, NATURALBC, NATURALBC, psilim[0],
        psilim[1], THETAMIN - padding, THETAMAX + padding);

    boozer->nrz = nrz;
    boozer->rlim = (real *)malloc(nrz * sizeof(real));
    boozer->zlim = (real *)malloc(nrz * sizeof(real));
    for (size_t i = 0; i < nrz; i++)
    {
        boozer->rlim[i] = rs[i];
        boozer->zlim[i] = zs[i];
    }
    return err;
}

void Boozer_free(Boozer *boozer)
{
    free(boozer->rlim);
    free(boozer->zlim);
    free(boozer->nu.c);
    free(boozer->theta.c);
}

void Boozer_offload(Boozer *boozer) { SUPPRESS_UNUSED_WARNING(boozer); }

err_t Boozer_map_coordinates(
    real psithetazeta[12], int *isinside, real r, real phi, real z, real t,
    Boozer *boozer, Bfield *bfield)
{
    err_t err = 0;

    /* These should be defined both inside and outside Boozer region */
    real psi_dpsi[4], rho[2], axisrz[2];
    err = Bfield_eval_psi_dpsi(psi_dpsi, r, phi, z, t, bfield);
    err = err ? err : Bfield_eval_rho(rho, psi_dpsi[0], bfield);
    err = err ? err : Bfield_eval_axis_rz(axisrz, bfield, phi);

    isinside[0] =
        psi_dpsi[0] < boozer->theta.y_max &&
        psi_dpsi[0] > boozer->theta.y_min &&
        math_point_in_polygon(r, z, boozer->rlim, boozer->zlim, boozer->nrz);

    int interperr = 0;
    real thgeo, theta_dtheta[6], nu[6];
    thgeo = fmod(atan2(z - axisrz[1], r - axisrz[0]) + CONST_2PI, CONST_2PI);
    interperr +=
        interp2Dcomp_eval_df(theta_dtheta, &boozer->theta, psi_dpsi[0], thgeo);
    interperr +=
        interp2Dcomp_eval_df(nu, &boozer->nu, psi_dpsi[0], theta_dtheta[0]);

    /* No need to raise error if we are outside the Boozer region */
    err = ERROR_CHECK(
        err, interperr && isinside[0], ERR_INTERPOLATED_OUTSIDE_RANGE,
        DATA_BOOZER_C);

    /* Some helper variables for calucating derivatives */
    real asq, dthgeo_dr, dthgeo_dz;
    asq = (r - axisrz[0]) * (r - axisrz[0]) + (z - axisrz[1]) * (z - axisrz[1]);
    dthgeo_dr = -(z - axisrz[1]) / asq;
    dthgeo_dz = (r - axisrz[0]) / asq;

    real *psi = &psithetazeta[0], *theta = &psithetazeta[4],
         *zeta = &psithetazeta[8];
    psi[0] = psi_dpsi[0];
    psi[1] = psi_dpsi[1];
    psi[2] = psi_dpsi[2];
    psi[3] = psi_dpsi[3];

    theta[0] = theta_dtheta[0];
    theta[1] = theta_dtheta[1] * psi[1] + theta_dtheta[2] * dthgeo_dr;
    theta[2] = 0;
    theta[3] = theta_dtheta[1] * psi[3] + theta_dtheta[2] * dthgeo_dz;

    zeta[0] = fmod(fmod(phi + nu[0], CONST_2PI) + CONST_2PI, CONST_2PI);
    zeta[1] = nu[1] * psi[1] + nu[2] * theta_dtheta[1];
    zeta[2] = 1.0;
    zeta[3] = nu[1] * psi[3] + nu[2] * theta_dtheta[3];

    return err;
}
