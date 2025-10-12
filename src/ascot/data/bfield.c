/**
 * Implements bfield.h.
 */
#include "bfield.h"
#include "bfield_analytical.h"
#include "bfield_cartesian.h"
#include "bfield_spline2d.h"
#include "bfield_spline3d.h"
#include "bfield_stellarator.h"
#include "defines.h"
#include <stdio.h>

void Bfield_free(Bfield *bfield)
{
    switch (bfield->type)
    {
    case BFIELD_CARTESIAN:
        BfieldCartesian_free(bfield->cartesian);
        break;
    case BFIELD_ANALYTICAL:
        BfieldAnalytical_free(bfield->analytical);
        break;
    case BFIELD_SPLINE2D:
        BfieldSpline2D_free(bfield->spline2d);
        break;
    case BFIELD_SPLINE3D:
        BfieldSpline3D_free(bfield->spline3d);
        break;
    case BFIELD_STELLARATOR:
        BfieldStellarator_free(bfield->stellarator);
        break;
    }
}

/**
 * @brief Offload data to the accelerator.
 *
 * @param data pointer to the data struct
 */
void Bfield_offload(Bfield *bfield)
{
    switch (bfield->type)
    {
    case BFIELD_CARTESIAN:
        BfieldCartesian_offload(bfield->cartesian);
        break;
    case BFIELD_ANALYTICAL:
        BfieldAnalytical_offload(bfield->analytical);
        break;
    case BFIELD_SPLINE2D:
        BfieldSpline2D_offload(bfield->spline2d);
        break;
    case BFIELD_SPLINE3D:
        BfieldSpline3D_offload(bfield->spline3d);
        break;
    case BFIELD_STELLARATOR:
        BfieldStellarator_offload(bfield->stellarator);
        break;
    }
}

err_t Bfield_eval_psi(
    real psi[1], real r, real phi, real z, real t, Bfield *bfield)
{
    (void)t; /* Unused until dynamic bfield is implemented. */
    err_t err = 0;
    switch (bfield->type)
    {
    case BFIELD_CARTESIAN:
        err = BfieldCartesian_eval_psi(psi, bfield->cartesian);
        break;
    case BFIELD_ANALYTICAL:
        err = BfieldAnalytical_eval_psi(psi, r, z, bfield->analytical);
        break;
    case BFIELD_SPLINE2D:
        err = BfieldSpline2D_eval_psi(psi, r, z, bfield->spline2d);
        break;
    case BFIELD_SPLINE3D:
        err = BfieldSpline3D_eval_psi(psi, r, z, bfield->spline3d);
        break;
    case BFIELD_STELLARATOR:
        err = BfieldStellarator_eval_psi(psi, r, phi, z, bfield->stellarator);
        break;
    }

    psi[0] = err ? 1.0 : psi[0];

    return err;
}

err_t Bfield_eval_psi_dpsi(
    real psi_dpsi[4], real r, real phi, real z, real t, Bfield *bfield)
{
    (void)t; /* Unused until dynamic bfield is implemented. */
    err_t err = 0;
    switch (bfield->type)
    {
    case BFIELD_CARTESIAN:
        err = BfieldCartesian_eval_psi_dpsi(psi_dpsi, bfield->cartesian);
        break;
    case BFIELD_ANALYTICAL:
        err =
            BfieldAnalytical_eval_psi_dpsi(psi_dpsi, r, z, bfield->analytical);
        break;
    case BFIELD_SPLINE2D:
        err = BfieldSpline2D_eval_psi_dpsi(psi_dpsi, r, z, bfield->spline2d);
        break;
    case BFIELD_SPLINE3D:
        err = BfieldSpline3D_eval_psi_dpsi(psi_dpsi, r, z, bfield->spline3d);
        break;
    case BFIELD_STELLARATOR:
        err = BfieldStellarator_eval_psi_dpsi(
            psi_dpsi, r, phi, z, bfield->stellarator);
        break;
    }

    psi_dpsi[0] = err ? 1.0 : psi_dpsi[0];
    psi_dpsi[1] = err ? 0.0 : psi_dpsi[1];
    psi_dpsi[2] = err ? 0.0 : psi_dpsi[2];
    psi_dpsi[3] = err ? 0.0 : psi_dpsi[3];

    return err;
}

err_t Bfield_eval_rho(real rho[2], real psi, Bfield *bfield)
{
    real psi0 = 0.0, psi1 = 1.0;
    switch (bfield->type)
    {
    case BFIELD_CARTESIAN:
        psi0 = bfield->cartesian->psival -
               bfield->cartesian->rhoval * bfield->cartesian->rhoval;
        psi1 = bfield->cartesian->psival -
               bfield->cartesian->rhoval * bfield->cartesian->rhoval + 1.0;
        break;
    case BFIELD_ANALYTICAL:
        psi0 = bfield->analytical->psilimits[0];
        psi1 = bfield->analytical->psilimits[1];
        break;
    case BFIELD_SPLINE2D:
        psi0 = bfield->spline2d->psilimits[0];
        psi1 = bfield->spline2d->psilimits[1];
        break;
    case BFIELD_SPLINE3D:
        psi0 = bfield->spline3d->psilimits[0];
        psi1 = bfield->spline3d->psilimits[1];
        break;
    case BFIELD_STELLARATOR:
        psi0 = bfield->stellarator->psilimits[0];
        psi1 = bfield->stellarator->psilimits[1];
        break;
    }

    err_t err = 0;
    real rho_squared = (psi - psi0) / (psi1 - psi0);
    err = ERROR_CHECK(err, rho_squared < 0, ERR_IMAGINARY_RHO, DATA_BFIELD_C);
    rho[0] = sqrt(fabs(rho_squared));
    rho[1] = 1.0 / (2 * (psi1 - psi0) * (rho[0] + __DBL_EPSILON__));

    rho[0] = err ? 1.0 : rho[0];
    rho[1] = err ? 0.0 : rho[1];

    return err;
}

err_t Bfield_eval_rho_drho(real rho_drho[4], real psi_dpsi[4], Bfield *bfield)
{
    err_t err = 0;
    real psi0, delta;
    switch (bfield->type)
    {
    case BFIELD_CARTESIAN:
        psi0 = bfield->cartesian->rhoval;
        delta = 1.0;
        break;
    case BFIELD_ANALYTICAL:
        psi0 = bfield->analytical->psilimits[0];
        delta = bfield->analytical->psilimits[1] - psi0;
        break;
    case BFIELD_SPLINE2D:
        psi0 = bfield->spline2d->psilimits[0];
        delta = bfield->spline2d->psilimits[1] - psi0;
        break;
    case BFIELD_SPLINE3D:
        psi0 = bfield->spline3d->psilimits[0];
        delta = bfield->spline3d->psilimits[1] - psi0;
        break;
    case BFIELD_STELLARATOR:
        psi0 = bfield->stellarator->psilimits[0];
        delta = bfield->stellarator->psilimits[1] - psi0;
        break;
    }

    real rho_squared = (psi_dpsi[0] - psi0) / delta;
    err = ERROR_CHECK(err, rho_squared < 0, ERR_IMAGINARY_RHO, DATA_BFIELD_C);
    rho_drho[0] = err ? 1.0 : sqrt(fabs(rho_squared));
    rho_drho[1] = rho_drho[0] ? psi_dpsi[1] / (2 * delta * rho_drho[0]) : 0.0;
    rho_drho[2] = rho_drho[0] ? psi_dpsi[2] / (2 * delta * rho_drho[0]) : 0.0;
    rho_drho[3] = rho_drho[0] ? psi_dpsi[2] / (2 * delta * rho_drho[0]) : 0.0;

    rho_drho[1] = err ? 0.0 : rho_drho[1];
    rho_drho[2] = err ? 0.0 : rho_drho[2];
    rho_drho[3] = err ? 0.0 : rho_drho[3];

    return err;
}

err_t Bfield_eval_b(real b[3], real r, real phi, real z, real t, Bfield *bfield)
{
    (void)t; /* Unused until dynamic bfield is implemented. */

    err_t err = 0;
    switch (bfield->type)
    {
    case BFIELD_CARTESIAN:
        err = BfieldCartesian_eval_b(b, r, phi, z, bfield->cartesian);
        break;
    case BFIELD_ANALYTICAL:
        err = BfieldAnalytical_eval_b(b, r, phi, z, bfield->analytical);
        break;
    case BFIELD_SPLINE2D:
        err = BfieldSpline2D_eval_b(b, r, z, bfield->spline2d);
        break;
    case BFIELD_SPLINE3D:
        err = BfieldSpline3D_eval_b(b, r, phi, z, bfield->spline3d);
        break;
    case BFIELD_STELLARATOR:
        err = BfieldStellarator_eval_b(b, r, phi, z, bfield->stellarator);
        break;
    }

    b[0] = err ? 1.0 : b[0];
    b[1] = err ? 1.0 : b[1];
    b[2] = err ? 1.0 : b[2];

    return err;
}

err_t Bfield_eval_b_db(
    real b_db[15], real r, real phi, real z, real t, Bfield *bfield)
{
    (void)t; /* Unused until dynamic bfield is implemented. */

    err_t err = 0;
    switch (bfield->type)
    {
    case BFIELD_CARTESIAN:
        err = BfieldCartesian_eval_b_db(b_db, r, phi, z, bfield->cartesian);
        break;
    case BFIELD_ANALYTICAL:
        err = BfieldAnalytical_eval_b_db(b_db, r, phi, z, bfield->analytical);
        break;
    case BFIELD_SPLINE2D:
        err = BfieldSpline2D_eval_b_db(b_db, r, z, bfield->spline2d);
        break;
    case BFIELD_SPLINE3D:
        err = BfieldSpline3D_eval_b_db(b_db, r, phi, z, bfield->spline3d);
        break;
    case BFIELD_STELLARATOR:
        err = BfieldStellarator_eval_b_db(b_db, r, phi, z, bfield->stellarator);
        break;
    }

    b_db[0] = err ? 1.0 : b_db[0];
    b_db[1] = err ? 0.0 : b_db[1];
    b_db[2] = err ? 0.0 : b_db[2];
    b_db[3] = err ? 0.0 : b_db[3];
    b_db[4] = err ? 1.0 : b_db[4];
    b_db[5] = err ? 0.0 : b_db[5];
    b_db[6] = err ? 0.0 : b_db[6];
    b_db[7] = err ? 0.0 : b_db[7];
    b_db[8] = err ? 1.0 : b_db[8];
    b_db[9] = err ? 0.0 : b_db[9];
    b_db[10] = err ? 0.0 : b_db[10];
    b_db[11] = err ? 0.0 : b_db[11];

    return err;
}

err_t Bfield_eval_axis_rz(real axisrz[2], Bfield *bfield, real phi)
{
    err_t err = 0;
    switch (bfield->type)
    {
    case BFIELD_CARTESIAN:
        err = BfieldCartesian_eval_axisrz(axisrz, bfield->cartesian);
        break;
    case BFIELD_ANALYTICAL:
        err = BfieldAnalytical_eval_axisrz(axisrz, bfield->analytical);
        break;
    case BFIELD_SPLINE2D:
        err = BfieldSpline2D_eval_axisrz(axisrz, bfield->spline2d);
        break;
    case BFIELD_SPLINE3D:
        err = BfieldSpline3D_eval_axisrz(axisrz, bfield->spline3d);
        break;
    case BFIELD_STELLARATOR:
        err = BfieldStellarator_eval_axisrz(axisrz, phi, bfield->stellarator);
        break;
    }

    axisrz[0] = err ? 1.0 : axisrz[0];
    axisrz[1] = err ? 0.0 : axisrz[0];

    return err;
}
