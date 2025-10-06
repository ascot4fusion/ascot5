/**
 * Implements B_field.h.
 */
#include "bfield.h"
#include "defines.h"
#include "bfield_analytical.h"
#include "bfield_cartesian.h"
#include "bfield_spline2d.h"
#include "bfield_spline3d.h"
#include "bfield_stellarator.h"
#include "errors.h"
#include <stdio.h>

void B_field_free(Bfield *bfield)
{
    switch (bfield->type)
    {
    case B_field_type_analytical:
        BfieldAnalytical_free(bfield->analytical);
        break;

    case B_field_type_spline2d:
        BfieldSpline2D_free(bfield->spline2d);
        break;

    case B_field_type_spline3d:
        BfieldSpline3D_free(bfield->spline3d);
        break;

    case B_field_type_stellarator:
        BfieldStellarator_free(bfield->stellarator);
        break;

    case B_field_type_cartesian:
        BfieldCartesian_free(bfield->cartesian);
        break;
    }
}

/**
 * @brief Offload data to the accelerator.
 *
 * @param data pointer to the data struct
 */
void B_field_offload(Bfield *bfield)
{
    switch (bfield->type)
    {
    case B_field_type_analytical:
        BfieldAnalytical_offload(bfield->analytical);
        break;

    case B_field_type_spline2d:
        BfieldSpline2D_offload(bfield->spline2d);
        break;

    case B_field_type_spline3d:
        BfieldSpline3D_offload(bfield->spline3d);
        break;

    case B_field_type_stellarator:
        BfieldStellarator_offload(bfield->stellarator);
        break;

    case B_field_type_cartesian:
        BfieldCartesian_offload(bfield->cartesian);
        break;
    }
}

a5err B_field_eval_psi(
    real psi[1], real r, real phi, real z, real t, Bfield *bfield)
{
    a5err err = 0;
    (void)t; /* Unused until dynamic bfield is implemented. */

    switch (bfield->type)
    {
    case B_field_type_analytical:
        err = BfieldAnalytical_eval_psi(psi, r, z, bfield->analytical);
        break;

    case B_field_type_spline2d:
        err = BfieldSpline2D_eval_psi(psi, r, z, bfield->spline2d);
        break;

    case B_field_type_spline3d:
        err = BfieldSpline3D_eval_psi(psi, r, z, bfield->spline3d);
        break;

    case B_field_type_stellarator:
        err = BfieldStellarator_eval_psi(psi, r, phi, z, bfield->stellarator);
        break;

    case B_field_type_cartesian:
        err = BfieldCartesian_eval_psi(psi, bfield->cartesian);
        break;

    default:
        /* Unregonized input. Produce error. */
        err = error_raise(ERR_UNKNOWN_INPUT, __LINE__, EF_B_FIELD);
        break;
    }

    if (err)
    {
        /* In case of error, return some reasonable value to avoid further
           complications */
        psi[0] = 1;
    }

    return err;
}

a5err B_field_eval_psi_dpsi(
    real psi_dpsi[4], real r, real phi, real z, real t, Bfield *bfield)
{
    a5err err = 0;
    (void)t; /* Unused until dynamic bfield is implemented. */

    switch (bfield->type)
    {
    case B_field_type_analytical:
        err =
            BfieldAnalytical_eval_psi_dpsi(psi_dpsi, r, z, bfield->analytical);
        break;

    case B_field_type_spline2d:
        err = BfieldSpline2D_eval_psi_dpsi(psi_dpsi, r, z, bfield->spline2d);
        break;

    case B_field_type_spline3d:
        err = BfieldSpline3D_eval_psi_dpsi(psi_dpsi, r, z, bfield->spline3d);
        break;

    case B_field_type_stellarator:
        err = BfieldStellarator_eval_psi_dpsi(
            psi_dpsi, r, phi, z, bfield->stellarator);
        break;

    case B_field_type_cartesian:
        err = BfieldCartesian_eval_psi_dpsi(psi_dpsi, bfield->cartesian);
        break;

    default:
        /* Unregonized input. Produce error. */
        err = error_raise(ERR_UNKNOWN_INPUT, __LINE__, EF_B_FIELD);
        break;
    }

    if (err)
    {
        /* In case of error, return some reasonable values to avoid further
           complications */
        psi_dpsi[0] = 1;
        for (int k = 1; k < 4; k++)
        {
            psi_dpsi[k] = 0;
        }
    }

    return err;
}

a5err B_field_eval_rho(real rho[2], real psi, Bfield *bfield)
{
    a5err err = 0;

    real psi0 = 0.0, psi1 = 1.0;
    switch (bfield->type)
    {
    case B_field_type_analytical:
        psi0 = bfield->analytical->psilimits[0];
        psi1 = bfield->analytical->psilimits[1];
        break;

    case B_field_type_spline2d:
        psi0 = bfield->spline2d->psilimits[0];
        psi1 = bfield->spline2d->psilimits[1];
        break;

    case B_field_type_spline3d:
        psi0 = bfield->spline3d->psilimits[0];
        psi1 = bfield->spline3d->psilimits[1];
        break;

    case B_field_type_stellarator:
        psi0 = bfield->stellarator->psilimits[0];
        psi1 = bfield->stellarator->psilimits[1];
        break;

    case B_field_type_cartesian:
        psi0 = bfield->cartesian->psival -
               bfield->cartesian->rhoval * bfield->cartesian->rhoval;
        psi1 = bfield->cartesian->psival -
               bfield->cartesian->rhoval * bfield->cartesian->rhoval + 1.0;
        break;

    default:
        /* Unregonized input. Produce error. */
        err = error_raise(ERR_UNKNOWN_INPUT, __LINE__, EF_B_FIELD);
        break;
    }

    real delta = (psi1 - psi0);
    if ((psi - psi0) / delta < 0)
    {
        err = error_raise(ERR_INPUT_UNPHYSICAL, __LINE__, EF_B_FIELD);
    }
    else
    {
        rho[0] = sqrt((psi - psi0) / delta);
        rho[1] = 1.0 / (2 * delta * rho[0]);
    }

    if (err)
    {
        /* In case of error, return some reasonable value to avoid further
           complications */
        rho[0] = 1;
        rho[1] = 0;
    }

    return err;
}

a5err B_field_eval_rho_drho(
    real rho_drho[4], real r, real phi, real z, Bfield *bfield)
{
    a5err err = 0;

    switch (bfield->type)
    {
    case B_field_type_analytical:
        err =
            BfieldAnalytical_eval_rho_drho(rho_drho, r, z, bfield->analytical);
        break;

    case B_field_type_spline2d:
        err = BfieldSpline2D_eval_rho_drho(rho_drho, r, z, bfield->spline2d);
        break;

    case B_field_type_spline3d:
        err = BfieldSpline3D_eval_rho_drho(rho_drho, r, z, bfield->spline3d);
        break;

    case B_field_type_stellarator:
        err = BfieldStellarator_eval_rho_drho(
            rho_drho, r, phi, z, bfield->stellarator);
        break;

    case B_field_type_cartesian:
        err = BfieldCartesian_eval_rho_drho(rho_drho, bfield->cartesian);
        break;

    default:
        /* Unregonized input. Produce error. */
        err = error_raise(ERR_UNKNOWN_INPUT, __LINE__, EF_B_FIELD);
        break;
    }

    if (err)
    {
        /* In case of error, return some reasonable values to avoid further
           complications */
        rho_drho[0] = 1;
        for (int k = 1; k < 4; k++)
        {
            rho_drho[k] = 0;
        }
    }

    return err;
}

a5err B_field_eval_B(
    real b[3], real r, real phi, real z, real t, Bfield *bfield)
{
    a5err err = 0;
    (void)t; /* Unused until dynamic bfield is implemented. */

    switch (bfield->type)
    {
    case B_field_type_analytical:
        err = BfieldAnalytical_eval_b(b, r, phi, z, bfield->analytical);
        break;

    case B_field_type_spline2d:
        err = BfieldSpline2D_eval_b(b, r, z, bfield->spline2d);
        break;

    case B_field_type_spline3d:
        err = BfieldSpline3D_eval_b(b, r, phi, z, bfield->spline3d);
        break;

    case B_field_type_stellarator:
        err = BfieldStellarator_eval_b(b, r, phi, z, bfield->stellarator);
        break;

    case B_field_type_cartesian:
        err = BfieldCartesian_eval_b(b, r, phi, z, bfield->cartesian);
        break;

    default:
        /* Unregonized input. Produce error. */
        err = error_raise(ERR_UNKNOWN_INPUT, __LINE__, EF_B_FIELD);
        break;
    }

    if (err)
    {
        /* In case of error, return some reasonable values to avoid further
           complications */
        b[0] = 1;
        for (int k = 1; k < 3; k++)
        {
            b[k] = 0;
        }
    }

    return err;
}

a5err B_field_eval_B_dB(
    real b_db[15], real r, real phi, real z, real t, Bfield *bfield)
{
    a5err err = 0;
    (void)t; /* Unused until dynamic bfield is implemented. */

    switch (bfield->type)
    {
    case B_field_type_analytical:
        err = BfieldAnalytical_eval_b_db(b_db, r, phi, z, bfield->analytical);
        break;

    case B_field_type_spline2d:
        err = BfieldSpline2D_eval_b_db(b_db, r, z, bfield->spline2d);
        break;

    case B_field_type_spline3d:
        err = BfieldSpline3D_eval_b_db(b_db, r, phi, z, bfield->spline3d);
        break;

    case B_field_type_stellarator:
        err = BfieldStellarator_eval_b_db(b_db, r, phi, z, bfield->stellarator);
        break;

    case B_field_type_cartesian:
        err = BfieldCartesian_eval_b_db(b_db, r, phi, z, bfield->cartesian);
        break;

    default:
        /* Unregonized input. Produce error. */
        err = error_raise(ERR_UNKNOWN_INPUT, __LINE__, EF_B_FIELD);
        break;
    }

    if (err)
    {
        /* In case of error, return some reasonable values to avoid further
           complications */
        b_db[0] = 1;
        for (int k = 1; k < 12; k++)
        {
            b_db[k] = 0;
        }
    }

    return err;
}

a5err B_field_get_axis_rz(real axisrz[2], Bfield *bfield, real phi)
{
    a5err err = 0;

    switch (bfield->type)
    {
    case B_field_type_analytical:
        err = BfieldAnalytical_eval_axisrz(axisrz, bfield->analytical);
        break;

    case B_field_type_spline2d:
        err = BfieldSpline2D_eval_axisrz(axisrz, bfield->spline2d);
        break;

    case B_field_type_spline3d:
        err = BfieldSpline3D_eval_axisrz(axisrz, bfield->spline3d);
        break;

    case B_field_type_stellarator:
        err = BfieldStellarator_eval_axisrz(axisrz, phi, bfield->stellarator);
        break;

    case B_field_type_cartesian:
        err = BfieldCartesian_eval_axisrz(axisrz, bfield->cartesian);
        break;

    default:
        /* Unregonized input. Produce error. */
        err = error_raise(ERR_UNKNOWN_INPUT, __LINE__, EF_B_FIELD);
        break;
    }

    if (err)
    {
        /* In case of error, return some reasonable values to avoid further
           complications */
        axisrz[0] = 1.0;
        axisrz[1] = 0.0;
    }

    return err;
}
