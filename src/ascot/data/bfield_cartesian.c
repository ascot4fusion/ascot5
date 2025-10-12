/**
 * Implements bfield_cartesian.h.
 */
#include "bfield_cartesian.h"
#include "bfield.h"
#include "defines.h"
#include "utils/mathlib.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int BfieldCartesian_init(
    BfieldCartesian *bfield, real psival, real rhoval, real axisrz[2],
    real bxyz[3], real jacobian[9])
{

    bfield->psival = psival;
    bfield->rhoval = rhoval;
    bfield->bxyz[0] = bxyz[0];
    bfield->bxyz[1] = bxyz[1];
    bfield->bxyz[2] = bxyz[2];
    bfield->axisrz[0] = axisrz[0];
    bfield->axisrz[1] = axisrz[1];
    for (size_t i = 0; i < 9; i++)
    {
        bfield->jacobian[i] = jacobian[i];
    }
    return 0;
}

void BfieldCartesian_free(BfieldCartesian *bfield)
{
    (void)bfield;
    // No resources were dynamically allocated
}

void BfieldCartesian_offload(BfieldCartesian *bfield)
{
    (void)bfield;
    GPU_MAP_TO_DEVICE(data->bxyz [0:3], data->jacobian [0:9])
}

err_t BfieldCartesian_eval_psi(real psi[1], BfieldCartesian *bfield)
{
    psi[0] = bfield->psival;
    return 0;
}

err_t BfieldCartesian_eval_psi_dpsi(real psi_dpsi[4], BfieldCartesian *bfield)
{
    psi_dpsi[0] = bfield->psival;
    psi_dpsi[1] = 0;
    psi_dpsi[2] = 0;
    psi_dpsi[3] = 0;

    return 0;
}

err_t BfieldCartesian_eval_b(
    real b[3], real r, real phi, real z, BfieldCartesian *bfield)
{
    real rpz[3] = {r, phi, z}, xyz[3];
    math_rpz2xyz(rpz, xyz);

    real bxyz[3], *jacobian = bfield->jacobian;
    bxyz[0] = jacobian[0] + jacobian[0] * xyz[0] + jacobian[1] * xyz[1] +
              jacobian[2] * xyz[2];
    bxyz[1] = jacobian[1] + jacobian[3] * xyz[0] + jacobian[4] * xyz[1] +
              jacobian[5] * xyz[2];
    bxyz[2] = jacobian[2] + jacobian[6] * xyz[0] + jacobian[7] * xyz[1] +
              jacobian[8] * xyz[2];
    math_vec_xyz2rpz(bxyz, b, phi);

    return 0;
}

err_t BfieldCartesian_eval_b_db(
    real b_db[12], real r, real phi, real z, BfieldCartesian *bfield)
{
    real rpz[3] = {r, phi, z}, xyz[3];
    math_rpz2xyz(rpz, xyz);

    real bxyz[3], *jacobian = bfield->jacobian;
    bxyz[0] = jacobian[0] + jacobian[0] * xyz[0] + jacobian[1] * xyz[1] +
              jacobian[2] * xyz[2];
    bxyz[1] = jacobian[1] + jacobian[3] * xyz[0] + jacobian[4] * xyz[1] +
              jacobian[5] * xyz[2];
    bxyz[2] = jacobian[2] + jacobian[6] * xyz[0] + jacobian[7] * xyz[1] +
              jacobian[8] * xyz[2];
    real b[3];
    math_vec_xyz2rpz(bxyz, b, phi);

    real bxyz_db[12] = {jacobian[0], jacobian[1], jacobian[2],
                        jacobian[3], jacobian[4], jacobian[5],
                        jacobian[6], jacobian[7], jacobian[8]};
    math_jac_xyz2rpz(bxyz_db, b_db, r, phi);
    b_db[0] = b[0];
    b_db[4] = b[1];
    b_db[8] = b[2];

    return 0;
}

err_t BfieldCartesian_eval_axisrz(real axisrz[2], BfieldCartesian *bfield)
{
    axisrz[0] = bfield->axisrz[0];
    axisrz[1] = bfield->axisrz[1];
    return 0;
}
