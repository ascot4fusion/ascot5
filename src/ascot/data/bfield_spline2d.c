/**
 * Implements bfield_spline2d.h.
 */
#include "bfield_spline2d.h"
#include "bfield.h"
#include "defines.h"
#include "utils/interp.h"
#include "utils/mathlib.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int BfieldSpline2D_init(
    BfieldSpline2D *bfield, size_t nr, size_t nz, real rlim[2], real zlim[2],
    real axisrz[2], real psilimits[2], real psi[nr * nz], real br[nr * nz],
    real bz[nr * nz], real bphi[nr * nz])
{
    int err = 0;
    bfield->axisrz[0] = axisrz[0];
    bfield->axisrz[1] = axisrz[1];
    bfield->psilimits[0] = psilimits[0];
    bfield->psilimits[1] = psilimits[1];

    err += Spline2D_init(
        &bfield->psi, nr, nz, NATURALBC, NATURALBC, rlim, zlim, psi);
    err += Spline2D_init(
        &bfield->br, nr, nz, NATURALBC, NATURALBC, rlim, zlim, br);
    err += Spline2D_init(
        &bfield->bz, nr, nz, NATURALBC, NATURALBC, rlim, zlim, bz);
    err += Spline2D_init(
        &bfield->bphi, nr, nz, NATURALBC, NATURALBC, rlim, zlim, bphi);
    return err;
}

void BfieldSpline2D_free(BfieldSpline2D *bfield)
{
    free(bfield->psi.c);
    free(bfield->br.c);
    free(bfield->bphi.c);
    free(bfield->bz.c);
}

void BfieldSpline2D_offload(BfieldSpline2D *bfield)
{
    (void)bfield;
    GPU_MAP_TO_DEVICE(
        bfield->psi, bfield->br, bfield->bphi, bfield->bz,
        bfield->psi.c [0:bfield->psi.nx * bfield->psi.ny * NSIZE_COMP2D],
        bfield->br.c [0:bfield->br.nx * bfield->br.ny * NSIZE_COMP2D],
        bfield->bphi.c [0:bfield->bphi.nx * bfield->bphi.ny * NSIZE_COMP2D],
        bfield->bz.c [0:bfield->bz.nx * bfield->bz.ny * NSIZE_COMP2D])
}

err_t BfieldSpline2D_eval_psi(
    real psi[1], real r, real z, BfieldSpline2D *bfield)
{
    err_t err = 0;
    int interperr = Spline2D_eval_f(&psi[0], &bfield->psi, r, z);

    err = ERROR_CHECK(
        err, interperr, ERR_INTERPOLATED_OUTSIDE_RANGE, DATA_BFIELD_SPLINE2D_C);
    return err;
}

err_t BfieldSpline2D_eval_psi_dpsi(
    real psi_dpsi[4], real r, real z, BfieldSpline2D *bfield)
{
    err_t err = 0;
    real psi_dpsi_temp[6];
    int interperr = Spline2D_eval_f_df(psi_dpsi_temp, &bfield->psi, r, z);
    psi_dpsi[0] = psi_dpsi_temp[0];
    psi_dpsi[1] = psi_dpsi_temp[1];
    psi_dpsi[2] = 0;
    psi_dpsi[3] = psi_dpsi_temp[2];

    err = ERROR_CHECK(
        err, interperr, ERR_INTERPOLATED_OUTSIDE_RANGE, DATA_BFIELD_SPLINE2D_C);
    return err;
}

err_t BfieldSpline2D_eval_b(real b[3], real r, real z, BfieldSpline2D *bfield)
{
    err_t err = 0;
    int interperr = 0;
    interperr += Spline2D_eval_f(&b[0], &bfield->br, r, z);
    interperr += Spline2D_eval_f(&b[1], &bfield->bphi, r, z);
    interperr += Spline2D_eval_f(&b[2], &bfield->bz, r, z);

    err = ERROR_CHECK(
        err, interperr, ERR_INTERPOLATED_OUTSIDE_RANGE, DATA_BFIELD_SPLINE2D_C);

    real psi_dpsi[6];
    interperr = Spline2D_eval_f_df(psi_dpsi, &bfield->psi, r, z);
    b[0] -= psi_dpsi[2] / r;
    b[2] += psi_dpsi[1] / r;

    err = ERROR_CHECK(
        err, interperr, ERR_INTERPOLATED_OUTSIDE_RANGE, DATA_BFIELD_SPLINE2D_C);
    return err;
}

err_t BfieldSpline2D_eval_b_db(
    real b_db[12], real r, real z, BfieldSpline2D *bfield)
{
    err_t err = 0;
    real b_db_temp[6];
    int interperr = 0;

    interperr += Spline2D_eval_f_df(b_db_temp, &bfield->br, r, z);
    b_db[0] = b_db_temp[0];
    b_db[3] = b_db_temp[1];
    b_db[4] = 0;
    b_db[5] = b_db_temp[2];

    interperr += Spline2D_eval_f_df(b_db_temp, &bfield->bphi, r, z);
    b_db[1] = b_db_temp[0];
    b_db[6] = b_db_temp[1];
    b_db[7] = 0;
    b_db[8] = b_db_temp[2];

    interperr += Spline2D_eval_f_df(b_db_temp, &bfield->bz, r, z);
    b_db[2] = b_db_temp[0];
    b_db[9] = b_db_temp[1];
    b_db[10] = 0;
    b_db[11] = b_db_temp[2];

    err = ERROR_CHECK(
        err, interperr, ERR_INTERPOLATED_OUTSIDE_RANGE, DATA_BFIELD_SPLINE2D_C);

    real psi_dpsi[6];
    interperr = Spline2D_eval_f_df(psi_dpsi, &bfield->psi, r, z);
    b_db[0] += -psi_dpsi[2] / r;
    b_db[3] += psi_dpsi[2] / (r * r) - psi_dpsi[5] / r;
    b_db[5] += -psi_dpsi[4] / r;
    b_db[2] += psi_dpsi[1] / r;
    b_db[9] += -psi_dpsi[1] / (r * r) + psi_dpsi[3] / r;
    b_db[11] += psi_dpsi[5] / r;

    err = ERROR_CHECK(
        err, interperr, ERR_INTERPOLATED_OUTSIDE_RANGE, DATA_BFIELD_SPLINE2D_C);
    return err;
}

err_t BfieldSpline2D_eval_axisrz(real axisrz[2], BfieldSpline2D *bfield)
{
    axisrz[0] = bfield->axisrz[0];
    axisrz[1] = bfield->axisrz[1];
    return 0;
}
