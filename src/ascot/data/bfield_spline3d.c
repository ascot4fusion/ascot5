/**
 * Implements bfield_spline3d.h.
 */
#include "bfield_spline3d.h"
#include "bfield.h"
#include "defines.h"
#include "utils/interp.h"
#include "utils/mathlib.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int BfieldSpline3D_init(
    BfieldSpline3D *bfield, size_t pnr, size_t pnz, size_t bnr, size_t bnz,
    size_t bnphi, real prlim[2], real pzlim[2], real brlim[2], real bzlim[2],
    real bphilim[2], real axisrz[2], real psilimits[2], real psi[pnr * pnz],
    real br[bnr * bnz * bnphi], real bz[bnr * bnz * bnphi],
    real bphi[bnr * bnz * bnphi])
{
    int err = 0;
    bfield->axisrz[0] = axisrz[0];
    bfield->axisrz[1] = axisrz[1];
    bfield->psilimits[0] = psilimits[0];
    bfield->psilimits[1] = psilimits[1];

    err += Spline2D_init(
        &bfield->psi, pnr, pnz, NATURALBC, NATURALBC, prlim,
        pzlim, psi);
    err += Spline3D_init(
        &bfield->br, bnr, bnphi, bnz, NATURALBC, PERIODICBC, NATURALBC,
        brlim, bphilim, bzlim, br);
    err += Spline3D_init(
        &bfield->bz, bnr, bnphi, bnz, NATURALBC, PERIODICBC, NATURALBC,
        brlim, bphilim, bzlim, bz);
    err += Spline3D_init(
        &bfield->bphi, bnr, bnphi, bnz, NATURALBC, PERIODICBC, NATURALBC,
        brlim, bphilim, bzlim, bphi);

    return err;
}

void BfieldSpline3D_free(BfieldSpline3D *bfield)
{
    free(bfield->psi.c);
    free(bfield->br.c);
    free(bfield->bphi.c);
    free(bfield->bz.c);
}

void BfieldSpline3D_offload(BfieldSpline3D *bfield)
{
    (void)bfield;
    GPU_MAP_TO_DEVICE(
        bfield->psi, bfield->br, bfield->bphi, bfield->bz,
        bfield->psi.c [0:bfield->psi.n_x * bfield->psi.n_y * NSIZE_COMP2D],
        bfield->br.c
        [0:bfield->br.n_x * bfield->br.n_y * bfield->br.n_z * NSIZE_COMP3D],
        bfield->bphi.c [0:bfield->bphi.n_x * bfield->bphi.n_y *
                          bfield->bphi.n_z * NSIZE_COMP3D],
        bfield->bz.c
        [0:bfield->bz.n_x * bfield->bz.n_y * bfield->bz.n_z * NSIZE_COMP3D])
}

err_t BfieldSpline3D_eval_psi(
    real psi[1], real r, real z, BfieldSpline3D *bfield)
{
    err_t err = 0;
    int interperr = Spline2D_eval_f(&psi[0], &bfield->psi, r, z);

    err = ERROR_CHECK(
        err, interperr, ERR_INTERPOLATED_OUTSIDE_RANGE, DATA_BFIELD_SPLINE3D_C);
    return err;
}

err_t BfieldSpline3D_eval_psi_dpsi(
    real psi_dpsi[4], real r, real z, BfieldSpline3D *bfield)
{
    err_t err = 0;
    real psi_dpsi_temp[6];
    int interperr = Spline2D_eval_f_df(psi_dpsi_temp, &bfield->psi, r, z);
    psi_dpsi[0] = psi_dpsi_temp[0];
    psi_dpsi[1] = psi_dpsi_temp[1];
    psi_dpsi[2] = 0;
    psi_dpsi[3] = psi_dpsi_temp[2];

    err = ERROR_CHECK(
        err, interperr, ERR_INTERPOLATED_OUTSIDE_RANGE, DATA_BFIELD_SPLINE3D_C);
    return err;
}

err_t BfieldSpline3D_eval_b(
    real b[3], real r, real phi, real z, BfieldSpline3D *bfield)
{
    err_t err = 0;
    int interperr = 0;

    interperr += Spline3D_eval_f(&b[0], &bfield->br, r, phi, z);
    interperr += Spline3D_eval_f(&b[1], &bfield->bphi, r, phi, z);
    interperr += Spline3D_eval_f(&b[2], &bfield->bz, r, phi, z);

    err = ERROR_CHECK(
        err, interperr, ERR_INTERPOLATED_OUTSIDE_RANGE, DATA_BFIELD_SPLINE3D_C);

    real psi_dpsi[6];
    interperr = Spline2D_eval_f_df(psi_dpsi, &bfield->psi, r, z);
    b[0] -= psi_dpsi[2] / r;
    b[2] += psi_dpsi[1] / r;

    err = ERROR_CHECK(
        err, interperr, ERR_INTERPOLATED_OUTSIDE_RANGE, DATA_BFIELD_SPLINE3D_C);
    return err;
}

err_t BfieldSpline3D_eval_b_db(
    real b_db[12], real r, real phi, real z, BfieldSpline3D *bfield)
{
    err_t err = 0;
    int interperr = 0;
    real b_db_temp[10];

    interperr += Spline3D_eval_f_df(b_db_temp, &bfield->br, r, phi, z);
    b_db[0] = b_db_temp[0];
    b_db[2] = b_db_temp[1];
    b_db[3] = b_db_temp[2];
    b_db[4] = b_db_temp[3];

    interperr += Spline3D_eval_f_df(b_db_temp, &bfield->bphi, r, phi, z);
    b_db[1] = b_db_temp[0];
    b_db[6] = b_db_temp[1];
    b_db[7] = b_db_temp[2];
    b_db[8] = b_db_temp[3];

    interperr += Spline3D_eval_f_df(b_db_temp, &bfield->bz, r, phi, z);
    b_db[2] = b_db_temp[0];
    b_db[9] = b_db_temp[1];
    b_db[10] = b_db_temp[2];
    b_db[11] = b_db_temp[3];

    err = ERROR_CHECK(
        err, interperr, ERR_INTERPOLATED_OUTSIDE_RANGE, DATA_BFIELD_SPLINE3D_C);

    real psi_dpsi[6];
    interperr = Spline2D_eval_f_df(psi_dpsi, &bfield->psi, r, z);

    b_db[0] -= psi_dpsi[2] / r;
    b_db[3] += psi_dpsi[2] / (r * r) - psi_dpsi[5] / r;
    b_db[5] -= psi_dpsi[4] / r;
    b_db[2] += psi_dpsi[1] / r;
    b_db[9] -= psi_dpsi[1] / (r * r) + psi_dpsi[3] / r;
    b_db[11] += psi_dpsi[5] / r;

    err = ERROR_CHECK(
        err, interperr, ERR_INTERPOLATED_OUTSIDE_RANGE, DATA_BFIELD_SPLINE3D_C);
    return err;
}

err_t BfieldSpline3D_eval_axisrz(real axisrz[2], BfieldSpline3D *bfield)
{
    axisrz[0] = bfield->axisrz[0];
    axisrz[1] = bfield->axisrz[1];
    return 0;
}
