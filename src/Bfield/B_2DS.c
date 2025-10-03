/**
 * Implements B_2DS.h.
 */
#include "B_2DS.h"
#include "ascot5.h"
#include "error.h"
#include "interp.h"
#include "math.h"
#include <math.h>
#include <stdlib.h>

int BfieldSpline2D_init(
    BfieldSpline2D *bfield, int nr, int nz, real rlim[2], real zlim[2],
    real axisrz[2], real psilimits[2], real psi[nr * nz], real br[nr * nz],
    real bz[nr * nz], real bphi[nr * nz])
{

    int err = 0;
    bfield->axisrz[0] = axisrz[0];
    bfield->axisrz[1] = axisrz[1];
    bfield->psilimits[0] = psilimits[0];
    bfield->psilimits[1] = psilimits[1];

    err += interp2Dcomp_setup(
        &bfield->psi, psi, nr, nz, NATURALBC, NATURALBC, rlim[0], rlim[1],
        zlim[0], zlim[1]);
    err += interp2Dcomp_setup(
        &bfield->br, br, nr, nz, NATURALBC, NATURALBC, rlim[0], rlim[1],
        zlim[0], zlim[1]);
    err += interp2Dcomp_setup(
        &bfield->bz, bz, nr, nz, NATURALBC, NATURALBC, rlim[0], rlim[1],
        zlim[0], zlim[1]);
    err += interp2Dcomp_setup(
        &bfield->bphi, bphi, nr, nz, NATURALBC, NATURALBC, rlim[0], rlim[1],
        zlim[0], zlim[1]);
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
        bfield->psi.c [0:bfield->psi.n_x * bfield->psi.n_y * NSIZE_COMP2D],
        bfield->br.c [0:bfield->br.n_x * bfield->br.n_y * NSIZE_COMP2D],
        bfield->bphi.c [0:bfield->bphi.n_x * bfield->bphi.n_y * NSIZE_COMP2D],
        bfield->bz.c [0:bfield->bz.n_x * bfield->bz.n_y * NSIZE_COMP2D])
}

a5err BfieldSpline2D_eval_psi(
    real psi[1], real r, real phi, real z, BfieldSpline2D *bfield)
{

    int interperr = 0;
    interperr += interp2Dcomp_eval_f(&psi[0], &bfield->psi, r, z);

    a5err err = 0;
    if (interperr)
    {
        err = error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_B_2DS);
    }
    return err;
}

a5err BfieldSpline2D_eval_psi_dpsi(
    real psi_dpsi[4], real r, real phi, real z, BfieldSpline2D *bfield)
{

    int interperr = 0;
    real psi_dpsi_temp[6];

    interperr += interp2Dcomp_eval_df(psi_dpsi_temp, &bfield->psi, r, z);
    psi_dpsi[0] = psi_dpsi_temp[0];
    psi_dpsi[1] = psi_dpsi_temp[1];
    psi_dpsi[2] = 0;
    psi_dpsi[3] = psi_dpsi_temp[2];

    a5err err = 0;
    if (interperr)
    {
        err = error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_B_2DS);
    }

    return err;
}

a5err BfieldSpline2D_eval_rho_drho(
    real rho_drho[4], real r, real phi, real z, BfieldSpline2D *bfield)
{
    int interperr = 0;
    real psi_dpsi[6];

    interperr += interp2Dcomp_eval_df(psi_dpsi, &bfield->psi, r, z);

    a5err err = 0;
    if (interperr)
    {
        err = error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_B_2DS);
    }

    /* Check that the values seem valid */
    real delta = bfield->psilimits[1] - bfield->psilimits[0];
    if (!err && (psi_dpsi[0] - bfield->psilimits[0]) / delta < 0)
    {
        err = error_raise(ERR_INPUT_UNPHYSICAL, __LINE__, EF_B_2DS);
    }

    /* Normalize psi to get rho */
    rho_drho[0] = sqrt((psi_dpsi[0] - bfield->psilimits[0]) / delta);
    rho_drho[1] = psi_dpsi[1] / (2 * delta * rho_drho[0]);
    rho_drho[2] = 0;
    rho_drho[3] = psi_dpsi[2] / (2 * delta * rho_drho[0]);

    return err;
}

a5err BfieldSpline2D_eval_b(
    real b[3], real r, real phi, real z, BfieldSpline2D *bfield)
{
    a5err err = 0;
    int interperr = 0;

    interperr += interp2Dcomp_eval_f(&b[0], &bfield->br, r, z);
    interperr += interp2Dcomp_eval_f(&b[1], &bfield->bphi, r, z);
    interperr += interp2Dcomp_eval_f(&b[2], &bfield->bz, r, z);

    if (interperr)
    {
        err = error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_B_2DS);
    }

    if (!err)
    {
        real psi_dpsi[6];
        interperr += interp2Dcomp_eval_df(psi_dpsi, &bfield->psi, r, z);
        b[0] = b[0] - psi_dpsi[2] / r;
        b[2] = b[2] + psi_dpsi[1] / r;

        if (interperr)
        {
            err = error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_B_2DS);
        }
    }

    return err;
}

a5err BfieldSpline2D_eval_b_db(
    real b_db[12], real r, real phi, real z, BfieldSpline2D *bfield)
{
    a5err err = 0;
    int interperr = 0;
    real b_db_temp[6];

    interperr += interp2Dcomp_eval_df(b_db_temp, &bfield->br, r, z);

    b_db[0] = b_db_temp[0];
    b_db[1] = b_db_temp[1];
    b_db[2] = 0;
    b_db[3] = b_db_temp[2];

    interperr += interp2Dcomp_eval_df(b_db_temp, &bfield->bphi, r, z);

    b_db[4] = b_db_temp[0];
    b_db[5] = b_db_temp[1];
    b_db[6] = 0;
    b_db[7] = b_db_temp[2];

    interperr += interp2Dcomp_eval_df(b_db_temp, &bfield->bz, r, z);

    b_db[8] = b_db_temp[0];
    b_db[9] = b_db_temp[1];
    b_db[10] = 0;
    b_db[11] = b_db_temp[2];

    if (interperr)
    {
        err = error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_B_2DS);
    }

    real psi_dpsi[6];

    if (!err)
    {
        interperr += interp2Dcomp_eval_df(psi_dpsi, &bfield->psi, r, z);

        b_db[0] = b_db[0] - psi_dpsi[2] / r;
        b_db[1] = b_db[1] + psi_dpsi[2] / (r * r) - psi_dpsi[5] / r;
        b_db[3] = b_db[3] - psi_dpsi[4] / r;
        b_db[8] = b_db[8] + psi_dpsi[1] / r;
        b_db[9] = b_db[9] - psi_dpsi[1] / (r * r) + psi_dpsi[3] / r;
        b_db[11] = b_db[11] + psi_dpsi[5] / r;

        if (interperr)
        {
            err = error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_B_2DS);
        }
    }

    int check = 0;
    check += ((b_db[0] * b_db[0] + b_db[4] * b_db[4] + b_db[8] * b_db[8]) == 0);
    if (!err && check)
    {
        err = error_raise(ERR_INPUT_UNPHYSICAL, __LINE__, EF_B_2DS);
    }

    return err;
}

a5err BfieldSpline2D_eval_axisrz(real axisrz[2], BfieldSpline2D *bfield)
{
    axisrz[0] = bfield->axisrz[0];
    axisrz[1] = bfield->axisrz[1];

    return 0;
}
