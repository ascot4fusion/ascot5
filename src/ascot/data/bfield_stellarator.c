/**
 * Implements bfield_stellarator.h.
 */
#include "bfield_stellarator.h"
#include "bfield.h"
#include "consts.h"
#include "defines.h"
#include "utils/interp.h"
#include "utils/mathlib.h"
#include <math.h>
#include <stdlib.h>

int BfieldStellarator_init(
    BfieldStellarator *bfield, size_t pnr, size_t pnz, size_t pnphi, size_t bnr,
    size_t bnz, size_t bnphi, size_t naxis, real prlim[2], real pzlim[2],
    real pphilim[2], real brlim[2], real bzlim[2], real bphilim[2],
    real axislim[2], real axisr[naxis], real axisz[naxis], real psilimits[2],
    real psi[pnr * pnz * pnphi], real br[bnr * bnz * bnphi],
    real bz[bnr * bnz * bnphi], real bphi[bnr * bnz * bnphi])
{
    int err = 0;
    bfield->psilimits[0] = psilimits[0];
    bfield->psilimits[1] = psilimits[1];
    err += interp3Dcomp_setup(
        &bfield->psi, psi, pnr, pnphi, pnz, NATURALBC, PERIODICBC, NATURALBC,
        prlim[0], prlim[1], pphilim[0], pphilim[1], pzlim[0], pzlim[1]);
    err += interp3Dcomp_setup(
        &bfield->br, br, bnr, bnphi, bnz, NATURALBC, PERIODICBC, NATURALBC,
        brlim[0], brlim[1], bphilim[0], bphilim[1], bzlim[0], bzlim[1]);
    err += interp3Dcomp_setup(
        &bfield->bz, bz, bnr, bnphi, bnz, NATURALBC, PERIODICBC, NATURALBC,
        brlim[0], brlim[1], bphilim[0], bphilim[1], bzlim[0], bzlim[1]);
    err += interp3Dcomp_setup(
        &bfield->bphi, bphi, bnr, bnphi, bnz, NATURALBC, PERIODICBC, NATURALBC,
        brlim[0], brlim[1], bphilim[0], bphilim[1], bzlim[0], bzlim[1]);

    real *c1 = (real *)malloc(naxis * sizeof(real));
    real *c2 = (real *)malloc(naxis * sizeof(real));
    for (size_t i = 0; i < naxis; i++)
    {
        c1[i] = axisr[i];
        c2[i] = axisz[i];
    }
    linint1D_init(
        &bfield->axisr, c1, naxis, PERIODICBC, axislim[0], axislim[1]);
    linint1D_init(
        &bfield->axisz, c2, naxis, PERIODICBC, axislim[0], axislim[1]);
    return 0;
}

void BfieldStellarator_free(BfieldStellarator *bfield)
{
    free(bfield->psi.c);
    free(bfield->br.c);
    free(bfield->bphi.c);
    free(bfield->bz.c);
}

void BfieldStellarator_offload(BfieldStellarator *bfield)
{
    (void)bfield;
    GPU_MAP_TO_DEVICE(
        bfield->axisr, bfield->axisr.c [0:bfield->axisr.n_x], bfield->axisz,
        bfield->axisz.c [0:bfield->axisz.n_x], bfield->psi, bfield->br,
        bfield->bphi, bfield->bz,
        bfield->psi.c
        [0:bfield->psi.n_x * bfield->psi.n_y * bfield->psi.n_z * NSIZE_COMP3D],
        bfield->br.c
        [0:bfield->br.n_x * bfield->br.n_y * bfield->br.n_z * NSIZE_COMP3D],
        bfield->bphi.c [0:bfield->bphi.n_x * bfield->bphi.n_y *
                          bfield->bphi.n_z * NSIZE_COMP3D],
        bfield->bz.c
        [0:bfield->bz.n_x * bfield->bz.n_y * bfield->bz.n_z * NSIZE_COMP3D])
}

err_t BfieldStellarator_eval_psi(
    real psi[1], real r, real phi, real z, BfieldStellarator *bfield)
{
    err_t err = 0;
    int interperr = interp3Dcomp_eval_f(&psi[0], &bfield->psi, r, phi, z);

    err = ERROR_CHECK(
        err, interperr, ERR_INTERPOLATED_OUTSIDE_RANGE,
        DATA_BFIELD_STELLARATOR_C);
    return err;
}

err_t BfieldStellarator_eval_psi_dpsi(
    real psi_dpsi[4], real r, real phi, real z, BfieldStellarator *bfield)
{
    err_t err = 0;
    real psi_dpsi_temp[10];
    int interperr =
        interp3Dcomp_eval_df(psi_dpsi_temp, &bfield->psi, r, phi, z);
    psi_dpsi[0] = psi_dpsi_temp[0];
    psi_dpsi[1] = psi_dpsi_temp[1];
    psi_dpsi[2] = psi_dpsi_temp[2];
    psi_dpsi[3] = psi_dpsi_temp[3];

    err = ERROR_CHECK(
        err, interperr, ERR_INTERPOLATED_OUTSIDE_RANGE,
        DATA_BFIELD_STELLARATOR_C);
    return err;
}

err_t BfieldStellarator_eval_b(
    real b[3], real r, real phi, real z, BfieldStellarator *bfield)
{
    err_t err = 0;
    int interperr = 0;
    interperr += interp3Dcomp_eval_f(&b[0], &bfield->br, r, phi, z);
    interperr += interp3Dcomp_eval_f(&b[1], &bfield->bphi, r, phi, z);
    interperr += interp3Dcomp_eval_f(&b[2], &bfield->bz, r, phi, z);

    err = ERROR_CHECK(
        err, interperr, ERR_INTERPOLATED_OUTSIDE_RANGE,
        DATA_BFIELD_STELLARATOR_C);
    return 0;
}

err_t BfieldStellarator_eval_b_db(
    real b_db[12], real r, real phi, real z, BfieldStellarator *bfield)
{
    err_t err = 0;
    int interperr = 0;
    real b_db_temp[10];

    interperr += interp3Dcomp_eval_df(b_db_temp, &bfield->br, r, phi, z);
    b_db[0] = b_db_temp[0];
    b_db[1] = b_db_temp[1];
    b_db[2] = b_db_temp[2];
    b_db[3] = b_db_temp[3];

    interperr += interp3Dcomp_eval_df(b_db_temp, &bfield->bphi, r, phi, z);
    b_db[4] = b_db_temp[0];
    b_db[5] = b_db_temp[1];
    b_db[6] = b_db_temp[2];
    b_db[7] = b_db_temp[3];

    interperr += interp3Dcomp_eval_df(b_db_temp, &bfield->bz, r, phi, z);
    b_db[8] = b_db_temp[0];
    b_db[9] = b_db_temp[1];
    b_db[10] = b_db_temp[2];
    b_db[11] = b_db_temp[3];

    err = ERROR_CHECK(
        err, interperr, ERR_INTERPOLATED_OUTSIDE_RANGE,
        DATA_BFIELD_STELLARATOR_C);
    return 0;
}

err_t BfieldStellarator_eval_axisrz(
    real axisrz[2], real phi, BfieldStellarator *bfield)
{
    err_t err = 0;
    int interperr = 0;
    interperr += linint1D_eval_f(&axisrz[0], &bfield->axisr, phi);
    interperr += linint1D_eval_f(&axisrz[1], &bfield->axisz, phi);

    err = ERROR_CHECK(
        err, interperr, ERR_INTERPOLATED_OUTSIDE_RANGE,
        DATA_BFIELD_STELLARATOR_C);
    return err;
}
