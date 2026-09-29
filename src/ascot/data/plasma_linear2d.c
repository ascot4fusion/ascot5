/**
 * Implements "plasma_linear2d.h".
 */
#include "plasma_linear2d.h"
#include "consts.h"
#include "defines.h"
#include "plasma.h"
#include "utils/interp.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int PlasmaLinear2D_init(
    PlasmaLinear2D *plasma, size_t nr, size_t nz, size_t nion, real rlim[2],
    real zlim[2], int anum[nion], int znum[nion], real mass[nion],
    real charge[nion], real Te[nr * nz], real Ti[nr * nz], real ne[nr * nz],
    real ni[nr * nz * nion], real vtor[nr * nz])
{
    if (nr * nz == 0)
        return -1;
    plasma->nspecies = nion + 1;

    plasma->anum = (int *)malloc(nion * sizeof(int));
    plasma->znum = (int *)malloc(nion * sizeof(int));
    plasma->mass = (real *)malloc((nion + 1) * sizeof(real));
    plasma->charge = (real *)malloc((nion + 1) * sizeof(real));

    plasma->mass[0] = CONST_M_E;
    plasma->charge[0] = -CONST_E;
    for (size_t i = 0; i < nion; i++)
    {
        plasma->znum[i] = znum[i];
        plasma->anum[i] = anum[i];
        plasma->mass[i + 1] = mass[i];
        plasma->charge[i + 1] = charge[i];
    }
    plasma->density = (Linear2D *)malloc((nion + 1) * sizeof(Linear2D));
    real *c = malloc(nr * nz * sizeof(real));
    for (size_t j = 0; j < (nr * nz); j++)
        c[j] = ne[j];

    Linear2D_init(
        &plasma->density[0], nr, nz, NATURALBC, NATURALBC, rlim, zlim, c);

    c = (real *)malloc(nr * nz * sizeof(real));
    for (size_t j = 0; j < nr * nz; j++)
    {
        c[j] = Te[j];
    }
    Linear2D_init(
        &plasma->temperature[0], nr, nz, NATURALBC, NATURALBC, rlim, zlim, c);

    c = (real *)malloc(nr * nz * sizeof(real));
    for (size_t j = 0; j < nr * nz; j++)
    {
        c[j] = Ti[j];
    }
    Linear2D_init(
        &plasma->temperature[1], nr, nz, NATURALBC, NATURALBC, rlim, zlim, c);

    c = (real *)malloc(nr * nz * sizeof(real));
    for (size_t j = 0; j < nr * nz; j++)
    {
        c[j] = vtor[j];
    }
    Linear2D_init(&plasma->vtor, nr, nz, NATURALBC, NATURALBC, rlim, zlim, c);

    for (size_t i = 0; i < nion; i++)
    {

        c = (real *)malloc(nr * nz * sizeof(real));
        for (size_t j = 0; j < nr * nz; j++)
        {
            c[j] = ni[i * (nr * nz) + j];
        }
        Linear2D_init(
            &plasma->density[1 + i], nr, nz, NATURALBC, NATURALBC, rlim, zlim,
            c);
    }
    return 0;
}

void PlasmaLinear2D_free(PlasmaLinear2D *plasma)
{
    free(plasma->mass);
    free(plasma->charge);
    free(plasma->anum);
    free(plasma->znum);
    free(plasma->vtor.c);
    free(plasma->temperature[0].c);
    free(plasma->temperature[1].c);
    for (size_t i = 0; i < plasma->nspecies; i++)
        free(plasma->density[i].c);

    free(plasma->density);
}

void PlasmaLinear2D_offload(PlasmaLinear2D *plasma)
{
    SUPPRESS_UNUSED_WARNING(plasma);
    // TODO
}

err_t PlasmaLinear2D_eval_temperature(
    real temperature[1], real r, real z, size_t i_species,
    PlasmaLinear2D *plasma)
{
    err_t err = 0;
    int interperr = 0;
    int ision = i_species > 0;
    interperr +=
        Linear2D_eval_f(temperature, &plasma->temperature[ision], r, z);
    err = ERROR_CHECK(
        err, interperr, ERR_INTERPOLATED_OUTSIDE_RANGE, DATA_PLASMA_LINEAR2D_C);
    return err;
}

err_t PlasmaLinear2D_eval_density(
    real density[1], real r, real z, size_t i_species, PlasmaLinear2D *plasma)
{
    err_t err = 0;
    int interperr = 0;
    interperr += Linear2D_eval_f(density, &plasma->density[i_species], r, z);
    err = ERROR_CHECK(
        err, interperr, ERR_INTERPOLATED_OUTSIDE_RANGE, DATA_PLASMA_LINEAR2D_C);
    return err;
}

err_t PlasmaLinear2D_eval_nT(
    real *density, real *temperature, real r, real z, PlasmaLinear2D *plasma)
{
    err_t err = 0;
    int interperr = 0;
    for (size_t i = 0; i < plasma->nspecies; i++)
    {
        int ision = i > 0;
        interperr +=
            Linear2D_eval_f(&temperature[i], &plasma->temperature[ision], r, z);
        interperr += Linear2D_eval_f(&density[i], &plasma->density[i], r, z);
    }
    err = ERROR_CHECK(
        err, interperr, ERR_INTERPOLATED_OUTSIDE_RANGE, DATA_PLASMA_LINEAR2D_C);
    return err;
}

err_t PlasmaLinear2D_eval_flow(
    real vflow[1], real r, real z, PlasmaLinear2D *plasma)
{
    err_t err = 0;
    int interperr = Linear2D_eval_f(vflow, &plasma->vtor, r, z);
    err = ERROR_CHECK(
        err, interperr, ERR_INTERPOLATED_OUTSIDE_RANGE, DATA_PLASMA_LINEAR2D_C);
    *vflow *= r;
    return err;
}
