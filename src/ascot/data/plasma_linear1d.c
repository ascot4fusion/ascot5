/**
 * Implements plasma_linear1d.h.
 */
#include "plasma_linear1d.h"
#include "consts.h"
#include "defines.h"
#include "plasma.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int PlasmaLinear1D_init(
    PlasmaLinear1D *plasma, size_t nrho, size_t nion, int anum[nion],
    int znum[nion], real mass[nion], real charge[nion], real rho[nrho],
    real Te[nrho], real Ti[nrho], real ne[nrho], real ni[nrho * nion],
    real vtor[nrho])
{

    plasma->nrho = nrho;
    plasma->nspecies = nion + 1;

    plasma->anum = (int *)malloc(nion * sizeof(int));
    plasma->znum = (int *)malloc(nion * sizeof(int));
    plasma->mass = (real *)malloc((nion + 1) * sizeof(real));
    plasma->charge = (real *)malloc((nion + 1) * sizeof(real));
    for (size_t i = 0; i < plasma->nspecies; i++)
    {
        if (i < nion)
        {
            plasma->znum[i] = znum[i];
            plasma->anum[i] = anum[i];
        }
        plasma->mass[i] = mass[i];
        plasma->charge[i] = charge[i];
    }
    plasma->rho = (real *)malloc(nrho * sizeof(real));
    plasma->vtor = (real *)malloc(nrho * sizeof(real));
    plasma->temp = (real *)malloc(2 * nrho * sizeof(real));
    plasma->dens = (real *)malloc((nion + 1) * nrho * sizeof(real));
    for (size_t i = 0; i < plasma->nrho; i++)
    {
        plasma->rho[i] = rho[i];
        plasma->vtor[i] = vtor[i];
        plasma->temp[i] = Te[i];
        plasma->temp[nrho + i] = Ti[i];
        plasma->dens[i] = ne[i];
        for (size_t j = 0; j < nion; j++)
        {
            plasma->dens[(j + 1) * nrho + i] = ni[j * nrho + i];
        }
    }
    return 0;
}

void PlasmaLinear1D_free(PlasmaLinear1D *plasma)
{
    free(plasma->mass);
    free(plasma->charge);
    free(plasma->anum);
    free(plasma->znum);
    free(plasma->rho);
    free(plasma->temp);
    free(plasma->dens);
}

void PlasmaLinear1D_offload(PlasmaLinear1D *plasma)
{
    SUPPRESS_UNUSED_WARNING(plasma);
    GPU_MAP_TO_DEVICE(
        plasma->rho [0:plasma->nrho], plasma->vtor [0:plasma->nrho],
        plasma->temp [0:2 * plasma->nrho], plasma->mass [0:plasma->nspecies],
        plasma->charge [0:plasma->nspecies],
        plasma->anum [0:plasma->nspecies - 1],
        plasma->znum [0:plasma->nspecies - 1],
        plasma->dens [0:plasma->nrho * plasma->nspecies])
}

err_t PlasmaLinear1D_eval_temp(
    real *temp, real rho, size_t species, PlasmaLinear1D *plasma)
{

    err_t err = 0;
    err = ERROR_CHECK(
        err, rho < plasma->rho[0] | rho >= plasma->rho[plasma->nrho - 1],
        ERR_INTERPOLATED_OUTSIDE_RANGE, DATA_PLASMA_LINEAR1D_C);

    size_t i_rho = 0;
    for (size_t i = 1; i < plasma->nrho; ++i)
        i_rho += (plasma->rho[i] <= rho);
    i_rho--;

    real t_rho = (rho - plasma->rho[i_rho]) /
                 (plasma->rho[i_rho + 1] - plasma->rho[i_rho]);

    real p1 = plasma->temp[(species > 0) * plasma->nrho + i_rho];
    real p2 = plasma->temp[(species > 0) * plasma->nrho + i_rho + 1];
    temp[0] = p1 + t_rho * (p2 - p1);

    return err;
}

err_t PlasmaLinear1D_eval_dens(
    real *dens, real rho, size_t species, PlasmaLinear1D *plasma)
{

    err_t err = 0;
    err = ERROR_CHECK(
        err, rho < plasma->rho[0] | rho >= plasma->rho[plasma->nrho - 1],
        ERR_INTERPOLATED_OUTSIDE_RANGE, DATA_PLASMA_LINEAR1D_C);

    size_t i_rho = 0;
    for (size_t i = 1; i < plasma->nrho; ++i)
        i_rho += (plasma->rho[i] <= rho);
    i_rho--;

    real t_rho = (rho - plasma->rho[i_rho]) /
                 (plasma->rho[i_rho + 1] - plasma->rho[i_rho]);

    real p1 = plasma->dens[species * plasma->nrho + i_rho];
    real p2 = plasma->dens[species * plasma->nrho + i_rho + 1];
    dens[0] = p1 + t_rho * (p2 - p1);

    return err;
}

err_t PlasmaLinear1D_eval_densandtemp(
    real *dens, real *temp, real rho, PlasmaLinear1D *plasma)
{
    err_t err = 0;
    err_t err = 0;
    err = ERROR_CHECK(
        err, rho < plasma->rho[0] | rho >= plasma->rho[plasma->nrho - 1],
        ERR_INTERPOLATED_OUTSIDE_RANGE, DATA_PLASMA_LINEAR1D_C);

    size_t i_rho = 0;
    for (size_t i = 1; i < plasma->nrho; ++i)
        i_rho += (plasma->rho[i] <= rho);
    i_rho--;

    real t_rho = (rho - plasma->rho[i_rho]) /
                 (plasma->rho[i_rho + 1] - plasma->rho[i_rho]);

    real p1, p2;
    for (size_t i = 0; i < plasma->nspecies; i++)
    {
        p1 = plasma->dens[i * plasma->nrho + i_rho];
        p2 = plasma->dens[i * plasma->nrho + i_rho + 1];
        dens[i] = p1 + t_rho * (p2 - p1);

        if (i < 2)
        {
            /* Electron and ion temperature */
            p1 = plasma->temp[i * plasma->nrho + i_rho];
            p2 = plasma->temp[i * plasma->nrho + i_rho + 1];
            temp[i] = p1 + t_rho * (p2 - p1);
        }
        else
        {
            /* Temperature is same for all ion species */
            temp[i] = temp[1];
        }
    }

    return err;
}

err_t PlasmaLinear1D_eval_flow(
    real vflow[1], real rho, real r, PlasmaLinear1D *plasma)
{
    err_t err = 0;
    err = ERROR_CHECK(
        err, rho < plasma->rho[0] | rho >= plasma->rho[plasma->nrho - 1],
        ERR_INTERPOLATED_OUTSIDE_RANGE, DATA_PLASMA_LINEAR1D_C);

    size_t i_rho = 0;
    for (size_t i = 1; i < plasma->nrho; ++i)
        i_rho += (plasma->rho[i] <= rho);
    i_rho--;

    real t_rho = (rho - plasma->rho[i_rho]) /
                 (plasma->rho[i_rho + 1] - plasma->rho[i_rho]);
    real p1 = plasma->vtor[i_rho];
    real p2 = plasma->vtor[i_rho + 1];
    vflow[0] = (p1 + t_rho * (p2 - p1)) * r;

    return err;
}
