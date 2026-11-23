/**
 * Implements "plasma_dynamic1d.h".
 */
#include "plasma_dynamic1d.h"
#include "consts.h"
#include "defines.h"
#include "plasma.h"
#include <stdio.h>
#include <stdlib.h>

int PlasmaDynamic1D_init(
    PlasmaDynamic1D *plasma, size_t nrho, size_t ntime, size_t nion,
    int anum[nion], int znum[nion], int charge[nion], real mass[nion],
    real rho[nrho], real time[ntime], real Te[nrho * ntime],
    real Ti[nrho * ntime], real ne[nrho * ntime], real ni[nrho * ntime * nion],
    real vtor[nrho * ntime])
{

    plasma->nrho = nrho;
    plasma->ntime = ntime;
    plasma->nspecies = nion + 1;

    plasma->anum = (int *)malloc(nion * sizeof(int));
    plasma->znum = (int *)malloc(nion * sizeof(int));
    plasma->mass = (real *)malloc((nion + 1) * sizeof(real));
    plasma->charge = (int *)malloc((nion + 1) * sizeof(int));

    plasma->mass[0] = CONST_M_E;
    plasma->charge[0] = -1;
    for (size_t i = 0; i < nion; i++)
    {
        plasma->znum[i] = znum[i];
        plasma->anum[i] = anum[i];
        plasma->mass[i + 1] = mass[i];
        plasma->charge[i + 1] = charge[i];
    }
    plasma->rho = (real *)malloc(nrho * sizeof(real));
    for (size_t i = 0; i < nrho; i++)
        plasma->rho[i] = rho[i];

    plasma->time = (real *)malloc(ntime * sizeof(real));
    for (size_t i = 0; i < ntime; i++)
        plasma->time[i] = time[i];

    plasma->vtor = (real *)malloc(nrho * ntime * sizeof(real));
    plasma->density = (real *)malloc((nion + 1) * nrho * ntime * sizeof(real));
    plasma->temperature = (real *)malloc(2 * nrho * ntime * sizeof(real));
    for (size_t i = 0; i < nrho; i++)
    {
        for (size_t j = 0; j < ntime; j++)
        {
            plasma->vtor[i * ntime + j] = vtor[i * ntime + j];
            plasma->density[i * ntime + j] = ne[i * ntime + j];
            plasma->temperature[i * ntime + j] = Te[i * ntime + j];
            plasma->temperature[nrho * ntime + i * ntime + j] =
                Ti[i * ntime + j];
            for (size_t k = 0; k < nion; k++)
            {
                plasma->density[(k + 1) * nrho * ntime + i * ntime + j] =
                    ni[k * nrho * ntime + i * ntime + j];
            }
        }
    }
    return 0;
}

void PlasmaDynamic1D_free(PlasmaDynamic1D *plasma)
{
    free(plasma->mass);
    free(plasma->charge);
    free(plasma->anum);
    free(plasma->znum);
    free(plasma->rho);
    free(plasma->time);
    free(plasma->temperature);
    free(plasma->density);
}

void PlasmaDynamic1D_offload(PlasmaDynamic1D *plasma)
{
    SUPPRESS_UNUSED_WARNING(plasma);
    GPU_MAP_TO_DEVICE(
        plasma->mass [0:plasma->nspecies], plasma->charge [0:plasma->nspecies],
        plasma->anum [0:plasma->nspecies - 1],
        plasma->znum [0:plasma->nspecies - 1], plasma->rho [0:plasma->nrho],
        plasma->time [0:plasma->ntime],
        plasma->vtor [0:plasma->ntime * plasma->nrho],
        plasma->temperature [0:plasma->ntime * plasma->nrho * 2],
        plasma->density [0:plasma->nrho * plasma->nspecies * plasma->ntime])
}

err_t PlasmaDynamic1D_eval_temperature(
    real temperature[1], real rho, real t, size_t i_species,
    PlasmaDynamic1D *plasma)
{
    real n[MAX_SPECIES], T[MAX_SPECIES];
    err_t err = PlasmaDynamic1D_eval_nT(n, T, rho, t, plasma);
    temperature[0] = T[i_species];

    return err;
}

err_t PlasmaDynamic1D_eval_density(
    real density[1], real rho, real t, size_t i_species,
    PlasmaDynamic1D *plasma)
{
    real n[MAX_SPECIES], T[MAX_SPECIES];
    err_t err = PlasmaDynamic1D_eval_nT(n, T, rho, t, plasma);
    density[0] = n[i_species];

    return err;
}

err_t PlasmaDynamic1D_eval_nT(
    real *density, real *temperature, real rho, real t, PlasmaDynamic1D *plasma)
{
    err_t err = 0;
    err = ERROR_CHECK(
        err,
        (rho < plasma->rho[0]) || (rho >= plasma->rho[plasma->nrho - 1]) ||
            (t < plasma->time[0]) || (t >= plasma->time[plasma->ntime - 1]),
        ERR_INTERPOLATED_OUTSIDE_RANGE, DATA_PLASMA_DYNAMIC1D_C);

    size_t i_rho = 0;
    for (size_t i = 0; i < plasma->nrho; ++i)
        i_rho += (plasma->rho[i] <= rho);
    i_rho--;

    size_t i_time = 0;
    for (size_t i = 0; i < plasma->ntime; ++i)
        i_time += (plasma->time[i] <= t);
    i_time--;

    i_rho = err ? 0 : i_rho;
    i_time = err ? 0 : i_time;

    real t_rho = (rho - plasma->rho[i_rho]) /
                 (plasma->rho[i_rho + 1] - plasma->rho[i_rho]);
    real t_time = (t - plasma->time[i_time]) /
                  (plasma->time[i_time + 1] - plasma->time[i_time]);

    for (size_t i = 0; i < plasma->nspecies; i++)
    {
        real p11, p12, p21, p22, p1, p2;

        /* Density: [species][rho][time] */
        size_t offset =
            i * plasma->nrho * plasma->ntime + i_rho * plasma->ntime + i_time;

        p11 = plasma->density[offset];
        p12 = plasma->density[offset + 1];
        p21 = plasma->density[offset + plasma->ntime];
        p22 = plasma->density[offset + plasma->ntime + 1];

        p1 = p11 + t_time * (p12 - p11);
        p2 = p21 + t_time * (p22 - p21);

        density[i] = p1 + t_rho * (p2 - p1);

        /* Temperature: [component][rho][time] */
        size_t temp_component = (i == 0) ? 0 : 1;

        offset = temp_component * plasma->nrho * plasma->ntime +
                 i_rho * plasma->ntime + i_time;

        p11 = plasma->temperature[offset];
        p12 = plasma->temperature[offset + 1];
        p21 = plasma->temperature[offset + plasma->ntime];
        p22 = plasma->temperature[offset + plasma->ntime + 1];

        p1 = p11 + t_time * (p12 - p11);
        p2 = p21 + t_time * (p22 - p21);

        temperature[i] = p1 + t_rho * (p2 - p1);
    }

    return err;
}

err_t PlasmaDynamic1D_eval_flow(
    real vflow[1], real rho, real t, real r, PlasmaDynamic1D *plasma)
{
    err_t err = 0;
    err = ERROR_CHECK(
        err,
        (rho < plasma->rho[0]) || (rho >= plasma->rho[plasma->nrho - 1]) ||
            (t < plasma->time[0]) || (t >= plasma->time[plasma->ntime - 1]),
        ERR_INTERPOLATED_OUTSIDE_RANGE, DATA_PLASMA_DYNAMIC1D_C);

    size_t i_rho = 0;
    for (size_t i = 0; i < plasma->nrho; ++i)
        i_rho += (plasma->rho[i] <= rho);
    i_rho--;

    size_t i_time = 0;
    for (size_t i = 0; i < plasma->ntime; ++i)
        i_time += (plasma->time[i] <= t);
    i_time--;

    i_rho = err ? 0 : i_rho;
    i_time = err ? 0 : i_time;

    real t_rho = (rho - plasma->rho[i_rho]) /
                 (plasma->rho[i_rho + 1] - plasma->rho[i_rho]);

    real t_time = (t - plasma->time[i_time]) /
                  (plasma->time[i_time + 1] - plasma->time[i_time]);

    real p11, p12, p21, p22, p1, p2;
    p11 = plasma->vtor[i_rho * plasma->ntime + i_time];
    p12 = plasma->vtor[(i_rho + 1) * plasma->ntime + i_time];
    p21 = plasma->vtor[i_rho * plasma->ntime + i_time + 1];
    p22 = plasma->vtor[(i_rho + 1) * plasma->ntime + i_time + 1];

    p1 = p11 + t_rho * (p12 - p11);
    p2 = p21 + t_rho * (p22 - p21);
    vflow[0] = p1 + t_time * (p2 - p1);

    vflow[0] *= r;
    return err;
}
