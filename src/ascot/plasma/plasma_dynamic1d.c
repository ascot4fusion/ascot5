/**
 * Implements plasma_dynamic1d.h.
 */
#include "plasma_dynamic1d.h"
#include "defines.h"
#include "consts.h"
#include "errors.h"
#include "plasma.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int PlasmaDynamic1D_init(
    PlasmaDynamic1D *plasma, int nrho, int ntime, int nion, int *anum,
    int *znum, real *rho, real *time, real *mass, real *charge, real *Te,
    real *Ti, real *ne, real *ni, real *vtor)
{

    plasma->nrho = nrho;
    plasma->ntime = ntime;
    plasma->nspecies = nion + 1;

    plasma->anum = (int *)malloc(nion * sizeof(int));
    plasma->znum = (int *)malloc(nion * sizeof(int));
    plasma->mass = (real *)malloc((nion + 1) * sizeof(real));
    plasma->charge = (real *)malloc((nion + 1) * sizeof(real));
    for (int i = 0; i < plasma->nspecies; i++)
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
    for (int i = 0; i < nrho; i++)
    {
        plasma->rho[i] = rho[i];
    }
    plasma->time = (real *)malloc(ntime * sizeof(real));
    for (int i = 0; i < ntime; i++)
    {
        plasma->time[i] = time[i];
    }
    plasma->vtor = (real *)malloc(nrho * ntime * sizeof(real));
    plasma->temp = (real *)malloc(2 * nrho * ntime * sizeof(real));
    plasma->dens = (real *)malloc((nion + 1) * nrho * ntime * sizeof(real));
    for (int i = 0; i < nrho; i++)
    {
        for (int j = 0; j < ntime; j++)
        {
            plasma->vtor[j * nrho + i] = vtor[j * nrho + i];
            plasma->temp[j * 2 * nrho + i] = Te[j * nrho + i];
            plasma->temp[(j * 2 + 1) * nrho + i] = Ti[j * nrho + i];
            plasma->dens[j * nrho + i] = ne[j * nrho + i];
            for (int k = 0; k < nion; k++)
            {
                plasma->dens[(k + 1) * nrho * ntime + j * nrho + i] =
                    ni[k * nrho * ntime + j * nrho + i];
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
    free(plasma->temp);
    free(plasma->dens);
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
        plasma->temp [0:plasma->ntime * plasma->nrho * 2],
        plasma->dens [0:plasma->nrho * plasma->nspecies * plasma->ntime])
}

a5err PlasmaDynamic1D_eval_temp(
    real *temp, real rho, real t, int species, PlasmaDynamic1D *plasma)
{

    real temp_dens[MAX_SPECIES], temp_temp[MAX_SPECIES];

    a5err err =
        PlasmaDynamic1D_eval_densandtemp(temp_dens, temp_temp, rho, t, plasma);

    *temp = temp_temp[species];

    return err;
}

a5err PlasmaDynamic1D_eval_dens(
    real *dens, real rho, real t, int species, PlasmaDynamic1D *plasma)
{
    real temp_dens[MAX_SPECIES], temp_temp[MAX_SPECIES];

    a5err err =
        PlasmaDynamic1D_eval_densandtemp(temp_dens, temp_temp, rho, t, plasma);

    *dens = temp_dens[species];

    return err;
}

a5err PlasmaDynamic1D_eval_densandtemp(
    real *dens, real *temp, real rho, real t, PlasmaDynamic1D *plasma)
{

    a5err err = 0;
    if (rho < plasma->rho[0])
    {
        err = error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_PLASMA_1D);
    }
    else if (rho >= plasma->rho[plasma->nrho - 1])
    {
        err = error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_PLASMA_1D);
    }
    else
    {
        int i_rho = 0;
        while (i_rho < plasma->nrho - 1 && plasma->rho[i_rho] <= rho)
        {
            i_rho++;
        }
        i_rho--;

        real t_rho = (rho - plasma->rho[i_rho]) /
                     (plasma->rho[i_rho + 1] - plasma->rho[i_rho]);

        int i_time = 0;
        while (i_time < plasma->ntime - 1 && plasma->time[i_time] <= t)
        {
            i_time++;
        }
        i_time--;

        real t_time = (t - plasma->time[i_time]) /
                      (plasma->time[i_time + 1] - plasma->time[i_time]);

        if (i_time < 0)
        {
            /* time < t[0], use first profile */
            i_time = 0;
            t_time = 0;
        }
        else if (i_time >= plasma->ntime - 2)
        {
            /* time > t[ntime-1], use last profile */
            i_time = plasma->ntime - 2;
            t_time = 1;
        }

        for (int i = 0; i < plasma->nspecies; i++)
        {
            real p11, p12, p21, p22, p1, p2;

            p11 = plasma->dens
                      [i_time * plasma->nspecies * plasma->nrho +
                       i * plasma->nrho + i_rho];
            p12 = plasma->dens
                      [i_time * plasma->nspecies * plasma->nrho +
                       i * plasma->nrho + i_rho + 1];
            p21 = plasma->dens
                      [(i_time + 1) * plasma->nspecies * plasma->nrho +
                       i * plasma->nrho + i_rho];
            p22 = plasma->dens
                      [(i_time + 1) * plasma->nspecies * plasma->nrho +
                       i * plasma->nrho + i_rho + 1];

            p1 = p11 + t_rho * (p12 - p11);
            p2 = p21 + t_rho * (p22 - p21);

            dens[i] = p1 + t_time * (p2 - p1);

            if (i < 2)
            {
                /* Electron and ion temperature */
                p11 =
                    plasma->temp
                        [i_time * 2 * plasma->nrho + i * plasma->nrho + i_rho];
                p12 = plasma->temp
                          [i_time * 2 * plasma->nrho + i * plasma->nrho +
                           i_rho + 1];
                p21 = plasma->temp
                          [(i_time + 1) * 2 * plasma->nrho + i * plasma->nrho +
                           i_rho];
                p22 = plasma->temp
                          [(i_time + 1) * 2 * plasma->nrho + i * plasma->nrho +
                           i_rho + 1];

                p1 = p11 + t_rho * (p12 - p11);
                p2 = p21 + t_rho * (p22 - p21);

                temp[i] = p1 + t_time * (p2 - p1);
            }
            else
            {
                /* Temperature is same for all ion species */
                temp[i] = temp[1];
            }
        }
    }

    return err;
}

a5err PlasmaDynamic1D_eval_flow(
    real *vflow, real rho, real t, real r, PlasmaDynamic1D *plasma)
{
    a5err err = 0;
    if (rho < plasma->rho[0])
    {
        err = error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_PLASMA_1D);
    }
    else if (rho >= plasma->rho[plasma->nrho - 1])
    {
        err = error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_PLASMA_1D);
    }
    else
    {
        int i_rho = 0;
        while (i_rho < plasma->nrho - 1 && plasma->rho[i_rho] <= rho)
        {
            i_rho++;
        }
        i_rho--;

        real t_rho = (rho - plasma->rho[i_rho]) /
                     (plasma->rho[i_rho + 1] - plasma->rho[i_rho]);

        int i_time = 0;
        while (i_time < plasma->ntime - 1 && plasma->time[i_time] <= t)
        {
            i_time++;
        }
        i_time--;

        real t_time = (t - plasma->time[i_time]) /
                      (plasma->time[i_time + 1] - plasma->time[i_time]);

        if (i_time < 0)
        {
            /* time < t[0], use first profile */
            i_time = 0;
            t_time = 0;
        }
        else if (i_time >= plasma->ntime - 2)
        {
            /* time > t[ntime-1], use last profile */
            i_time = plasma->ntime - 2;
            t_time = 1;
        }

        real p11, p12, p21, p22, p1, p2;

        p11 = plasma->vtor[i_time * plasma->nrho + i_rho];
        p12 = plasma->vtor[(i_time + 1) * plasma->nrho + i_rho];
        p21 = plasma->vtor[i_time * plasma->nrho + i_rho + 1];
        p22 = plasma->vtor[(i_time + 1) * plasma->nrho + i_rho + 1];

        p1 = p11 + t_rho * (p12 - p11);
        p2 = p21 + t_rho * (p22 - p21);

        *vflow = p1 + t_time * (p2 - p1);
    }
    *vflow *= r;
    return err;
}
