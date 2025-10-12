/**
 * Implements neutral_radial.h.
 */
#include "defines.h"
#include "neutral.h"
#include "utils/interp.h"
#include "utils/mathlib.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int NeutralRadial_init(
    NeutralRadial *neutral, size_t n, size_t nrho, int *anum, int *znum,
    real rholim[2], real *density, real *temperature)
{

    neutral->n = n;
    neutral->anum = (int *)malloc(n * sizeof(int));
    neutral->znum = (int *)malloc(n * sizeof(int));
    neutral->n0 = (Linear1D *)malloc(n * sizeof(Linear1D));
    neutral->T0 = (Linear1D *)malloc(n * sizeof(Linear1D));
    int err = 0;
    for (size_t i = 0; i < neutral->n; i++)
    {
        neutral->anum[i] = anum[i];
        neutral->znum[i] = znum[i];

        real *c = (real *)malloc(nrho * sizeof(real));
        err += c == NULL ? 1 : 0;
        for (size_t i = 0; i < nrho; i++)
        {
            c[i] = density[i];
        }
        linint1D_init(
            &neutral->n0[i], c, nrho, NATURALBC, rholim[0], rholim[1]);
        c = (real *)malloc(nrho * sizeof(real));
        err += c == NULL ? 1 : 0;
        for (size_t i = 0; i < nrho; i++)
        {
            c[i] = temperature[i];
        }
        linint1D_init(
            &neutral->T0[i], c, nrho, NATURALBC, rholim[0], rholim[1]);
    }
    return err;
}

void NeutralRadial_free(NeutralRadial *neutral)
{
    free(neutral->anum);
    free(neutral->znum);
    for (size_t i = 0; i < neutral->n; i++)
    {
        free(neutral->n0[i].c);
        free(neutral->T0[i].c);
    }
    free(neutral->n0);
    free(neutral->T0);
}

void NeutralRadial_offload(NeutralRadial *neutral) { SUPPRESS_UNUSED_WARNING(neutral); }

err_t NeutralRadial_eval_n0(real *n0, real rho, NeutralRadial *neutral)
{
    err_t err = 0;
    int interperr = 0;
    for (size_t i = 0; i < neutral->n; i++)
    {
        interperr += linint1D_eval_f(&n0[i], &neutral->n0[i], rho);
    }

    if (interperr)
    {
        return error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_N0_1D);
    }

    return err;
}

err_t NeutralRadial_eval_T0(real *T0, real rho, NeutralRadial *neutral)
{
    err_t err = 0;
    int interperr = 0;
    for (size_t i = 0; i < neutral->n; i++)
    {
        interperr += linint1D_eval_f(&T0[i], &neutral->T0[i], rho);
    }

    if (interperr)
    {
        return error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_N0_1D);
    }

    return err;
}
