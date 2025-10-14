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
    NeutralRadial *neutral, size_t n, size_t nrho, real rholim[2],
    real density[n * nrho], real temperature[n * nrho])
{
    neutral->n = n;
    neutral->density = (Linear1D *)malloc(n * sizeof(Linear1D));
    neutral->temperature = (Linear1D *)malloc(n * sizeof(Linear1D));
    int err = 0;
    for (size_t i = 0; i < neutral->n; i++)
    {
        real *c = (real *)malloc(nrho * sizeof(real));
        err += c == NULL ? 1 : 0;
        for (size_t i = 0; i < nrho; i++)
        {
            c[i] = density[i];
        }
        linint1D_init(
            &neutral->density[i], c, nrho, NATURALBC, rholim[0], rholim[1]);
        c = (real *)malloc(nrho * sizeof(real));
        err += c == NULL ? 1 : 0;
        for (size_t i = 0; i < nrho; i++)
        {
            c[i] = temperature[i];
        }
        linint1D_init(
            &neutral->temperature[i], c, nrho, NATURALBC, rholim[0], rholim[1]);
    }
    return err;
}

void NeutralRadial_free(NeutralRadial *neutral)
{
    for (size_t i = 0; i < neutral->n; i++)
    {
        free(neutral->density[i].c);
        free(neutral->temperature[i].c);
    }
    free(neutral->density);
    free(neutral->temperature);
}

void NeutralRadial_offload(NeutralRadial *neutral)
{
    SUPPRESS_UNUSED_WARNING(neutral);
}

err_t NeutralRadial_eval_density(
    real *density, real rho, NeutralRadial *neutral)
{
    err_t err = 0;
    int interperr = 0;
    for (size_t i = 0; i < neutral->n; i++)
    {
        interperr += linint1D_eval_f(&density[i], &neutral->density[i], rho);
    }

    err = ERROR_CHECK(
        err, interperr, ERR_INTERPOLATED_OUTSIDE_RANGE, DATA_NEUTRAL_RADIAL_C);
    return err;
}

err_t NeutralRadial_eval_temperature(
    real *temperature, real rho, NeutralRadial *neutral)
{
    err_t err = 0;
    int interperr = 0;
    for (size_t i = 0; i < neutral->n; i++)
    {
        interperr +=
            linint1D_eval_f(&temperature[i], &neutral->temperature[i], rho);
    }

    err = ERROR_CHECK(
        err, interperr, ERR_INTERPOLATED_OUTSIDE_RANGE, DATA_NEUTRAL_RADIAL_C);
    return err;
}
