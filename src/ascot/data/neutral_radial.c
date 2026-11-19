/**
 * Implements "neutral_radial.h".
 */
#include "defines.h"
#include "neutral.h"
#include "utils/interp.h"
#include "utils/mathlib.h"
#include <stdlib.h>
#include <string.h>

int NeutralRadial_init(
    NeutralRadial *neutral, size_t nspecies, size_t nrho, int anum[nspecies],
    int znum[nspecies], real rholim[2], real density[nspecies * nrho],
    real temperature[nspecies * nrho])
{
    neutral->nspecies = nspecies;

    neutral->anum = (int *)malloc(nspecies * sizeof(int));
    neutral->znum = (int *)malloc(nspecies * sizeof(int));
    neutral->density = (Linear1D *)malloc(nspecies * sizeof(Linear1D));
    neutral->temperature = (Linear1D *)malloc(nspecies * sizeof(Linear1D));

    int err = 0;
    for (size_t i = 0; i < neutral->nspecies; i++)
    {
        neutral->znum[i] = znum[i];
        neutral->anum[i] = anum[i];

        real *c = (real *)malloc(nrho * sizeof(real));
        err += c == NULL ? 1 : 0;
        for (size_t j = 0; j < nrho; j++)
            c[j] = density[i * nrho + j];

        Linear1D_init(&neutral->density[i], nrho, NATURALBC, rholim, c);
        c = (real *)malloc(nrho * sizeof(real));
        err += c == NULL ? 1 : 0;
        for (size_t j = 0; j < nrho; j++)
            c[j] = temperature[i * nrho + j];

        Linear1D_init(&neutral->temperature[i], nrho, NATURALBC, rholim, c);
    }
    return err;
}

void NeutralRadial_free(NeutralRadial *neutral)
{
    for (size_t i = 0; i < neutral->nspecies; i++)
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
    real density[1], real rho, size_t i_species, NeutralRadial *neutral)
{
    err_t err = 0;
    int interperr = Linear1D_eval_f(density, &neutral->density[i_species], rho);

    err = ERROR_CHECK(
        err, interperr, ERR_INTERPOLATED_OUTSIDE_RANGE, DATA_NEUTRAL_RADIAL_C);
    return err;
}

err_t NeutralRadial_eval_temperature(
    real temperature[1], real rho, size_t i_species, NeutralRadial *neutral)
{
    err_t err = 0;
    int interperr =
        Linear1D_eval_f(temperature, &neutral->temperature[i_species], rho);

    err = ERROR_CHECK(
        err, interperr, ERR_INTERPOLATED_OUTSIDE_RANGE, DATA_NEUTRAL_RADIAL_C);
    return err;
}

err_t NeutralRadial_eval_nT(
    real *density, real *temperature, real rho, NeutralRadial *neutral)
{
    err_t err = 0;
    int interperr = 0;
    for (size_t i = 0; i < neutral->nspecies; i++)
    {
        interperr += Linear1D_eval_f(&density[i], &neutral->density[i], rho);
        interperr +=
            Linear1D_eval_f(&temperature[i], &neutral->temperature[i], rho);
    }

    err = ERROR_CHECK(
        err, interperr, ERR_INTERPOLATED_OUTSIDE_RANGE, DATA_NEUTRAL_RADIAL_C);
    return err;
}
