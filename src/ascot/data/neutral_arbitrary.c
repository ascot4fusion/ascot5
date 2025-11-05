/**
 * Implements neutral_arbitrary.h.
 */
#include "defines.h"
#include "neutral.h"
#include "utils/interp.h"
#include "utils/mathlib.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int NeutralArbitrary_init(
    NeutralArbitrary *neutral, size_t n, size_t nr, size_t nphi, size_t nz,
    real rlim[2], real philim[2], real zlim[2],
    real density[n * nr * nphi * nz], real temperature[n * nr * nphi * nz])
{
    neutral->n = n;
    neutral->density = (Linear3D *)malloc(n * sizeof(Linear3D));
    neutral->temperature = (Linear3D *)malloc(n * sizeof(Linear3D));
    int err = 0;
    for (size_t i = 0; i < neutral->n; i++)
    {
        real *c = (real *)malloc(nr * nphi * nz * sizeof(real));
        err += c == NULL ? 1 : 0;
        for (size_t i = 0; i < nr * nphi * nz; i++)
        {
            c[i] = density[i];
        }
        Linear3D_init(
            &neutral->density[i], nr, nphi, nz, NATURALBC, PERIODICBC,
            NATURALBC, rlim, philim, zlim, c);
        c = (real *)malloc(nr * nphi * nz * sizeof(real));
        err += c == NULL ? 1 : 0;
        for (size_t i = 0; i < nr * nphi * nz; i++)
        {
            c[i] = temperature[i];
        }
        Linear3D_init(
            &neutral->temperature[i], nr, nphi, nz, NATURALBC, PERIODICBC,
            NATURALBC, rlim, philim, zlim, c);
    }
    return err;
}

void NeutralArbitrary_free(NeutralArbitrary *neutral)
{
    for (size_t i = 0; i < neutral->n; i++)
    {
        free(neutral->density[i].c);
        free(neutral->temperature[i].c);
    }
    free(neutral->density);
    free(neutral->temperature);
}

void NeutralArbitrary_offload(NeutralArbitrary *neutral)
{
    SUPPRESS_UNUSED_WARNING(neutral);
}

err_t NeutralArbitrary_eval_density(
    real *density, real r, real phi, real z, NeutralArbitrary *neutral)
{
    err_t err = 0;
    int interperr = 0;
    for (size_t i = 0; i < neutral->n; i++)
    {
        interperr +=
            Linear3D_eval_f(&density[i], &neutral->density[i], r, phi, z);
    }

    err = ERROR_CHECK(
        err, interperr, ERR_INTERPOLATED_OUTSIDE_RANGE,
        DATA_NEUTRAL_ARBITRARY_C);
    return err;
}

err_t NeutralArbitrary_eval_temperature(
    real *temperature, real r, real phi, real z, NeutralArbitrary *neutral)
{
    err_t err = 0;
    int interperr = 0;
    for (size_t i = 0; i < neutral->n; i++)
    {
        interperr += Linear3D_eval_f(
            &temperature[i], &neutral->temperature[i], r, phi, z);
    }

    err = ERROR_CHECK(
        err, interperr, ERR_INTERPOLATED_OUTSIDE_RANGE,
        DATA_NEUTRAL_ARBITRARY_C);
    return err;
}
