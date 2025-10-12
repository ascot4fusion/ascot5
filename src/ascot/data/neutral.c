/**
 * Implements neutral.h.
 */
#include "neutral.h"
#include "defines.h"
#include "neutral_radial.h"
#include "neutral_arbitrary.h"
#include <stdio.h>

void Neutral_free(Neutral *data)
{
    switch (data->type)
    {
    case NEUTRAL_RADIAL:
        NeutralRadial_free(data->radial);
        break;
    case NEUTRAL_ARBITRARY:
        NeutralArbitrary_free(data->arbitrary);
        break;
    }
}

void Neutral_offload(Neutral *data)
{
    switch (data->type)
    {
    case NEUTRAL_RADIAL:
        NeutralRadial_offload(data->radial);
        break;
    case NEUTRAL_ARBITRARY:
        NeutralArbitrary_offload(data->arbitrary);
        break;
    }
}

err_t Neutral_eval_n0(
    real *n0, real rho, real r, real phi, real z, real t, Neutral *neutral)
{
    (void)t; // Unused until dynamic neutral is implemented.
    size_t n = 0;
    err_t err = 0;
    switch (neutral->type)
    {
    case NEUTRAL_RADIAL:
        err = NeutralRadial_eval_n0(n0, rho, neutral->radial);
        n = neutral->radial->n;
        break;
    case NEUTRAL_ARBITRARY:
        err = NeutralArbitrary_eval_n0(n0, r, phi, z, neutral->arbitrary);
        n = neutral->arbitrary->n;
        break;
    }

    for(size_t i = 0; i < n; i++)
        n0[i] = err ? 1.0 : n0[i];

    return err;
}

err_t Neutral_eval_T0(
    real *T0, real rho, real r, real phi, real z, real t, Neutral *neutral)
{
    (void)t; // Unused until dynamic neutral is implemented.
    size_t n = 0;
    err_t err = 0;
    switch (neutral->type)
    {
    case NEUTRAL_RADIAL:
        err = NeutralRadial_eval_T0(T0, rho, neutral->radial);
        n = neutral->radial->n;
        break;
    case NEUTRAL_ARBITRARY:
        err = NeutralArbitrary_eval_T0(T0, r, phi, z, neutral->arbitrary);
        n = neutral->arbitrary->n;
        break;
    }

    for(size_t i = 0; i < n; i++)
        T0[i] = err ? 1.0 : T0[i];

    return err;
}

size_t Neutral_get_n_species(Neutral *neutral)
{
    size_t n = 0;
    switch (neutral->type)
    {
    case NEUTRAL_RADIAL:
        n = neutral->radial->n;
        break;
    case NEUTRAL_ARBITRARY:
        n = neutral->arbitrary->n;
        break;
    }

    return n;
}
