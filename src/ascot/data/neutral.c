/**
 * Implements neutral.h.
 */
#include "neutral.h"
#include "defines.h"
#include "neutral_arbitrary.h"
#include "neutral_radial.h"
#include <stdio.h>

void Neutral_free(Neutral *neutral)
{
    switch (neutral->type)
    {
    case NEUTRAL_RADIAL:
        NeutralRadial_free(neutral->radial);
        break;
    case NEUTRAL_ARBITRARY:
        NeutralArbitrary_free(neutral->arbitrary);
        break;
    }
}

void Neutral_offload(Neutral *neutral)
{
    switch (neutral->type)
    {
    case NEUTRAL_RADIAL:
        NeutralRadial_offload(neutral->radial);
        break;
    case NEUTRAL_ARBITRARY:
        NeutralArbitrary_offload(neutral->arbitrary);
        break;
    }
}

err_t Neutral_eval_density(
    real *density, real rho, real r, real phi, real z, real t, Neutral *neutral)
{
    (void)t; // Unused until dynamic neutral is implemented.
    size_t n = 0;
    err_t err = 0;
    switch (neutral->type)
    {
    case NEUTRAL_RADIAL:
        err = NeutralRadial_eval_density(density, rho, neutral->radial);
        n = neutral->radial->n;
        break;
    case NEUTRAL_ARBITRARY:
        err = NeutralArbitrary_eval_density(
            density, r, phi, z, neutral->arbitrary);
        n = neutral->arbitrary->n;
        break;
    }

    for (size_t i = 0; i < n; i++)
        density[i] = err ? 1.0 : density[i];

    return err;
}

err_t Neutral_eval_temperature(
    real *temperature, real rho, real r, real phi, real z, real t,
    Neutral *neutral)
{
    (void)t; // Unused until dynamic neutral is implemented.
    size_t n = 0;
    err_t err = 0;
    switch (neutral->type)
    {
    case NEUTRAL_RADIAL:
        err = NeutralRadial_eval_temperature(temperature, rho, neutral->radial);
        n = neutral->radial->n;
        break;
    case NEUTRAL_ARBITRARY:
        err = NeutralArbitrary_eval_temperature(
            temperature, r, phi, z, neutral->arbitrary);
        n = neutral->arbitrary->n;
        break;
    }

    for (size_t i = 0; i < n; i++)
        temperature[i] = err ? 1.0 : temperature[i];

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
