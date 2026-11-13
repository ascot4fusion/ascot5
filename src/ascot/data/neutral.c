/**
 * Implements "neutral.h".
 */
#include "neutral.h"
#include "defines.h"
#include "neutral_arbitrary.h"
#include "neutral_radial.h"

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
    real density[1], real rho, real r, real phi, real z, real t,
    size_t i_species, Neutral *neutral)
{
    (void)t; // Unused until dynamic neutral is implemented.
    err_t err = 0;
    switch (neutral->type)
    {
    case NEUTRAL_RADIAL:
        err = NeutralRadial_eval_density(
            density, rho, i_species, neutral->radial);
        break;
    case NEUTRAL_ARBITRARY:
        err = NeutralArbitrary_eval_density(
            density, r, phi, z, i_species, neutral->arbitrary);
        break;
    }

    density[0] = err ? 1.0 : density[0];

    return err;
}

err_t Neutral_eval_temperature(
    real temperature[1], real rho, real r, real phi, real z, real t,
    size_t i_species, Neutral *neutral)
{
    (void)t; // Unused until dynamic neutral is implemented.
    err_t err = 0;
    switch (neutral->type)
    {
    case NEUTRAL_RADIAL:
        err = NeutralRadial_eval_temperature(
            temperature, rho, i_species, neutral->radial);
        break;
    case NEUTRAL_ARBITRARY:
        err = NeutralArbitrary_eval_temperature(
            temperature, r, phi, z, i_species, neutral->arbitrary);
        break;
    }

    temperature[0] = err ? 1.0 : temperature[0];

    return err;
}

err_t Neutral_eval_nT(
    real *density, real *temperature, real rho, real r, real phi, real z,
    real t, Neutral *neutral)
{
    (void)t; // Unused until dynamic neutral is implemented.
    size_t n = 0;
    err_t err = 0;
    switch (neutral->type)
    {
    case NEUTRAL_RADIAL:
        err = NeutralRadial_eval_nT(density, temperature, rho, neutral->radial);
        n = neutral->radial->nspecies;
        break;
    case NEUTRAL_ARBITRARY:
        err = NeutralArbitrary_eval_nT(
            density, temperature, r, phi, z, neutral->arbitrary);
        n = neutral->arbitrary->nspecies;
        break;
    }

    for (size_t i = 0; i < n; i++)
    {
        density[i] = err ? 1.0 : density[i];
        temperature[i] = err ? 1.0 : temperature[i];
    }

    return err;
}

size_t Neutral_get_n_species(Neutral *neutral)
{
    size_t n = 0;
    switch (neutral->type)
    {
    case NEUTRAL_RADIAL:
        n = neutral->radial->nspecies;
        break;
    case NEUTRAL_ARBITRARY:
        n = neutral->arbitrary->nspecies;
        break;
    }

    return n;
}

const int *Neutral_get_species_znum(Neutral *neutral)
{
    const int *znum = NULL;
    switch (neutral->type)
    {
    case NEUTRAL_RADIAL:
        znum = neutral->radial->znum;
        break;
    case NEUTRAL_ARBITRARY:
        znum = neutral->arbitrary->znum;
        break;
    }

    return znum;
}

const int *Neutral_get_species_anum(Neutral *neutral)
{
    const int *anum = NULL;
    switch (neutral->type)
    {
    case NEUTRAL_RADIAL:
        anum = neutral->radial->anum;
        break;
    case NEUTRAL_ARBITRARY:
        anum = neutral->arbitrary->anum;
        break;
    }

    return anum;
}
