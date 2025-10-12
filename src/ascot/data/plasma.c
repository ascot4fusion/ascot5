/**
 * Implements plasma.h.
 */
#include "plasma.h"
#include "consts.h"
#include "defines.h"
#include "plasma_dynamic1d.h"
#include "plasma_linear1d.h"
#include <stdio.h>
#include <stdlib.h>

void Plasma_free(Plasma *plasma)
{
    switch (plasma->type)
    {
    case PLASMA_LINEAR1D:
        PlasmaLinear1D_free(plasma->linear1d);
        break;
    case PLASMA_DYNAMIC1D:
        PlasmaDynamic1D_free(plasma->dynamic1d);
        break;
    }
}

void Plasma_offload(Plasma *plasma)
{
    switch (plasma->type)
    {
    case PLASMA_LINEAR1D:
        PlasmaLinear1D_offload(plasma->linear1d);
        break;
    case PLASMA_DYNAMIC1D:
        PlasmaDynamic1D_offload(plasma->dynamic1d);
        break;
    }
}

err_t Plasma_eval_temp(
    real temperature[1], real rho, real r, real phi, real z, real t,
    size_t i_species, Plasma *plasma)
{
    /* Unused until 3D plasma is implemented */
    (void)r;
    (void)phi;
    (void)z;
    err_t err = 0;
    switch (plasma->type)
    {
    case PLASMA_LINEAR1D:
        err = PlasmaLinear1D_eval_temp(
            temperature, rho, i_species, plasma->linear1d);
        break;
    case PLASMA_DYNAMIC1D:
        err = PlasmaDynamic1D_eval_temp(
            temperature, rho, t, i_species, plasma->dynamic1d);
        break;
    }
    if (err)
    {
        /* In case of error, return some reasonable value to avoid further
           complications */
        temp[0] = 1e20;
    }

    return err;
}

err_t Plasma_eval_dens(
    real density[1], real rho, real r, real phi, real z, real t,
    size_t i_species, Plasma *plasma)
{
    /* Unused until 3D plasma is implemented */
    (void)r;
    (void)phi;
    (void)z;
    err_t err = 0;
    switch (plasma->type)
    {
    case PLASMA_LINEAR1D:
        err =
            PlasmaLinear1D_eval_dens(density, rho, i_species, plasma->linear1d);
        break;
    case PLASMA_DYNAMIC1D:
        err = PlasmaDynamic1D_eval_dens(
            density, rho, t, i_species, plasma->dynamic1d);
        break;
    }

    if (err)
    {
        /* In case of error, return some reasonable value to avoid further
           complications */
        dens[0] = 1e20;
    }
    return err;
}

err_t Plasma_eval_densandtemp(
    real *dens, real *temp, real rho, real r, real phi, real z, real t,
    Plasma *plasma)
{
    err_t err = 0;

    /* Unused until 3D plasma is implemented */
    (void)r;
    (void)phi;
    (void)z;

    switch (plasma->type)
    {
    case PLASMA_LINEAR1D:
        err =
            PlasmaLinear1D_eval_densandtemp(dens, temp, rho, plasma->linear1d);
        break;
    case PLASMA_DYNAMIC1D:
        err = PlasmaDynamic1D_eval_densandtemp(
            dens, temp, rho, t, plasma->dynamic1d);
        break;
    }

    if (err)
    {
        /* In case of error, return some reasonable values to avoid further
           complications */
        for (int i = 0; i < MAX_SPECIES; i++)
        {
            dens[i] = 1e20;
            temp[i] = 1e3;
        }
    }

    return err;
}

err_t Plasma_eval_flow(
    real vflow[1], real rho, real r, real phi, real z, real t, Plasma *plasma)
{
    /* Unused until 3D plasma is implemented */
    (void)phi;
    (void)z;
    err_t err = 0;
    switch (plasma->type)
    {
    case PLASMA_LINEAR1D:
        err = PlasmaLinear1D_eval_flow(vflow, rho, r, plasma->linear1d);
        break;
    case PLASMA_DYNAMIC1D:
        err = PlasmaDynamic1D_eval_flow(vflow, rho, t, r, plasma->dynamic1d);
        break;
    }

        vflow[0] = err ? 0.0 : vflow[0];
    return err;
}

size_t Plasma_get_n_species(Plasma *plasma)
{
    size_t nspecies = 0;
    switch (plasma->type)
    {
    case PLASMA_LINEAR1D:
        nspecies = plasma->linear1d->nspecies;
        break;
    case PLASMA_DYNAMIC1D:
        nspecies = plasma->dynamic1d->nspecies;
        break;
    }

    return nspecies;
}

const real *Plasma_get_species_mass(Plasma *plasma)
{
    const real *mass = NULL;
    switch (plasma->type)
    {
    case PLASMA_LINEAR1D:
        mass = plasma->linear1d->mass;
        break;
    case PLASMA_DYNAMIC1D:
        mass = plasma->dynamic1d->mass;
        break;
    }

    return mass;
}

const real *Plasma_get_species_charge(Plasma *plasma)
{
    const real *charge = NULL;
    switch (plasma->type)
    {
    case PLASMA_LINEAR1D:
        charge = plasma->linear1d->charge;
        break;
    case PLASMA_DYNAMIC1D:
        charge = plasma->dynamic1d->charge;
        break;
    }

    return charge;
}

const int *Plasma_get_species_znum(Plasma *plasma)
{
    const int *znum = NULL;
    switch (plasma->type)
    {
    case PLASMA_LINEAR1D:
        znum = plasma->linear1d->znum;
        break;
    case PLASMA_DYNAMIC1D:
        znum = plasma->dynamic1d->znum;
        break;
    }

    return znum;
}

const int *Plasma_get_species_anum(Plasma *plasma)
{
    const int *anum = NULL;
    switch (plasma->type)
    {
    case PLASMA_LINEAR1D:
        anum = plasma->linear1d->anum;
        break;
    case PLASMA_DYNAMIC1D:
        anum = plasma->dynamic1d->anum;
        break;
    }

    return anum;
}
