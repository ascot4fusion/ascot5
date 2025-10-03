/**
 * Implements plasma.h.
 */
#include "plasma.h"
#include "ascot5.h"
#include "consts.h"
#include "error.h"
#include "plasma/plasma_1D.h"
#include "plasma/plasma_1Dt.h"
#include <stdio.h>
#include <stdlib.h>

void plasma_free(plasma_data *plasma)
{
    switch (plasma->type)
    {
    case plasma_type_1D:
        PlasmaLinear1D_free(plasma->linear1d);
        break;
    case plasma_type_1Dt:
        PlasmaDynamic1D_free(plasma->dynamic1d);
        break;
    }
}

void plasma_offload(plasma_data *plasma)
{
    switch (plasma->type)
    {
    case plasma_type_1D:
        PlasmaLinear1D_offload(plasma->linear1d);
        break;
    case plasma_type_1Dt:
        PlasmaDynamic1D_offload(plasma->dynamic1d);
        break;
    }
}

a5err plasma_eval_temp(
    real *temp, real rho, real r, real phi, real z, real t, int species,
    plasma_data *plasma)
{
    a5err err = 0;

    switch (plasma->type)
    {
    case plasma_type_1D:
        err = PlasmaLinear1D_eval_temp(temp, rho, species, plasma->linear1d);
        break;
    case plasma_type_1Dt:
        err =
            PlasmaDynamic1D_eval_temp(temp, rho, t, species, plasma->dynamic1d);
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

a5err plasma_eval_dens(
    real *dens, real rho, real r, real phi, real z, real t, int species,
    plasma_data *plasma)
{
    a5err err = 0;

    switch (plasma->type)
    {
    case plasma_type_1D:
        err = PlasmaLinear1D_eval_dens(dens, rho, species, plasma->linear1d);
        break;

    case plasma_type_1Dt:
        err =
            PlasmaDynamic1D_eval_dens(dens, rho, t, species, plasma->dynamic1d);
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

a5err plasma_eval_densandtemp(
    real *dens, real *temp, real rho, real r, real phi, real z, real t,
    plasma_data *plasma)
{
    a5err err = 0;

    switch (plasma->type)
    {
    case plasma_type_1D:
        err =
            PlasmaLinear1D_eval_densandtemp(dens, temp, rho, plasma->linear1d);
        break;
    case plasma_type_1Dt:
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

a5err plasma_eval_flow(
    real *vflow, real rho, real r, real phi, real z, real t,
    plasma_data *plasma)
{
    a5err err = 0;

    switch (plasma->type)
    {
    case plasma_type_1D:
        err = PlasmaLinear1D_eval_flow(vflow, rho, r, plasma->linear1d);
        break;
    case plasma_type_1Dt:
        err = PlasmaDynamic1D_eval_flow(vflow, rho, t, r, plasma->dynamic1d);
        break;
    }

    if (err)
    {
        /* In case of error, return some reasonable values to avoid further
           complications */
        *vflow = 0;
    }

    return err;
}

int plasma_get_n_species(plasma_data *plasma)
{
    int n = 0;
    switch (plasma->type)
    {
    case plasma_type_1D:
        n = plasma->linear1d->nspecies;
        break;
    case plasma_type_1Dt:
        n = plasma->dynamic1d->nspecies;
        break;
    }

    return n;
}

const real *plasma_get_species_mass(plasma_data *plasma)
{
    const real *mass = NULL;
    switch (plasma->type)
    {
    case plasma_type_1D:
        mass = plasma->linear1d->mass;
        break;
    case plasma_type_1Dt:
        mass = plasma->dynamic1d->mass;
        break;
    }

    return mass;
}

const real *plasma_get_species_charge(plasma_data *plasma)
{
    const real *charge = NULL;
    switch (plasma->type)
    {
    case plasma_type_1D:
        charge = plasma->linear1d->charge;
        break;
    case plasma_type_1Dt:
        charge = plasma->dynamic1d->charge;
        break;
    }

    return charge;
}

const int *plasma_get_species_znum(plasma_data *plasma)
{
    const int *znum = NULL;
    switch (plasma->type)
    {
    case plasma_type_1D:
        znum = plasma->linear1d->znum;
        break;
    case plasma_type_1Dt:
        znum = plasma->dynamic1d->znum;
        break;
    }

    return znum;
}

const int *plasma_get_species_anum(plasma_data *plasma)
{
    const int *anum = NULL;
    switch (plasma->type)
    {
    case plasma_type_1D:
        anum = plasma->linear1d->anum;
        break;
    case plasma_type_1Dt:
        anum = plasma->dynamic1d->anum;
        break;
    }

    return anum;
}
