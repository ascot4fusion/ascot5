/**
 * Implements neutral.h.
 */
#include "neutral.h"
#include "ascot5.h"
#include "error.h"
#include "neutral/N0_1D.h"
#include "neutral/N0_3D.h"
#include <stdio.h>

void neutral_free(neutral_data *data)
{
    switch (data->type)
    {
    case neutral_type_1D:
        N0_1D_free(data->N01D);
        break;
    case neutral_type_3D:
        N0_3D_free(data->N03D);
        break;
    }
}

void neutral_offload(neutral_data *data)
{
    switch (data->type)
    {
    case neutral_type_1D:
        N0_1D_offload(data->N01D);
        break;
    case neutral_type_3D:
        N0_3D_offload(data->N03D);
        break;
    }
}

a5err neutral_eval_n0(
    real *n0, real rho, real r, real phi, real z, real t, neutral_data *ndata)
{
    a5err err = 0;

    switch (ndata->type)
    {
    case neutral_type_1D:
        err = N0_1D_eval_n0(n0, rho, ndata->N01D);
        break;
    case neutral_type_3D:
        err = N0_3D_eval_n0(n0, r, phi, z, ndata->N03D);
        break;
    default:
        /* Unregonized input. Produce error. */
        err = error_raise(ERR_UNKNOWN_INPUT, __LINE__, EF_NEUTRAL);
        break;
    }

    if (err)
    {
        /* Return some reasonable values to avoid further errors */
        n0[0] = 0;
    }

    return err;
}

a5err neutral_eval_t0(
    real *t0, real rho, real r, real phi, real z, real t, neutral_data *ndata)
{
    a5err err = 0;

    switch (ndata->type)
    {
    case neutral_type_1D:
        err = N0_1D_eval_t0(t0, rho, ndata->N01D);
        break;
    case neutral_type_3D:
        err = N0_3D_eval_t0(t0, r, phi, z, ndata->N03D);
        break;
    default:
        /* Unregonized input. Produce error. */
        err = error_raise(ERR_UNKNOWN_INPUT, __LINE__, EF_NEUTRAL);
        break;
    }

    if (err)
    {
        /* Return some reasonable values to avoid further errors */
        t0[0] = 1;
    }

    return err;
}

int neutral_get_n_species(neutral_data *ndata)
{
    int n = 0;
    switch (ndata->type)
    {
    case neutral_type_1D:
        n = N0_1D_get_n_species(ndata->N01D);
        break;
    case neutral_type_3D:
        n = N0_3D_get_n_species(ndata->N03D);
        break;
    }

    return n;
}
