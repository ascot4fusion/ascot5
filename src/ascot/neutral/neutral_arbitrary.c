/**
 * Implements N0_3D.h.
 */
#include "neutral.h"
#include "defines.h"
#include "errors.h"
#include "linint.h"
#include "mathlib.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int N0_3D_init(
    N0_3D_data *data, int n_r, real r_min, real r_max, int n_phi, real phi_min,
    real phi_max, int n_z, real z_min, real z_max, int n_species, int *anum,
    int *znum, real *density, real *temperature)
{
    data->n_species = n_species;
    data->anum = (int *)malloc(n_species * sizeof(int));
    data->znum = (int *)malloc(n_species * sizeof(int));
    data->n0 = (linint3D_data *)malloc(n_species * sizeof(linint3D_data));
    data->t0 = (linint3D_data *)malloc(n_species * sizeof(linint3D_data));
    int err = 0;
    for (int i = 0; i < data->n_species; i++)
    {
        data->anum[i] = anum[i];
        data->znum[i] = znum[i];

        real *c = (real *)malloc(n_r * n_phi * n_z * sizeof(real));
        err += c == NULL ? 1 : 0;
        for (int i = 0; i < n_r * n_phi * n_z; i++)
        {
            c[i] = density[i];
        }
        linint3D_init(
            &data->n0[i], c, n_r, n_phi, n_z, NATURALBC, PERIODICBC, NATURALBC,
            r_min, r_max, phi_min, phi_max, z_min, z_max);
        c = (real *)malloc(n_r * n_phi * n_z * sizeof(real));
        err += c == NULL ? 1 : 0;
        for (int i = 0; i < n_r * n_phi * n_z; i++)
        {
            c[i] = temperature[i];
        }
        linint3D_init(
            &data->t0[i], c, n_r, n_phi, n_z, NATURALBC, PERIODICBC, NATURALBC,
            r_min, r_max, phi_min, phi_max, z_min, z_max);
    }
    return err;
}

void N0_3D_free(N0_3D_data *data)
{
    free(data->anum);
    free(data->znum);
    for (int i = 0; i < data->n_species; i++)
    {
        free(data->n0[i].c);
        free(data->t0[i].c);
    }
    free(data->n0);
    free(data->t0);
}

void N0_3D_offload(N0_3D_data *data)
{
    SUPPRESS_UNUSED_WARNING(data);
}

a5err N0_3D_eval_n0(real *n0, real r, real phi, real z, N0_3D_data *ndata)
{
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    for (int i = 0; i < ndata->n_species; i++)
    {
        interperr += linint3D_eval_f(&n0[i], &ndata->n0[i], r, phi, z);
    }

    if (interperr)
    {
        return error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_N0_3D);
    }

    return err;
}

a5err N0_3D_eval_t0(real *t0, real r, real phi, real z, N0_3D_data *ndata)
{
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    for (int i = 0; i < ndata->n_species; i++)
    {
        interperr += linint3D_eval_f(&t0[i], &ndata->t0[i], r, phi, z);
    }

    if (interperr)
    {
        return error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_N0_3D);
    }

    return err;
}

int N0_3D_get_n_species(N0_3D_data *ndata) { return ndata->n_species; }
