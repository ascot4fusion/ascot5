/**
 * @file E_2DS.c
 * @brief 2D spline electric field
 */
#include <stdlib.h>
#include "../ascot5.h"
#include "../offload.h"
#include "../print.h"
#include "../error.h"
#include "../spline/interp.h"
#include "E_2DS.h"

int E_2DS_init(E_2DS_data* data, int nr, real rmin, real rmax, int nz,
    real zmin, real zmax, real* vpot) {

    int err = 0;

    /* Set up the splines */
    err = interp2Dcomp_setup(&data->vpot, vpot, nr, nz, NATURALBC, NATURALBC,
                             rmin, rmax, zmin, zmax);

    if(err) {
        print_err("Error: Failed to initialize splines.\n");
        return 1;
    }

    return 0;
}


void E_2DS_free(E_2DS_data* data) {
    free(data->vpot.c);
}


void E_2DS_offload(E_2DS_data* data) {
    GPU_MAP_TO_DEVICE(
        data->vpot, data->vpot.c[0:data->vpot.n_x*data->vpot.n_y*NSIZE_COMP2D]
    )
}


a5err E_2DS_eval_E(real E[3], real r, real z, E_2DS_data* Edata) {
    real vdv[6];
    int interperr = 0;
    interperr += interp2Dcomp_eval_df(vdv, &Edata->vpot, r, z);

    a5err err = 0;
    if(interperr) {
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_E_2DS );
    }

    E[0] = -vdv[1];
    E[1] = 0;
    E[2] = -vdv[2];

    return err;
}
