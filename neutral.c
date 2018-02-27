/**
 * @file neutral.c
 * @brief Neutral density interface
 */
#include <stdio.h>
#include "ascot5.h"
#include "error.h"
#include "neutral.h"
#include "neutral/N0_3D.h"
#include "neutral/N0_ST.h"

void neutral_init_offload(neutral_offload_data* offload_data,
                          real** offload_array) {
    
    switch(offload_data->type) {
        case neutral_type_3D:
        N0_3D_init_offload(&(offload_data->N03D), offload_array);
        offload_data->offload_array_length = offload_data->N03D.offload_array_length;
        break;

        case neutral_type_ST:
        N0_ST_init_offload(&(offload_data->N0ST), offload_array);
        offload_data->offload_array_length = offload_data->N0ST.offload_array_length;
        break;
    }
}

void neutral_free_offload(neutral_offload_data* offload_data,
                          real** offload_array) {
    switch(offload_data->type) {
        case neutral_type_3D:
        N0_3D_free_offload(&(offload_data->N03D), offload_array);
        break;

        case neutral_type_ST:
        N0_ST_free_offload(&(offload_data->N0ST), offload_array);
        break;
    }
}

int neutral_init(neutral_data* ndata, neutral_offload_data* offload_data,
                  real* offload_array) {
    int err = 0;

    switch(offload_data->type) {
        case neutral_type_3D:
        err = N0_3D_init(&(ndata->N03D), &(offload_data->N03D), offload_array);
        break;

        case neutral_type_ST:
        err = N0_ST_init(&(ndata->N0ST), &(offload_data->N0ST), offload_array);
        break;
    }
    ndata->type = offload_data->type;

    return err;
}

a5err neutral_eval_n0(real n0[], real r, real phi, real z,
                      neutral_data* ndata) {
    a5err err = 0;

    switch(ndata->type) {
        case neutral_type_3D:
        err = N0_3D_eval_n0(n0, r, phi, z, &(ndata->N03D));
        break;

        case neutral_type_ST:
        err = N0_ST_eval_n0(n0, r, phi, z, &(ndata->N0ST));
        break;
    }

    if(err) {
	/* Return some reasonable values to avoid further errors */
	n0[0] = 0;
    }

    return err;
}


a5err neutral_eval_n0_SIMD(int i, real n0[NSIMD], real r, real phi, real z,
                      neutral_data* ndata) {
    a5err err = 0;

    switch(ndata->type) {
        case neutral_type_3D:
        N0_3D_eval_n0_SIMD(i, n0, r, phi, z, &(ndata->N03D));
        break;

        case neutral_type_ST:
        N0_ST_eval_n0_SIMD(i, n0, r, phi, z, &(ndata->N0ST));
        break;
    }

    if(err) {
	/* Return some reasonable values to avoid further errors */
	n0[i] = 0;
    }

    return err;
}

