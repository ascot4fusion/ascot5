/**
 * @file E_field.c
 * @brief Electric field interface
 */
#include <stdio.h>
#include "ascot5.h"
#include "E_field.h"
#include "B_field.h"
#include "E_1D.h"

void E_field_init_offload(E_field_offload_data* offload_data,
                          real** offload_array) {
    offload_data->type = E_field_type_1D;

    switch(offload_data->type) {
    case E_field_type_1D:
        E_1D_init_offload(&(offload_data->E1D), offload_array);
        offload_data->offload_array_length = offload_data->E1D.offload_array_length;
        break;
    }
}

void E_field_free_offload(E_field_offload_data* offload_data,
                          real** offload_array) {
    switch(offload_data->type) {
    case E_field_type_1D:
        E_1D_free_offload(&(offload_data->E1D), offload_array);
        break;
    }
}

void E_field_init(E_field_data* Edata, E_field_offload_data* offload_data,
                  real* offload_array) {
    switch(offload_data->type) {
    case E_field_type_1D:
        E_1D_init(&(Edata->E1D), &(offload_data->E1D), offload_array);
        break;
    }
    Edata->type = offload_data->type;
}

/* E[0] = E_r, E[1] = E_phi, E[2] = E_z */
void E_field_eval_E(real E[], real r, real phi, real z, E_field_data* Edata, B_field_data* Bdata) {
    switch(Edata->type) {
    case E_field_type_1D:
	
        E_1D_eval_E(E, r, phi, z, &(Edata->E1D), Bdata);
        break;
    }
}
