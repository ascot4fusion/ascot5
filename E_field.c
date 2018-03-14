/**
 * @file E_field.c
 * @brief Electric field interface
 */
#include <stdio.h>
#include "ascot5.h"
#include "error.h"
#include "E_field.h"
#include "B_field.h"
#include "Efield/E_TC.h"
#include "Efield/E_1D.h"
#include "Efield/E_1DS.h"
#include "Efield/E_3D.h"


void E_field_init_offload(E_field_offload_data* offload_data,
                          real** offload_array) {
  /*offload_data->type = E_field_type_1D;*/

    switch(offload_data->type) {
    case E_field_type_1D:
        E_1D_init_offload(&(offload_data->E1D), offload_array);
        offload_data->offload_array_length = offload_data->E1D.offload_array_length;
        break;
    case E_field_type_1DS:
	// Do nothing
        break;
    case E_field_type_TC:
	// Do nothing
        break;
    case E_field_type_3D:
        E_3D_init_offload(&(offload_data->E3D), offload_array);
        offload_data->offload_array_length = offload_data->E3D.offload_array_length;
        break;
    }
}

void E_field_free_offload(E_field_offload_data* offload_data,
                          real** offload_array) {
    switch(offload_data->type) {
    case E_field_type_1D:
        E_1D_free_offload(&(offload_data->E1D), offload_array);
        break;
    case E_field_type_1DS:
        E_1DS_free_offload(&(offload_data->E1DS), offload_array);
        break;
    case E_field_type_TC:
        E_TC_free_offload(&(offload_data->ETC), offload_array);
        break;
    case E_field_type_3D:
      E_3D_free_offload(&(offload_data->E3D), offload_array);
      break;

    }
}

int E_field_init(E_field_data* Edata, E_field_offload_data* offload_data,
                  real* offload_array) {
    int err = 0;

    switch(offload_data->type) {
    case E_field_type_1D:
        E_1D_init(&(Edata->E1D), &(offload_data->E1D), offload_array);
        break;
    case E_field_type_1DS:
        E_1DS_init(&(Edata->E1DS), &(offload_data->E1DS), offload_array);
        break;
    case E_field_type_TC:
        E_TC_init(&(Edata->ETC), &(offload_data->ETC), offload_array);
        break;
    case E_field_type_3D:
        E_3D_init(&(Edata->E3D), &(offload_data->E3D), offload_array);
        break;
    }
    Edata->type = offload_data->type;

    return err;
}

/* E[0] = E_r, E[1] = E_phi, E[2] = E_z */
a5err E_field_eval_E(real E[], real r, real phi, real z, E_field_data* Edata, B_field_data* Bdata) {
    a5err err = 0;

    switch(Edata->type) {
    case E_field_type_1D:
        E_1D_eval_E(E, r, phi, z, &(Edata->E1D), Bdata);
        break;
    case E_field_type_1DS:
        E_1DS_eval_E(E, r, phi, z, &(Edata->E1DS), Bdata);
        break;
    case E_field_type_TC:
        E_TC_eval_E(E, r, phi, z, &(Edata->ETC), Bdata);
        break;
    case E_field_type_3D:
        E_3D_eval_E(E, r, phi, z, &(Edata->E3D), Bdata);
        break;
    }

    return err;
}

a5err E_field_eval_E_SIMD(int i, real E[3][NSIMD], real r, real phi, real z, E_field_data* Edata, B_field_data* Bdata) {
    a5err err = 0;

    switch(Edata->type) {
    case E_field_type_1D:
        E_1D_eval_E_SIMD(i, E, r, phi, z, &(Edata->E1D), Bdata);
        break;
    case E_field_type_1DS:
        E_1DS_eval_E_SIMD(i, E, r, phi, z, &(Edata->E1DS), Bdata);
        break;
    case E_field_type_TC:
        E_TC_eval_E_SIMD(i, E, r, phi, z, &(Edata->ETC), Bdata);
        break;
    case E_field_type_3D:
        E_3D_eval_E_SIMD(i, E, r, phi, z, &(Edata->E3D), Bdata);
        break;
    }

    return err;
}
