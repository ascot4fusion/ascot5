/**
 * @authot Konsta Sarkimaki konsta.sarkimaki@aalto.fi
 * @file E_TC.c
 * @brief Trivial Cartesian Electric field.
 *
 * Electric field that has constant x, y, and z components.
 * Intended for testing purposes.
 */
#include <stdlib.h>
#include "../ascot5.h"
#include "../math.h"
#include "../B_field.h"
#include "E_TC.h"


void E_TC_init_offload(E_TC_offload_data* offload_data,
			     real** offload_array) {
    // redundant function
}

void E_TC_free_offload(E_TC_offload_data* offload_data,
			     real** offload_array) {
    free(*offload_array);
    *offload_array = NULL;
}

void E_TC_init(E_TC_data* Edata, E_TC_offload_data* offload_data,
		     real* offload_array) {
    Edata->Exyz = offload_array;
}
 
void E_TC_eval_E(real* E, real r, real phi, real z, E_TC_data* Edata, B_field_data* Bdata) {
    math_vec_xyz2rpz(Edata->Exyz, E, phi);
}

void E_TC_eval_E_SIMD(int i, real E[3][NSIMD], real r, real phi, real z, E_TC_data* Edata, B_field_data* Bdata) {
    real E0[3];
    math_vec_xyz2rpz(Edata->Exyz, E0, phi);
    E[0][i] = E0[0];
    E[1][i] = E0[1];
    E[2][i] = E0[2];
}
