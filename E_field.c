/**
 * @file E_field.c
 * @brief Electric field interface
 */
#include <stdio.h>
#include "ascot5.h"
#include "E_field.h"

void E_field_init_offload(E_field_offload_data* offload_data,
                          real** offload_array) {
    FILE* f = fopen("input.efield", "r");
}

void E_field_free_offload(E_field_offload_data* offload_data,
                          real** offload_array) {
}

void E_field_init(E_field_data* Edata, E_field_offload_data* offload_data,
                  real* offload_array) {
}

void E_field_eval_E(real E[], real r, real phi, real z, E_field_data* Edata) {
  E[0] = 0.0;
  E[1] = 0.0;
  E[2] = 0.0;
}
