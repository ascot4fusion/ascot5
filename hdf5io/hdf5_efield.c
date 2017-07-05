/**
 * @file hdf5_efield.c
 * @brief HDF5 format 1d radial electric field input
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>
#include "hdf5_hl.h"
#include "../E_field.h"


void hdf5_efield_init_offload(hid_t f, E_field_offload_data* offload_data, real** offload_array) {
    herr_t err;

    err = H5LTfind_dataset(f, "/erad/");
    hdf5_efield_init_offload_1D(f, &(offload_data->E1D), offload_array);


}

void hdf5_efield_init_offload_1D(hid_t f, E_1D_offload_data* offload_data, real** offload_array) {
    
    herr_t err;

    err = H5LTget_attribute_int(f, "/erad/", "n_rho", &(offload_data->n_rho));
    err = H5LTget_attribute_double(f, "/erad/", "rho_min", &(offload_data->rho_min));
    err = H5LTget_attribute_double(f, "/erad/", "rho_max", &(offload_data->rho_max));
    
    /* Allocate n_rho space for both rho and dV/drho */
    offload_data->offload_array_length = 2*offload_data->n_rho;
    *offload_array = (real*) malloc(2*offload_data->n_rho*sizeof(real));

    /* Pointers to beginning of different data series to make code more
     * readable */
    real* rho = &(*offload_array)[0];
    real* dV = &(*offload_array)[offload_data->n_rho];
    
    err = H5LTread_dataset_double(f,"/bfield/rho",rho);
    err = H5LTread_dataset_double(f,"/bfield/dV_drho",dV);

    /* Effective minor radius */
    real r_eff;
    err = H5LTget_attribute_double(f, "/erad/", "r_eff", &(r_eff));    
    /* Scale derivatives by effective minor radius */
    for(int i = 0; i < offload_data->n_rho; i++) {
        dV[i] = r_eff * dV[i];
    }
}
