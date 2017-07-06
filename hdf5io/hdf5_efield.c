/**
 * @file hdf5_efield.c
 * @brief HDF5 format 1d radial electric field input
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include "hdf5_hl.h"
#include "../E_field.h"
#include "hdf5_efield.h"


int hdf5_efield_init_offload(hid_t f, E_field_offload_data* offload_data, real** offload_array) {
    herr_t err;

    err = H5LTfind_dataset(f, "/efield/");
    if(err < 0) {
        return -1;
    }

    char type[32];
    err = H5LTget_attribute_string(f, "/efield/", "type", type);
    if(err < 0) {
        return -1;
    }
    if(strncmp(type,"E_TC",4) == 0) {
        offload_data->type = E_field_type_TC;
	hdf5_efield_init_offload_TC(f, &(offload_data->ETC), offload_array);

	#if VERBOSE > 0
	    printf("\nLoaded trivial cartesian electric field (E_TC)\n");
	    printf("with parameters:\n");
	    printf("- Electric field Exyz = (%le, %le, %le)\n",
		   (*offload_array)[0],(*offload_array)[1],(*offload_array)[2]);
	#endif

        return 1;
    }
    if(strcmp(type,"erad") == 0) {
        offload_data->type = E_field_type_1D;
	hdf5_efield_init_offload_1D(f, &(offload_data->E1D), offload_array);
        return 1;
    }

    return -1;
}

void hdf5_efield_init_offload_1D(hid_t f, E_1D_offload_data* offload_data, real** offload_array) {
    
    herr_t err;

    err = H5LTget_attribute_int(f, "/efield/erad/", "n_rho", &(offload_data->n_rho));
    err = H5LTget_attribute_double(f, "/efield/erad/", "rho_min", &(offload_data->rho_min));
    err = H5LTget_attribute_double(f, "/efield/erad/", "rho_max", &(offload_data->rho_max));
    
    /* Allocate n_rho space for both rho and dV/drho */
    offload_data->offload_array_length = 2*offload_data->n_rho;
    *offload_array = (real*) malloc(2*offload_data->n_rho*sizeof(real));

    /* Pointers to beginning of different data series to make code more
     * readable */
    real* rho = &(*offload_array)[0];
    real* dV = &(*offload_array)[offload_data->n_rho];
    
    err = H5LTread_dataset_double(f,"/efield/erad/rho",rho);
    err = H5LTread_dataset_double(f,"/efield/erad/dV_drho",dV);

    /* Effective minor radius */
    real r_eff;
    err = H5LTget_attribute_double(f, "/efield/erad/", "r_eff", &(r_eff));    
    /* Scale derivatives by effective minor radius */
    for(int i = 0; i < offload_data->n_rho; i++) {
        dV[i] = r_eff * dV[i];
    }
}

void hdf5_efield_init_offload_TC(hid_t f, E_TC_offload_data* offload_data, real** offload_array) {
    
    herr_t err;

    *offload_array = (real*) malloc(3*sizeof(real));
    err = H5LTread_dataset_double(f,"/efield/etc/Exyz", *offload_array);

    offload_data->offload_array_length = 3;
}
