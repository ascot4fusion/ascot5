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
#include "../Efield/E_TC.h"
#include "../Efield/E_1D.h"
#include "../Efield/E_1DS.h"
#include "hdf5_efield.h"
#include "hdf5_helpers.h"


int hdf5_efield_init_offload(hid_t f, E_field_offload_data* offload_data, real** offload_array) {
    herr_t err;

    #if VERBOSE > 0
        printf("\nReading electric field input from the HDF5 file...\n");
    #endif
    
    err = hdf5_find_group(f, "/efield/");
    if(err < 0) {
        return -1;
    }
    
    char active[11];
    err = H5LTget_attribute_string(f, "/efield/", "active", active);
    if(err < 0) {
        return -1;
    }
    active[10] = '\0';

    #if VERBOSE > 0
        printf("Active qid is %s\n", active);
    #endif

    /* Go through all different input types and see which one the active qid corresponds to.
     * Then read this input. */
    char path[256];

    hdf5_generate_qid_path("/efield/E_TC-XXXXXXXXXX", active, path);
    if(hdf5_find_group(f, path) == 0) {	
        offload_data->type = E_field_type_TC;
	hdf5_efield_init_offload_TC(f, &(offload_data->ETC), offload_array, active);
	offload_data->offload_array_length = offload_data->ETC.offload_array_length;

	#if VERBOSE > 0
	    printf("\nLoaded trivial cartesian electric field (E_TC)\n");
	    printf("with parameters:\n");
	    printf("- Electric field Exyz = (%le, %le, %le)\n",
		   (*offload_array)[0],(*offload_array)[1],(*offload_array)[2]);
	#endif

        return 1;
    }

    hdf5_generate_qid_path("/efield/E_1DS-XXXXXXXXXX", active, path);
    if(hdf5_find_group(f, path) == 0) {
        offload_data->type = E_field_type_1DS;
	hdf5_efield_init_offload_1DS(f, &(offload_data->E1DS), offload_array, active);
	offload_data->offload_array_length = offload_data->E1DS.offload_array_length;

	#if VERBOSE > 0
	    printf("\nLoaded radial electric field (E_1DS)\n");
	    printf("with parameters:\n");
	    printf("(n_rho, rho_min, rho_max) = (%d, %le, %le)\n",
		   offload_data->E1DS.n_rho,offload_data->E1DS.rho_min,offload_data->E1DS.rho_max);
	#endif
            
        return 1;
    }

    return -1;
}

void hdf5_efield_init_offload_1D(hid_t f, E_1D_offload_data* offload_data, real** offload_array, char* qid) {
    
    herr_t err;
    char path[256];

    err = H5LTget_attribute_int(f, hdf5_generate_qid_path("/efield/E_1D-XXXXXXXXXX/", qid, path), "n_rho", &(offload_data->n_rho));
    err = H5LTget_attribute_double(f, hdf5_generate_qid_path("/efield/E_1D-XXXXXXXXXX/", qid, path), "rho_min", &(offload_data->rho_min));
    err = H5LTget_attribute_double(f, hdf5_generate_qid_path("/efield/E_1D-XXXXXXXXXX/", qid, path), "rho_max", &(offload_data->rho_max));
    
    /* Allocate n_rho space for both rho and dV/drho */
    offload_data->offload_array_length = 2*offload_data->n_rho;
    *offload_array = (real*) malloc(2*offload_data->n_rho*sizeof(real));

    /* Pointers to beginning of different data series to make code more
     * readable */
    real* rho = &(*offload_array)[0];
    real* dV = &(*offload_array)[offload_data->n_rho];
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/efield/E_1D-XXXXXXXXXX/rho", qid, path), rho);
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/efield/E_1D-XXXXXXXXXX/dV_drho", qid, path), dV);

    /* Effective minor radius */
    real r_eff;
    err = H5LTget_attribute_double(f, hdf5_generate_qid_path("/efield/E_1D-XXXXXXXXXX/", qid, path), "r_eff", &(r_eff));
    
    /* Scale derivatives by effective minor radius */
    for(int i = 0; i < offload_data->n_rho; i++) {
        dV[i] = r_eff * dV[i];
    }
}

void hdf5_efield_init_offload_1DS(hid_t f, E_1DS_offload_data* offload_data, real** offload_array, char* qid) {    
    herr_t err;
    char path[256];

    err = H5LTget_attribute_int(f, hdf5_generate_qid_path("/efield/E_1DS-XXXXXXXXXX/", qid, path), "n_rho", &(offload_data->n_rho));
    err = H5LTget_attribute_double(f, hdf5_generate_qid_path("/efield/E_1DS-XXXXXXXXXX/", qid, path), "rho_min", &(offload_data->rho_min));
    err = H5LTget_attribute_double(f, hdf5_generate_qid_path("/efield/E_1DS-XXXXXXXXXX/", qid, path), "rho_max", &(offload_data->rho_max));
    
    /* Allocate n_rho space for both rho and dV/drho */
    offload_data->offload_array_length = 2*offload_data->n_rho;
    *offload_array = (real*) malloc(2*offload_data->n_rho*sizeof(real));

    /* Pointers to beginning of different data series to make code more
     * readable */
    real* rho = &(*offload_array)[0];
    real* dV = &(*offload_array)[offload_data->n_rho];
    
    err = H5LTread_dataset_double(f,hdf5_generate_qid_path("/efield/E_1DS-XXXXXXXXXX/rho", qid, path),rho);
    err = H5LTread_dataset_double(f,hdf5_generate_qid_path("/efield/E_1DS-XXXXXXXXXX/dV_drho", qid, path),dV);

    /* Effective minor radius */
    real r_eff;
    err = H5LTget_attribute_double(f, hdf5_generate_qid_path("/efield/E_1DS-XXXXXXXXXX/", qid, path), "r_eff", &(r_eff));
    
    /* Scale derivatives by effective minor radius */
    for(int i = 0; i < offload_data->n_rho; i++) {
        dV[i] = r_eff * dV[i];
    }
}

void hdf5_efield_init_offload_TC(hid_t f, E_TC_offload_data* offload_data, real** offload_array, char* qid) {
    
    herr_t err;
    char path[256];
    
    *offload_array = (real*) malloc(3*sizeof(real));
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/efield/E_TC-XXXXXXXXXX/Exyz", qid, path), &(*offload_array)[0]);
    offload_data->offload_array_length = 3;
}
