/**
 * @file hdf5_efield.c
 * @brief HDF5 format 1d/3d radial electric field input
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include "hdf5_hl.h" //should I include ascot.5h and math.h?
#include "../math.h"
#include "../E_field.h"
#include "../Efield/E_TC.h"
#include "../Efield/E_1D.h"
#include "../Efield/E_1DS.h"
#include "hdf5_efield.h"
#include "hdf5_helpers.h"

/**
 * @brief Initialize electric field offload data from h5 file
 *
 * This function reads the electric field data from the input.h5 file,
 * fills the offload struct with parameters and
 * allocates and fills the offload array.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */

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

    hdf5_generate_qid_path("/efield/E_3D-XXXXXXXXXX", active, path);
    if(hdf5_find_group(f, path) == 0) {
      hdf5_efield_init_offload_3D(f, &(offload_data->E3D), offload_array, active);
      offload_data->type = E_field_type_3D;
      offload_data->offload_array_length=offload_data->E3D.offload_array_length;

        #if VERBOSE > 0
      printf("\nLoaded 3D electric field (E_3D)\n");
      printf("with parameters:\n");
      printf("- rmin, rmax, nr = %le, %le, %d\n",
	     offload_data->E3D.r_min,offload_data->E3D.r_max,offload_data->E3D.n_r);
      printf("- zmin, zmax, nz = %le, %le, %d\n",
	     offload_data->E3D.z_min,offload_data->E3D.z_max,offload_data->E3D.n_z);
      printf("- phimin, phimax, nphi = %le, %le, %d\n",
             offload_data->E3D.phi_min,offload_data->E3D.phi_max,offload_data->E3D.n_phi);

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


void hdf5_efield_init_offload_3D(hid_t f, E_3D_offload_data* offload_data, real** offload_array, char* qid) {
herr_t err;
char path[256];

/* Read and initialize magnetic field Rz-grid */
err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/efield/E_3D-XXXXXXXXXX/n_R", qid, path), &(offload_data->n_r));
if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/efield/E_3D-XXXXXXXXXX/n_phi", qid, path), &(offload_data->n_phi));
if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/efield/E_3D-XXXXXXXXXX/n_z", qid, path), &(offload_data->n_z));
if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/efield/E_3D-XXXXXXXXXX/R_min", qid, path), &(offload_data->r_min));
if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/efield/E_3D-XXXXXXXXXX/R_max", qid, path), &(offload_data->r_max));
if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/efield/E_3D-XXXXXXXXXX/phi_min", qid, path), &(offload_data->phi_min));
if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/efield/E_3D-XXXXXXXXXX/phi_max", qid, path), &(offload_data->phi_max));
if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/efield/E_3D-XXXXXXXXXX/z_min", qid, path), &(offload_data->z_min));
if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/efield/E_3D-XXXXXXXXXX/z_max", qid, path), &(offload_data->z_max));
if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}


/* Convert degrees to radians (input must be in deg but ASCOT internally works in rads) */


offload_data->phi_max = math_deg2rad(offload_data->phi_max);
offload_data->phi_min = math_deg2rad(offload_data->phi_min);


/* Calculate grid size */

offload_data->r_grid = (offload_data->r_max - offload_data->r_min) / (offload_data->n_r - 1);
offload_data->z_grid = (offload_data->z_max - offload_data->z_min) / (offload_data->n_z - 1);
offload_data->phi_grid = (offload_data->phi_max - offload_data->phi_min) / (offload_data->n_phi - 1);


/* Allocate offload_array */
int E_size = offload_data->n_r*offload_data->n_z*offload_data->n_phi;

*offload_array = (real*) malloc((3 * E_size) * sizeof(real));
offload_data->offload_array_length = 3 * E_size;

/* Read the magnetic field */
 err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/efield/E_3D-XXXXXXXXXX/E_R", qid, path), &(*offload_array)[0*E_size]);
 if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

 err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/efield/E_3D-XXXXXXXXXX/E_phi", qid, path), &(*offload_array)[1*E_size]);
 if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

 err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/efield/E_3D-XXXXXXXXXX/E_z", qid, path), &(*offload_array)[2*E_size]);
 if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

}


