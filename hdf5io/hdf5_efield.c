/**
 * @file hdf5_efield.c
<<<<<<< HEAD
 * @brief HDF5 format 1d/3d radial electric field input
=======
 * @brief Module for reading electric field data from HDF5 file
>>>>>>> develop
 *
 * Electric field  must be read by calling hdf5_efield_init_offload() contained
 * in this module. This module contains reading routines for all electric field
 * types.
 *
 * Note: the routines within this module that read the data from HDF5 file
 * may use the offload data struct and offload arrays as (temporary) storage.
 * However, the actual initialization is done at the specific
 * E_field_init_offload() function the electric field data corresponds to. Check
 * from that function what the offload data and the offload array are expected
 * to contain. As a rule of thumb, the reading routines here should only read
 * the data and maybe do some trivial computations but nothing complicated.
 */
#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>
<<<<<<< HEAD
<<<<<<< HEAD
#include "hdf5_hl.h" //should I include ascot.5h and math.h?
#include "../math.h"
=======
#include "hdf5_hl.h"
#include "../ascot5.h"
#include "../print.h"
>>>>>>> develop
=======
#include <hdf5_hl.h>
#include "../ascot5.h"
>>>>>>> develop
#include "../E_field.h"
#include "../Efield/E_TC.h"
#include "../Efield/E_1D.h"
#include "../Efield/E_1DS.h"
#include "../print.h"
#include "hdf5_helpers.h"
#include "hdf5_efield.h"

<<<<<<< HEAD
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
=======
int hdf5_efield_read_1D(hid_t f, E_1D_offload_data* offload_data,
                        real** offload_array, char* qid);
int hdf5_efield_read_1DS(hid_t f, E_1DS_offload_data* offload_data,
                         real** offload_array, char* qid);
int hdf5_efield_read_TC(hid_t f, E_TC_offload_data* offload_data,
                        real** offload_array, char* qid);
/**
 * @brief Read electric field data from HDF5 file
 *
 * This function reads electric field data with given qid while also
 * initializing offload data and allocating and filling offload array. The file
 * is opened and closed outside this function.
 *
 * The electric field data is stored under /efield/ group in ASCOT5 HDF5 file.
 * Several electric fields of same or different type maybe stored there in
 * their respective groups as long as the group name contains QID as an
 * identifier.
 *
 * @param f HDF5 file identifier for a file which is opened and closed outside
 *          of this function
 * @param offload_data pointer to offload data struct which is initialized here
 * @param offload_array pointer to offload array which is allocated and
 *                      initialized here
 * @param qid QID of the data that is to be read
 *
 * @return zero if reading and initialization succeeded
 */
int hdf5_efield_init_offload(hid_t f, E_field_offload_data* offload_data,
                             real** offload_array, char* qid) {
>>>>>>> develop
    char path[256];
    int err = 1;

    /* Read data the QID corresponds to */

    hdf5_gen_path("/efield/E_TC-XXXXXXXXXX", qid, path);
    if( !hdf5_find_group(f, path) ) {
        offload_data->type = E_field_type_TC;
        err = hdf5_efield_read_TC(f, &(offload_data->ETC),
                                  offload_array, qid);
    }

    hdf5_gen_path("/efield/E_1DS-XXXXXXXXXX", qid, path);
    if(hdf5_find_group(f, path) == 0) {
        offload_data->type = E_field_type_1DS;
        err = hdf5_efield_read_1DS(f, &(offload_data->E1DS),
                                   offload_array, qid);
    }

    /* Initialize if data was read succesfully */
    if(!err) {
        err = E_field_init_offload(offload_data, offload_array);
    }

<<<<<<< HEAD
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
=======
    return err;
>>>>>>> develop
}

/**
 * @brief Read E1D electric field data from HDF5 file.
 *
 * @param f HDF5 file from which data is read
 * @param offload_data pointer to offload data
 * @param offload_array pointer to offload array
 * @param qid QID of the data
 *
 * @return Zero if reading and initialization of data succeeded
 */
int hdf5_efield_read_1D(hid_t f, E_1D_offload_data* offload_data,
                        real** offload_array, char* qid) {

    herr_t err;
    char path[256];

    err = H5LTget_attribute_int(f, hdf5_generate_qid_path("/efield/E_1D-XXXXXXXXXX/", qid, path), "n_rho", &(offload_data->n_rho));
    if(err) {print_err("Error: Failed to read dataset."); return 1;}
    err = H5LTget_attribute_double(f, hdf5_generate_qid_path("/efield/E_1D-XXXXXXXXXX/", qid, path), "rho_min", &(offload_data->rho_min));
    if(err) {print_err("Error: Failed to read dataset."); return 1;}
    err = H5LTget_attribute_double(f, hdf5_generate_qid_path("/efield/E_1D-XXXXXXXXXX/", qid, path), "rho_max", &(offload_data->rho_max));
    if(err) {print_err("Error: Failed to read dataset."); return 1;}

    /* Allocate n_rho space for both rho and dV/drho */
    offload_data->offload_array_length = 2*offload_data->n_rho;
    *offload_array = (real*) malloc(2*offload_data->n_rho*sizeof(real));

    /* Pointers to beginning of different data series to make code more
     * readable */
    real* rho = &(*offload_array)[0];
    real* dV = &(*offload_array)[offload_data->n_rho];

    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/efield/E_1D-XXXXXXXXXX/rho", qid, path), rho);
    if(err) {print_err("Error: Failed to read dataset."); return 1;}
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/efield/E_1D-XXXXXXXXXX/dV_drho", qid, path), dV);
    if(err) {print_err("Error: Failed to read dataset."); return 1;}

    /* Effective minor radius */
    real r_eff;
    err = H5LTget_attribute_double(f, hdf5_generate_qid_path("/efield/E_1D-XXXXXXXXXX/", qid, path), "r_eff", &(r_eff));
    if(err) {print_err("Error: Failed to read dataset."); return 1;}

    /* Scale derivatives by effective minor radius */
    for(int i = 0; i < offload_data->n_rho; i++) {
        dV[i] = r_eff * dV[i];
    }

    return err;
}

/**
 * @brief Read E1DS electric field data from HDF5 file.
 *
 * @param f HDF5 file from which data is read
 * @param offload_data pointer to offload data
 * @param offload_array pointer to offload array
 * @param qid QID of the data
 *
 * @return Zero if reading succeeded
 */
int hdf5_efield_read_1DS(hid_t f, E_1DS_offload_data* offload_data,
                         real** offload_array, char* qid) {
    #undef EPATH
    #define EPATH "/efield/E_1DS-XXXXXXXXXX/"

    if( hdf5_read_int(EPATH "n_rho", &(offload_data->n_rho),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(EPATH "rho_min", &(offload_data->rho_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(EPATH "rho_max", &(offload_data->rho_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Allocate n_rho space for dV/drho */
    *offload_array = (real*) malloc(offload_data->n_rho*sizeof(real));

    if( hdf5_read_double(EPATH "dV_drho", *offload_array,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Effective minor radius */
    real r_eff;
    if( hdf5_read_double(EPATH "r_eff", &(r_eff),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Scale derivatives by effective minor radius */
    for(int i = 0; i < offload_data->n_rho; i++) {
        (*offload_array)[i] = r_eff * (*offload_array)[i];
    }

    return 0;
}

<<<<<<< HEAD
void hdf5_efield_init_offload_TC(hid_t f, E_TC_offload_data* offload_data, real** offload_array, char* qid) {
    
    herr_t err;
    char path[256];
    
    *offload_array = (real*) malloc(3*sizeof(real));
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/efield/E_TC-XXXXXXXXXX/Exyz", qid, path), &(*offload_array)[0]);
<<<<<<< HEAD
=======
    if(err) {print_err("Error: Failed to read dataset."); return;}
    
>>>>>>> develop
    offload_data->offload_array_length = 3;
=======
/**
 * @brief Read magnetic field data of type E_TC
 *
 * The E_TC data is stored in HDF5 file under the group
 * /efield/E_TC-XXXXXXXXXX/ where X's mark the QID.
 *
 * This function assumes the group holds the following datasets:
 *
 * - double Exyz Electric field values [E_x, E_y, E_z] [V/m]
 *
 * @param f HDF5 file identifier for a file which is opened and closed outside
 *          of this function
 * @param offload_data pointer to offload data struct which is allocated here
 * @param offload_array pointer to offload array but no data is stored there so
 *                      it is not allocated and NULL pointer is returned instead
 * @param qid QID of the B_TC field that is to be read
 *
 * @return zero if reading succeeded
 */
int hdf5_efield_read_TC(hid_t f, E_TC_offload_data* offload_data,
                        real** offload_array, char* qid) {
    #undef EPATH
    #define EPATH "/efield/E_TC-XXXXXXXXXX/"

    *offload_array = NULL;

    if( hdf5_read_double(EPATH "Exyz", offload_data->Exyz,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    return 0;
>>>>>>> develop
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


