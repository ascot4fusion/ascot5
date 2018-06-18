/**
 * @file hdf5_neutral.c
 * @brief HDF5 format neutral density input
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../ascot5.h"
#include "../neutral.h"
#include "../neutral/N0_3D.h"
#include "../consts.h"
#include "../math.h"
#include "hdf5_neutral.h"
#include "hdf5.h"
#include "hdf5_helpers.h"
#include "hdf5_hl.h"

int hdf5_neutral_init_offload(hid_t f, neutral_offload_data* offload_data,
                              real** offload_array) {

    herr_t err;

#if VERBOSE > 0
    printf("\nReading neutral input from the HDF5 file...\n");
#endif
    
    err = hdf5_find_group(f, "/neutral/");
    if(err < 0) {
        return -1;
    }
    
    char active[11];
    err = H5LTget_attribute_string(f, "/neutral/", "active", active);
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
	
    hdf5_generate_qid_path("/neutral/N0_3D-XXXXXXXXXX", active, path);
    if(hdf5_find_group(f, path) == 0) { 
	hdf5_neutral_init_offload_3D(f, &(offload_data->N03D), offload_array, active);
	offload_data->type = neutral_type_3D;
	offload_data->offload_array_length = offload_data->N03D.offload_array_length;
	
#if VERBOSE > 0
        printf("\nLoaded 3D neutral density (N0_3D)\n");
	printf("with parameters:\n");
        printf("- rmin, rmax, nr = %le, %le, %d\n",
               offload_data->N03D.r_min,offload_data->N03D.r_max,offload_data->N03D.n_r);
        printf("- phimin, phimax, nphi = %le, %le, %d\n",
               offload_data->N03D.phi_min,offload_data->N03D.phi_max,offload_data->N03D.n_phi);
        printf("- zmin, zmax, nz = %le, %le, %d\n",
               offload_data->N03D.z_min,offload_data->N03D.z_max,offload_data->N03D.n_z);
#endif
        return 1;
    }

    hdf5_generate_qid_path("/neutral/N0_ST-XXXXXXXXXX", active, path);
    if(hdf5_find_group(f, path) == 0) { 
	hdf5_neutral_init_offload_ST(f, &(offload_data->N0ST), offload_array, active);
	offload_data->type = neutral_type_ST;
	offload_data->offload_array_length = offload_data->N0ST.offload_array_length;
	
#if VERBOSE > 0
        printf("\nLoaded stellarator neutral density (N0_ST)\n");
        printf("with parameters:\n");
        printf("- number of toroidal periods = %d\n",
               offload_data->N0ST.periods);
        printf("- rmin, rmax, nr = %le, %le, %d\n",
               offload_data->N0ST.r_min,offload_data->N0ST.r_max,offload_data->N0ST.n_r);
        printf("- phimin, phimax, nphi = %le, %le, %d\n",
               offload_data->N0ST.phi_min,offload_data->N0ST.phi_max,offload_data->N0ST.n_phi);
        printf("- zmin, zmax, nz = %le, %le, %d\n",
               offload_data->N0ST.z_min,offload_data->N0ST.z_max,offload_data->N0ST.n_z);
#endif
        return 1;
    }

    return -1;
}

/**
 * @brief Load neutral data from HDF5 file and prepare parameters
 *
 * This function reads the 3D neutral data from file f, fills the
 * offload struct with parameters and allocates and fills the offload array.
 *
 * @param f hdf5 file identifier
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 *
 */
void hdf5_neutral_init_offload_3D(hid_t f, N0_3D_offload_data* offload_data,
                                  real** offload_array, char* qid) {
    herr_t err;
    char path[256];

    /* Read and initialize magnetic field Rz-grid */
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/neutral/N0_3D-XXXXXXXXXX/n_r", qid, path), &(offload_data->n_r));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/neutral/N0_3D-XXXXXXXXXX/n_phi", qid, path), &(offload_data->n_phi));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/neutral/N0_3D-XXXXXXXXXX/n_z", qid, path), &(offload_data->n_z));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/neutral/N0_3D-XXXXXXXXXX/r_min", qid, path), &(offload_data->r_min));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/neutral/N0_3D-XXXXXXXXXX/r_max", qid, path), &(offload_data->r_max));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/neutral/N0_3D-XXXXXXXXXX/phi_min", qid, path), &(offload_data->phi_min));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/neutral/N0_3D-XXXXXXXXXX/phi_max", qid, path), &(offload_data->phi_max));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/neutral/N0_3D-XXXXXXXXXX/z_min", qid, path), &(offload_data->z_min));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/neutral/N0_3D-XXXXXXXXXX/z_max", qid, path), &(offload_data->z_max));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

    /* Convert degrees to radians */
    offload_data->phi_max = math_deg2rad(offload_data->phi_max);
    offload_data->phi_min = math_deg2rad(offload_data->phi_min);
    
    offload_data->r_grid = (offload_data->r_max - offload_data->r_min)
        / (offload_data->n_r - 1);
    offload_data->z_grid = (offload_data->z_max - offload_data->z_min)
        / (offload_data->n_z - 1);
    offload_data->phi_grid = (offload_data->phi_max - offload_data->phi_min)
        / (offload_data->n_phi - 1);

    *offload_array = (real*) malloc(offload_data->n_r * offload_data->n_phi * offload_data->n_z * sizeof(real));
    offload_data->offload_array_length = offload_data->n_r * offload_data->n_phi * offload_data->n_z;
    
    /* Read the neutral density */
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/neutral/N0_3D-XXXXXXXXXX/n0", qid, path), *offload_array);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

    return;
}

/**
 * @brief Load neutral data from HDF5 file and prepare parameters
 *
 * This function reads the stellarator neutral data from file f, fills the
 * offload struct with parameters and allocates and fills the offload array.
 *
 * @param f hdf5 file identifier
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 *
 */
void hdf5_neutral_init_offload_ST(hid_t f, N0_ST_offload_data* offload_data,
                                  real** offload_array, char* qid) {
    herr_t err;
    char path[256];
    
    /* Number of toroidal periods */
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/neutral/N0_ST-XXXXXXXXXX/toroidalPeriods", qid, path), &(offload_data->periods));

    /* Read and initialize magnetic field Rz-grid */
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/neutral/N0_ST-XXXXXXXXXX/n_r", qid, path), &(offload_data->n_r));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/neutral/N0_ST-XXXXXXXXXX/n_phi", qid, path), &(offload_data->n_phi));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/neutral/N0_ST-XXXXXXXXXX/n_z", qid, path), &(offload_data->n_z));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/neutral/N0_ST-XXXXXXXXXX/r_min", qid, path), &(offload_data->r_min));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/neutral/N0_ST-XXXXXXXXXX/r_max", qid, path), &(offload_data->r_max));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/neutral/N0_ST-XXXXXXXXXX/phi_min", qid, path), &(offload_data->phi_min));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/neutral/N0_ST-XXXXXXXXXX/phi_max", qid, path), &(offload_data->phi_max));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/neutral/N0_ST-XXXXXXXXXX/z_min", qid, path), &(offload_data->z_min));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/neutral/N0_ST-XXXXXXXXXX/z_max", qid, path), &(offload_data->z_max));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

    /* Convert degrees to radians */
    offload_data->phi_max = math_deg2rad(offload_data->phi_max);
    offload_data->phi_min = math_deg2rad(offload_data->phi_min);
    
    offload_data->r_grid = (offload_data->r_max - offload_data->r_min)
        / (offload_data->n_r - 1);
    offload_data->z_grid = (offload_data->z_max - offload_data->z_min)
        / (offload_data->n_z - 1);
    offload_data->phi_grid = (offload_data->phi_max - offload_data->phi_min)
        / (offload_data->n_phi - 1);

    *offload_array = (real*) malloc(offload_data->n_r * offload_data->n_phi * offload_data->n_z * sizeof(real));
    offload_data->offload_array_length = offload_data->n_r * offload_data->n_phi * offload_data->n_z;
    
    /* Read the neutral density */
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/neutral/N0_ST-XXXXXXXXXX/n0", qid, path), *offload_array);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

    int temp_n0_size = offload_data->n_r*offload_data->n_z*offload_data->n_phi;
        
    real* temp_n0   = (real*) malloc(temp_n0_size*sizeof(real));
        
    /* Read the neutral density */
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/neutral/N0_ST-XXXXXXXXXX/n0", qid, path), temp_n0);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    /* We need to use stellarator symmetry here.
     * http://dx.doi.org/10.1016/S0167-2789(97)00216-9
     * The data is expected to include half a period.
     */
    int n0_size = offload_data->n_r*offload_data->n_z*(2*(offload_data->n_phi - 1));
    *offload_array = (real*) malloc(n0_size * sizeof(real));
    offload_data->offload_array_length = n0_size;
    
    int i_phi;
    int i_z;
    int i_r;
    int temp_ind, off_ind, sym_ind;
    for (i_phi = 0; i_phi < offload_data->n_phi; i_phi++) {
        for (i_z = 0; i_z < offload_data->n_z; i_z++) {
            for (i_r = 0; i_r < offload_data->n_r; i_r++) {
                /* Stellarator symmetry: i_phi        <=> 2*(n_phi-1)-i_phi
                 *                       i_z          <=> n_z-i_z-1
                 * So, a point at:      (i_r, i_phi, i_z)  <=> (i_r, 2*(n_phi-1)-i_phi, n_z-i_z-1)
                 *
                 * temp_n0 data is in the format: (i_r, i_phi, i_z) = temp_n0(i_z*n_phi*n_r + i_phi*n_r + i_r)
                 * offload_array data -"-"-"-"-  : (i_r, i_phi, i_z) = (*offload_array)[i_phi*n_z*n_r + i_z*n_r + i_r ]
                 * => (*offload_array)[i_phi*n_z*n_r + i_z*n_r + i_r ] = temp_n0(i_z*n_phi*n_r + i_phi*n_r + i_r);
                 * The values are: Sym[B_r, B_phi, B_z] = [-B_r, B_phi, B_z]
                 */

                /* Index of data point in temp arrays */
                temp_ind = i_z*offload_data->n_phi*offload_data->n_r + i_phi*offload_data->n_r + i_r;

                /* Index of data point in offload_array and corresponding stel.-symmetric index  */
                off_ind = i_phi*offload_data->n_z*offload_data->n_r + i_z*offload_data->n_r + i_r;
                sym_ind = (2*(offload_data->n_phi - 1) - i_phi)*offload_data->n_z*offload_data->n_r
                    + (offload_data->n_z - i_z - 1)*offload_data->n_r + i_r;

                (*offload_array)[off_ind] =  temp_n0[temp_ind];
                if (i_phi != 0) {
                    (*offload_array)[sym_ind] = -temp_n0[temp_ind];
                }
            }
        }
    }
    /* Phi data is now for one toroidal period */
    offload_data->n_phi = 2*(offload_data->n_phi - 1);
    offload_data->phi_max = 2*offload_data->phi_max;
    free(temp_n0);
    
}

