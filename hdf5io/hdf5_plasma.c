/**
 * @file hdf5_plasma_1d.c
 * @brief HDF5 format 1d plasma input
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../ascot5.h"
#include "../plasma.h"
#include "../plasma/plasma_1D.h"
#include "../plasma/plasma_1DS.h"
#include "../consts.h"
#include "hdf5.h"
#include "hdf5_helpers.h"
#include "hdf5_hl.h"


int hdf5_plasma_init_offload(hid_t f, plasma_offload_data* offload_data,
			     real** offload_array) {

    herr_t err;
    hsize_t dims[3];

    #if VERBOSE > 0
        printf("Reading plasma input from the HDF5 file...\n");
    #endif
    
    err = hdf5_find_group(f, "/plasma/");
    if(err < 0) {
        return -1;
    }
    
    char active[11];
    err = H5LTget_attribute_string(f, "/plasma/", "active", active);
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
	
    hdf5_generate_qid_path("/plasma/plasma_1D-XXXXXXXXXX", active, path);
    if(hdf5_find_group(f, path) == 0) { 
	hdf5_plasma_init_offload_1D(f, &(offload_data->plasma_1D), offload_array, active);
	offload_data->type = plasma_type_1D;
	offload_data->offload_array_length = offload_data->plasma_1D.offload_array_length;
        return 1;
    }

    return -1;
}

/**
 * @brief Load plasma data from HDF5 file and prepare parameters
 *
 * This function reads the 1D plasma data from file f, fills the
 * offload struct with parameters and allocates and fills the offload array.
 *
 * @param f hdf5 file identifier
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 *
 */
int hdf5_plasma_init_offload_1D(hid_t f, plasma_1D_offload_data* offload_data,
				real** offload_array, char* qid) {
    herr_t err;
    char path[256];
    hsize_t dims[3];

    int i, j;

    int n_ions;
    err = H5LTget_attribute_int(f, hdf5_generate_qid_path("/plasma/plasma_1D-XXXXXXXXXX", qid, path), "n_ions", &n_ions);
    offload_data->n_species = n_ions + 1; /* Include electrons */
    
    err = H5LTget_attribute_int(f, hdf5_generate_qid_path("/plasma/plasma_1D-XXXXXXXXXX", qid, path), "n_rho", &(offload_data->n_rho));
    int n_rho = offload_data->n_rho;
    
    /* Allocate space for rho + temperature for each ion species and electrons
     * + density for each ion species and electrons */
    offload_data->offload_array_length = n_rho + 2*offload_data->n_species*n_rho;
    *offload_array = (real*) malloc(sizeof(real)
				    * offload_data->offload_array_length);
    
    /* Pointers to beginning of different data series to make code more
     * readable */
    real* rho = &(*offload_array)[0];
    real* temp_e = &(*offload_array)[n_rho];
    real* temp_i = &(*offload_array)[n_rho*2];
    real* dens_e = &(*offload_array)[n_rho*2
				     + n_rho*n_ions];
    real* dens_i = &(*offload_array)[n_rho*2
				     + n_rho*n_ions
				     + n_rho];

    /* Read Znum and calculate charge */
    int Znum[n_ions];
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/plasma/plasma_1D-XXXXXXXXXX/Z_num", qid, path), Znum);
    /* Electron charge -1 */
    offload_data->charge[0] = -1 * CONST_E;
    for(i = 0; i < n_ions; i++) {
	offload_data->charge[i+1] = Znum[i] * CONST_E;
    }

    /* Read Amass and calculate mass */
    int Amass[n_ions];
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/plasma/plasma_1D-XXXXXXXXXX/A_mass", qid, path), Amass);
    offload_data->mass[0] = CONST_M_E;
    for(i = 0; i < n_ions; i++) {
	offload_data->mass[i+1] = Amass[i] * CONST_U;
    }
    
    /* Read actual data into array */
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/plasma/plasma_1D-XXXXXXXXXX/rho", qid, path), rho);
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/plasma/plasma_1D-XXXXXXXXXX/temp_e", qid, path), temp_e);
    for(i = 0; i < n_rho; i++) {
	temp_e[i] = temp_e[i] * CONST_E / CONST_KB;
    }
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/plasma/plasma_1D-XXXXXXXXXX/dens_e", qid, path), dens_e);

    /* All ions have same temperature */
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/plasma/plasma_1D-XXXXXXXXXX/temp_i", qid, path), temp_i);
    for(i = 0; i < n_rho; i++) {
	temp_i[i] = temp_i[i] * CONST_E / CONST_KB;
    }
    
    real temp_dens_i[n_ions*n_rho];
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/plasma/plasma_1D-XXXXXXXXXX/dens_i", qid, path), temp_dens_i);
    for(j = 0; j < n_ions; j++) {
	for(i = 0; i < n_rho; i++) {
	    dens_i[j*n_rho + i] = temp_dens_i[j*n_rho + i];
	}
    }

    #if VERBOSE > 0
    printf("\nLoaded 1D plasma profiles (P_1D)\n");
    printf("with parameters:\n");
    printf("- %d number of rho values ranging from %le to %le\n",
	   n_rho, rho[0], rho[n_rho-1]);
    printf("- Number of ion species %d:\n",
	   n_ions); 
    for(int k=0; k<n_ions; k++) {
	printf("  - Znum %d, Amass %d\n",Znum[k],Amass[k]);
    }
    //printf("- Number of neutral species %d, Znum %d, Anum %d\n");
    printf("- Central electron temperature %le eV and density %le m^-3\n",
	   temp_e[0] * CONST_KB / CONST_E, dens_e[0]);
    printf("- Central ion temperature %le eV and densities [m^-3]\n", 
	   temp_i[0] * CONST_KB / CONST_E);
    for(int k=0; k<n_ions; k++) {
	printf("  - %le\n",dens_i[k*n_rho]);
    }
    //printf("- Central neutral temperature %le and density %le\n");
    #endif

    return 0;
}

/**
 * @brief Load plasma data from HDF5 file and prepare parameters
 *
 * This function reads the 1D plasma data from file f, fills the
 * offload struct with parameters and allocates and fills the offload array.
 *
 * @param f hdf5 file identifier
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 *
 */
int hdf5_plasma_init_offload_1DS(hid_t f, plasma_1DS_offload_data* offload_data,
                                 real** offload_array, char* qid) {
    char path[256];
    a5err err = 0;
    int i, j;
    err = H5LTfind_dataset(f, "/plasma/P_1D");
    
    int n_ions;
    err = H5LTget_attribute_int(f, "/plasma/", "n_ions", &n_ions);
    offload_data->n_species = n_ions + 1; /* Include electrons */
    
    err = H5LTget_attribute_int(f, "/plasma/P_1D", "n_rho", &(offload_data->n_rho));
    int n_rho = offload_data->n_rho;
    
    /* Allocate space for rho + temperature for each ion species and electrons
     * + density for each ion species and electrons */
    offload_data->offload_array_length = n_rho + 2*offload_data->n_species*n_rho;
    *offload_array = (real*) malloc(sizeof(real)
                                    * offload_data->offload_array_length);
    
    /* Pointers to beginning of different data series to make code more
     * readable */
    real* rho = &(*offload_array)[0];
    real* temp_e = &(*offload_array)[n_rho];
    real* temp_i = &(*offload_array)[n_rho*2];
    real* dens_e = &(*offload_array)[n_rho*2
                                     + n_rho*n_ions];
    real* dens_i = &(*offload_array)[n_rho*2
                                     + n_rho*n_ions
                                     + n_rho];

    /* Read Znum and calculate charge */
    int Znum[n_ions];
    err = H5LTread_dataset_int(f,"/plasma/Z_num",Znum);
    /* Electron charge -1 */
    offload_data->charge[0] = -1 * CONST_E;
    for(i = 0; i < n_ions; i++) {
        offload_data->charge[i+1] = Znum[i] * CONST_E;
    }

    /* Read Amass and calculate mass */
    int Amass[n_ions];
    err = H5LTread_dataset_int(f,"/plasma/A_mass", Amass);
    offload_data->mass[0] = CONST_M_E;
    for(i = 0; i < n_ions; i++) {
        offload_data->mass[i+1] = Amass[i] * CONST_U;
    }
    
    /* Read actual data into array */
    err = H5LTread_dataset_double(f,"/plasma/P_1D/rho", rho);
    err = H5LTread_dataset_double(f,"/plasma/P_1D/temp_e", temp_e);
    for(i = 0; i < n_rho; i++) {
        temp_e[i] = temp_e[i] * CONST_E / CONST_KB;
    }
    err = H5LTread_dataset_double(f,"/plasma/P_1D/dens_e", dens_e);

    /* All ions have same temperature */
    err = H5LTread_dataset_double(f,"/plasma/P_1D/temp_i", temp_i);
    for(i = 0; i < n_rho; i++) {
        temp_i[i] = temp_i[i] * CONST_E / CONST_KB;
    }
    
    real temp_dens_i[n_ions*n_rho];
    err = H5LTread_dataset_double(f,"/plasma/P_1D/dens_i", temp_dens_i);
    for(j = 0; j < n_ions; j++) {
        for(i = 0; i < n_rho; i++) {
            dens_i[j*n_rho + i] = temp_dens_i[j*n_rho + i];
        }
    }
    
    return 0;
}
