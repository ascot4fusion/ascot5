/**
 * @file hdf5_bfield.c
 * @brief HDF5 format bfield reading
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <hdf5.h>
#include "hdf5_hl.h"
#include "../math.h"
#include "../ascot5.h"
#include "../B_field.h"
#include "../Bfield/B_2DS.h"
#include "../Bfield/B_3DS.h"
#include "../Bfield/B_STS.h"
#include "hdf5_helpers.h"
#include "hdf5_bfield.h"

/**
 * @brief Initialize magnetic field offload data from h5 file
 *
 * This function reads the magnetic field data from the input.h5 file,
 * fills the offload struct with parameters and
 * allocates and fills the offload array.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
int hdf5_bfield_init_offload(hid_t f, B_field_offload_data* offload_data, real** offload_array) {
    herr_t err;

    #if VERBOSE > 0
        printf("\nReading magnetic field input from the HDF5 file...\n");
    #endif
    
    err = hdf5_find_group(f, "/bfield/");
    if(err < 0) {
        return -1;
    }
    
    char active[11];
    err = H5LTget_attribute_string(f, "/bfield/", "active", active);
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
	
    hdf5_generate_qid_path("/bfield/B_TC-XXXXXXXXXX", active, path);
    if(hdf5_find_group(f, path) == 0) {
	hdf5_bfield_init_offload_TC(f, &(offload_data->BTC), offload_array, active);
	offload_data->type = B_field_type_TC;
	offload_data->offload_array_length = offload_data->BTC.offload_array_length;
	
	#if VERBOSE > 0
	    printf("\nLoaded trivial cartesian magnetic field (B_TC)\n");
	    printf("with parameters:\n");
	    printf("- magnetic axis at (R,z) = (%le,%le)\n",
		   offload_data->BTC.axisr,offload_data->BTC.axisz);
	    printf("- psi = %le and rho %le\n",
		   offload_data->BTC.psival,offload_data->BTC.rhoval);
	    printf("- magnetic field at origo Bxyz = (%le, %le, %le)\n",
		   (*offload_array)[0],(*offload_array)[1],(*offload_array)[2]);
	    printf("- magnetic field gradient\n");
	    printf("  dBxdx = %le, dBxdy = %le, dBxdz = %le\n",
		   (*offload_array)[3],(*offload_array)[4],(*offload_array)[5]);
	    printf("  dBydx = %le, dBydy = %le, dBydz = %le\n",
		   (*offload_array)[6],(*offload_array)[7],(*offload_array)[8]);
	    printf("  dBzdx = %le, dBzdy = %le, dBzdz = %le\n",
		   (*offload_array)[9],(*offload_array)[10],(*offload_array)[11]);
	#endif

	return 1;
    }
    
    hdf5_generate_qid_path("/bfield/B_GS-XXXXXXXXXX", active, path);
    if(hdf5_find_group(f, path ) == 0) {
	hdf5_bfield_init_offload_GS(f, &(offload_data->BGS), offload_array, active);
	offload_data->type = B_field_type_GS;
	offload_data->offload_array_length = offload_data->BGS.offload_array_length;

	#if VERBOSE > 0
	    printf("\nLoaded analytical tokamak magnetic field (B_GS)\n");
	    printf("with parameters:\n");
	    printf("- magnetic axis at (R,z) = (%le,%le)\n",
		   offload_data->BGS.R0,offload_data->BGS.z0);
	    printf("- psi axis = %le, psi separatrix %le, and psi multiplier %le\n",
		   offload_data->BGS.psi0,offload_data->BGS.psi1,offload_data->BGS.psi_mult);
	    printf("- Toroidal field on-axis = %le and beta %le\n",
		   offload_data->BGS.B_phi0,(*offload_array)[12]);
	    for(int i=0; i<12; i++) {
	        printf("- C%d %le\n",i+1,(*offload_array)[i]);
	    }
	#endif
	return 1;
    }

    hdf5_generate_qid_path("/bfield/B_2DS-XXXXXXXXXX", active, path);
    if(hdf5_find_group(f, path) == 0) {
        hdf5_bfield_init_offload_2DS(f, &(offload_data->B2DS), offload_array, active);
	offload_data->type = B_field_type_2DS;
	offload_data->offload_array_length=offload_data->B2DS.offload_array_length;

	#if VERBOSE > 0
	    printf("\nLoaded 2D magnetic field (B_2DS)\n");
	    printf("with parameters:\n");
	    printf("- magnetic axis at (R,z) = (%le,%le)\n",
		   offload_data->B2DS.axis_r,offload_data->B2DS.axis_z);
	    printf("- psi axis = %le and psi separatrix %le\n",
		   offload_data->B2DS.psi0,offload_data->B2DS.psi1);
	    printf("- rmin, rmax, nr = %le, %le, %d\n",
		   offload_data->B2DS.r_min,offload_data->B2DS.r_max,offload_data->B2DS.n_r);
	    printf("- zmin, zmax, nz = %le, %le, %d\n",
		   offload_data->B2DS.z_min,offload_data->B2DS.z_max,offload_data->B2DS.n_z);
	#endif

        return 1;
    }
    
    hdf5_generate_qid_path("/bfield/B_3DS-XXXXXXXXXX", active, path);
    if(hdf5_find_group(f, path) == 0) {
        hdf5_bfield_init_offload_3DS(f, &(offload_data->B3DS), offload_array, active);
	offload_data->type = B_field_type_3DS;
	offload_data->offload_array_length=offload_data->B3DS.offload_array_length;

	#if VERBOSE > 0
	    printf("\nLoaded 3D magnetic field (B_3DS)\n");
	    printf("with parameters:\n");
	    printf("- magnetic axis at (R,z) = (%le,%le)\n",
		   offload_data->B3DS.axis_r,offload_data->B3DS.axis_z);
	    printf("- psi axis = %le and psi separatrix %le\n",
		   offload_data->B3DS.psi0,offload_data->B3DS.psi1);
	    printf("- rmin, rmax, nr = %le, %le, %d\n",
		   offload_data->B3DS.psigrid_r_min,offload_data->B3DS.psigrid_r_max,offload_data->B3DS.psigrid_n_r);
	    printf("- zmin, zmax, nz = %le, %le, %d\n",
		   offload_data->B3DS.psigrid_z_min,offload_data->B3DS.psigrid_z_max,offload_data->B3DS.psigrid_n_z);
	#endif
        return 1;
    }

    hdf5_generate_qid_path("/bfield/B_STS-XXXXXXXXXX", active, path);
    if(hdf5_find_group(f, path) == 0) {
        hdf5_bfield_init_offload_STS(f, &(offload_data->BSTS), offload_array, active);
        offload_data->type = B_field_type_STS;
        offload_data->offload_array_length=offload_data->BSTS.offload_array_length;
	#if VERBOSE > 0
	    printf("\nLoaded stellarator magnetic field (B_STS)\n");
	    printf("with parameters:\n");
	    printf("- number of toroidal periods = %d\n",
		   offload_data->BSTS.periods);
	    printf("- rmin, rmax, nr = %le, %le, %d\n",
		   offload_data->BSTS.r_min,offload_data->BSTS.r_max,offload_data->BSTS.n_r);
	    printf("- phimin, phimax, nphi = %le, %le, %d\n",
		   offload_data->BSTS.phi_min,offload_data->BSTS.phi_max,offload_data->BSTS.n_phi);
	    printf("- zmin, zmax, nz = %le, %le, %d\n",
		   offload_data->BSTS.z_min,offload_data->BSTS.z_max,offload_data->BSTS.n_z);
	#endif
        return 1;
    }
    
    printf("\nFailed to load magnetic field.\n");
    return -1;
}

void hdf5_bfield_init_offload_2DS(hid_t f, B_2DS_offload_data* offload_data, real** offload_array, char* qid) {
    herr_t err;
    char path[256];
    
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/bfield/B_2DS-XXXXXXXXXX/n_R", qid, path), &(offload_data->n_r));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/bfield/B_2DS-XXXXXXXXXX/n_z", qid, path), &(offload_data->n_z));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_2DS-XXXXXXXXXX/R_min", qid, path), &(offload_data->r_min));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_2DS-XXXXXXXXXX/R_max", qid, path), &(offload_data->r_max));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_2DS-XXXXXXXXXX/z_min", qid, path), &(offload_data->z_min));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_2DS-XXXXXXXXXX/z_max", qid, path), &(offload_data->z_max));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

    offload_data->r_grid = (offload_data->r_max - offload_data->r_min)
                           / (offload_data->n_r - 1);
    offload_data->z_grid = (offload_data->z_max - offload_data->z_min)
                           / (offload_data->n_z - 1);

    /* Allocate offload_array; psi and each component (r,phi,z) is
     * size n_r*n_z */
    int B_size = offload_data->n_r * offload_data->n_z;
    *offload_array = (real*) malloc(4 * B_size * sizeof(real));
    offload_data->offload_array_length = 4 * B_size;

    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_2DS-XXXXXXXXXX/psi", qid, path),   &(*offload_array)[0]);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_2DS-XXXXXXXXXX/B_R", qid, path),   &(*offload_array)[B_size]);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_2DS-XXXXXXXXXX/B_phi", qid, path), &(*offload_array)[2*B_size]);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_2DS-XXXXXXXXXX/B_z", qid, path),   &(*offload_array)[3*B_size]);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

    /* Read the first two values; These are the poloidal flux (psi) values at
     * magnetic axis and at x point (that is, separatrix). */
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_2DS-XXXXXXXXXX/psi0", qid, path), &(offload_data->psi0));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_2DS-XXXXXXXXXX/psi1", qid, path), &(offload_data->psi1));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

    /* Read magnetic axis r and z coordinates */
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_2DS-XXXXXXXXXX/axis_R", qid, path), &(offload_data->axis_r));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_2DS-XXXXXXXXXX/axis_z", qid, path), &(offload_data->axis_z));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
}

void hdf5_bfield_init_offload_3DS(hid_t f, B_3DS_offload_data* offload_data, real** offload_array, char* qid) {
    herr_t err;
    char path[256];

    /* Read and initialize magnetic field Rz-grid */
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/bfield/B_3DS-XXXXXXXXXX/n_R", qid, path), &(offload_data->Bgrid_n_r));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/bfield/B_3DS-XXXXXXXXXX/n_phi", qid, path), &(offload_data->n_phi));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/bfield/B_3DS-XXXXXXXXXX/n_z", qid, path), &(offload_data->Bgrid_n_z));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_3DS-XXXXXXXXXX/R_min", qid, path), &(offload_data->Bgrid_r_min));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_3DS-XXXXXXXXXX/R_max", qid, path), &(offload_data->Bgrid_r_max));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_3DS-XXXXXXXXXX/phi_min", qid, path), &(offload_data->phi_min));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_3DS-XXXXXXXXXX/phi_max", qid, path), &(offload_data->phi_max));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

    // Convert to radians
    offload_data->phi_min = offload_data->phi_min*(math_pi/180);
    offload_data->phi_max = offload_data->phi_max*(math_pi/180);
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_3DS-XXXXXXXXXX/z_min", qid, path), &(offload_data->Bgrid_z_min));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_3DS-XXXXXXXXXX/z_max", qid, path), &(offload_data->Bgrid_z_max));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

    offload_data->Bgrid_r_grid = (offload_data->Bgrid_r_max - offload_data->Bgrid_r_min)
                           / (offload_data->Bgrid_n_r - 1);
    offload_data->Bgrid_z_grid = (offload_data->Bgrid_z_max - offload_data->Bgrid_z_min)
                           / (offload_data->Bgrid_n_z - 1);
    offload_data->phi_grid = (offload_data->phi_max - offload_data->phi_min)
                           / (offload_data->n_phi);

    /* Read and initialize psi field Rz-grid */
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/bfield/B_3DS-XXXXXXXXXX/psigrid_n_R", qid, path), &(offload_data->psigrid_n_r));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/bfield/B_3DS-XXXXXXXXXX/psigrid_n_z", qid, path), &(offload_data->psigrid_n_z));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_3DS-XXXXXXXXXX/psigrid_R_min", qid, path), &(offload_data->psigrid_r_min));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_3DS-XXXXXXXXXX/psigrid_R_max", qid, path), &(offload_data->psigrid_r_max));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_3DS-XXXXXXXXXX/psigrid_z_min", qid, path), &(offload_data->psigrid_z_min));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_3DS-XXXXXXXXXX/psigrid_z_max", qid, path), &(offload_data->psigrid_z_max));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

    offload_data->psigrid_r_grid = (offload_data->psigrid_r_max - offload_data->psigrid_r_min)
                           / (offload_data->psigrid_n_r - 1);
    offload_data->psigrid_z_grid = (offload_data->psigrid_z_max - offload_data->psigrid_z_min)
                           / (offload_data->psigrid_n_z - 1);

    /* Allocate offload_array */
    int psi_size = offload_data->psigrid_n_r*offload_data->psigrid_n_z;
    int B_size = offload_data->Bgrid_n_r*offload_data->Bgrid_n_z*offload_data->n_phi;
    
    *offload_array = (real*) malloc((psi_size + 3 * B_size) * sizeof(real));
    offload_data->offload_array_length = psi_size + 3 * B_size;

    /* Read psi */
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_3DS-XXXXXXXXXX/psi", qid, path), &(*offload_array)[3*B_size]);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

    /* Read the magnetic field */
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_3DS-XXXXXXXXXX/B_R", qid, path), &(*offload_array)[0*B_size]);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_3DS-XXXXXXXXXX/B_phi", qid, path), &(*offload_array)[1*B_size]);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_3DS-XXXXXXXXXX/B_z", qid, path), &(*offload_array)[2*B_size]);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    /* Read the first two values; These are the poloidal flux (psi) values at
     * magnetic axis and at x point (that is, separatrix). */
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_3DS-XXXXXXXXXX/psi0", qid, path), &(offload_data->psi0));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_3DS-XXXXXXXXXX/psi1", qid, path), &(offload_data->psi1));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

    /* Read magnetic axis r and z coordinates */
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_3DS-XXXXXXXXXX/axis_R", qid, path), &(offload_data->axis_r));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_3DS-XXXXXXXXXX/axis_z", qid, path), &(offload_data->axis_z));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
}

/**
 * @brief Load stellarator magnetic field data from h5 file for splines
 *
 * This function reads the magnetic field data from the input.h5 file,
 * fills the offload struct with parameters and
 * allocates and fills the offload array.
 *
 * @todo Error checking
 * @todo Reading the magnetic axis
 * @todo Reading volume profiles
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
void hdf5_bfield_init_offload_STS(hid_t f, B_STS_offload_data* offload_data, real** offload_array, char* qid) {
    herr_t err;
    char path[256];
    int periods;

    /* Number of toroidal periods */
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/bfield/B_STS-XXXXXXXXXX/toroidalPeriods", qid, path), &periods);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    offload_data->periods = periods;

    /* Read the coordinate data */
 
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/bfield/B_STS-XXXXXXXXXX/n_r", qid, path), &(offload_data->n_r));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/bfield/B_STS-XXXXXXXXXX/n_phi", qid, path), &(offload_data->n_phi));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/bfield/B_STS-XXXXXXXXXX/n_z", qid, path), &(offload_data->n_z));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_STS-XXXXXXXXXX/r_min", qid, path), &(offload_data->r_min));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_STS-XXXXXXXXXX/r_max", qid, path), &(offload_data->r_max));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_STS-XXXXXXXXXX/z_min", qid, path), &(offload_data->z_min));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_STS-XXXXXXXXXX/z_max", qid, path), &(offload_data->z_max));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_STS-XXXXXXXXXX/phi_min", qid, path), &(offload_data->phi_min));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_STS-XXXXXXXXXX/phi_max", qid, path), &(offload_data->phi_max));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

    offload_data->r_grid = (offload_data->r_max - offload_data->r_min)
                           / (offload_data->n_r - 1);
    /* Note! phi data in input.h5 in deg */
    offload_data->phi_min = math_deg2rad(offload_data->phi_min);
    offload_data->phi_max = math_deg2rad(offload_data->phi_max);

    offload_data->phi_grid = (offload_data->phi_max - offload_data->phi_min)
                           / (offload_data->n_phi - 1);

    offload_data->z_grid = (offload_data->z_max - offload_data->z_min)
                           / (offload_data->n_z - 1);
    
    /* Read the bfield data */
    /* Allocate enough space for half a period */
    int temp_B_size = offload_data->n_r*offload_data->n_z*offload_data->n_phi;

    real* temp_B_r   = (real*) malloc(temp_B_size*sizeof(real));
    real* temp_B_phi = (real*) malloc(temp_B_size*sizeof(real));
    real* temp_B_z   = (real*) malloc(temp_B_size*sizeof(real));
    real* temp_B_s   = (real*) malloc(temp_B_size*sizeof(real));

    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_STS-XXXXXXXXXX/B_r", qid, path), temp_B_r);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_STS-XXXXXXXXXX/B_phi", qid, path), temp_B_phi);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_STS-XXXXXXXXXX/B_z", qid, path), temp_B_z);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_STS-XXXXXXXXXX/s", qid, path), temp_B_s);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

    /* We need to use stellarator symmetry here.
     * http://dx.doi.org/10.1016/S0167-2789(97)00216-9
     * The data is expected to include half a period.
     */
    int B_size = offload_data->n_r*offload_data->n_z*(2*(offload_data->n_phi - 1));
    *offload_array = (real*) malloc(4 * B_size * sizeof(real));
    offload_data->offload_array_length = 4 * B_size;
    
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
                 * temp_B_x data is in the format: (i_r, i_phi, i_z) = temp_B_x(i_z*n_phi*n_r + i_phi*n_r + i_r)
                 * offload_array data -"-"-"-"-  : (i_r, i_phi, i_z) = (*offload_array)[i_phi*n_z*n_r + i_z*n_r + i_r ]
                 * => (*offload_array)[i_phi*n_z*n_r + i_z*n_r + i_r ] = temp_B_x(i_z*n_phi*n_r + i_phi*n_r + i_r);
                 * The values are: Sym[B_r, B_phi, B_z] = [-B_r, B_phi, B_z]
                 */

                /* Index of data point in temp arrays */
                temp_ind = i_z*offload_data->n_phi*offload_data->n_r + i_phi*offload_data->n_r + i_r;

                /* Index of data point in offload_array and corresponding stel.-symmetric index  */
                off_ind = i_phi*offload_data->n_z*offload_data->n_r + i_z*offload_data->n_r + i_r;
                sym_ind = (2*(offload_data->n_phi - 1) - i_phi)*offload_data->n_z*offload_data->n_r
                    + (offload_data->n_z - i_z - 1)*offload_data->n_r + i_r;

                /* B_r */
                (*offload_array)[off_ind] =  temp_B_r[temp_ind];
                if (i_phi != 0) {
                    (*offload_array)[sym_ind] = -temp_B_r[temp_ind];
                }

                // B_phi
                (*offload_array)[B_size + off_ind] = temp_B_phi[temp_ind];
                if (i_phi != 0) {
                    (*offload_array)[B_size + sym_ind] = temp_B_phi[temp_ind];
                }
                // B_z
                (*offload_array)[2*B_size + off_ind] = temp_B_z[temp_ind];
                if (i_phi != 0) {
                    (*offload_array)[2*B_size + sym_ind] = temp_B_z[temp_ind];
                }
                // B_s
                (*offload_array)[3*B_size + off_ind] = temp_B_s[temp_ind];
                if (i_phi != 0) {
                    (*offload_array)[3*B_size + sym_ind] = temp_B_s[temp_ind];
                }                
            }
        }
    }
    /* Phi data is now for one toroidal period */
    offload_data->n_phi = 2*(offload_data->n_phi - 1);
    offload_data->phi_max = 2*offload_data->phi_max;
        
    /* Dummy values for magnetic axis */
    offload_data->axis_r = 0;
    offload_data->axis_z = 0;

    free(temp_B_r);
    free(temp_B_phi);
    free(temp_B_z);
    free(temp_B_s);
}

/**
 * @brief Reads Trivial Cartesian magnetic field from hdf5 file
 * 
 * @param f hdf5 source file
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
void hdf5_bfield_init_offload_TC(hid_t f, B_TC_offload_data* offload_data, real** offload_array, char* qid) {
    herr_t err;
    char path[256];

    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_TC-XXXXXXXXXX/axisr", qid, path), &(offload_data->axisr));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_TC-XXXXXXXXXX/axisz", qid, path), &(offload_data->axisr));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_TC-XXXXXXXXXX/psival", qid, path), &(offload_data->psival));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_TC-XXXXXXXXXX/rhoval", qid, path), &(offload_data->rhoval));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

    offload_data->offload_array_length = 12;

    *offload_array = (real*) malloc(offload_data->offload_array_length*sizeof(real));
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_TC-XXXXXXXXXX/Bxyz", qid, path), &(*offload_array)[0]);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_TC-XXXXXXXXXX/J", qid, path), &(*offload_array)[3]);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
}

/**
 * @brief Reads analytical magnetic field from hdf5 file
 * 
 * @param f hdf5 source file
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
void hdf5_bfield_init_offload_GS(hid_t f, B_GS_offload_data* offload_data, real** offload_array, char* qid) {
    herr_t err;
    char path[256];

    /* Equilibrium */
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_GS-XXXXXXXXXX/R0", qid, path), &(offload_data->R0));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_GS-XXXXXXXXXX/z0", qid, path), &(offload_data->z0));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_GS-XXXXXXXXXX/B_phi0", qid, path), &(offload_data->B_phi0));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_GS-XXXXXXXXXX/psi0", qid, path), &(offload_data->psi0));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_GS-XXXXXXXXXX/psi1", qid, path), &(offload_data->psi1));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_GS-XXXXXXXXXX/psi_mult", qid, path), &(offload_data->psi_mult));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

    /* Ripple */
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_GS-XXXXXXXXXX/delta0", qid, path), &(offload_data->delta0));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_GS-XXXXXXXXXX/alpha0", qid, path), &(offload_data->alpha0));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_GS-XXXXXXXXXX/a0", qid, path), &(offload_data->a0));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/bfield/B_GS-XXXXXXXXXX/Nripple", qid, path), &(offload_data->Nripple));
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

    offload_data->offload_array_length = 13;

    *offload_array = (real*) malloc(offload_data->offload_array_length*sizeof(real));
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_GS-XXXXXXXXXX/psi_coeff", qid, path), &(*offload_array)[0]);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
    
}
