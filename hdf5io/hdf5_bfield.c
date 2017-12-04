/**
 * @file hdf5_bfield.c
 * @brief HDF5 format bfield reading
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../math.h"
#include "../ascot5.h"
#include "../B_field.h"
#include "../Bfield/B_2D.h"
#include "../Bfield/B_2DS.h"
#include "../Bfield/B_3D.h"
#include "../Bfield/B_3DS.h"
#include "../Bfield/B_ST.h"
#include "hdf5.h"
#include "hdf5_helpers.h"
#include "hdf5_hl.h"
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
    int version;
    int periods;
    hsize_t dims[3];

    err = H5LTfind_dataset(f, "/bfield/");
    if(err < 0) {
        return -1;
    }
    char type[32];

    err = H5LTget_attribute_string(f, "/bfield/", "type", type);
    if(err < 0) {
        return -1;
    }

    if(strncmp(type,"B_TC",4) == 0) {
	hdf5_bfield_init_offload_TC(f, &(offload_data->BTC), offload_array);
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
    else if(strncmp(type,"B_GS",4) == 0) {
	hdf5_bfield_init_offload_GS(f, &(offload_data->BGS), offload_array);
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
    else if(strncmp(type,"B_2D",4) == 0) {
        hdf5_bfield_init_offload_2D(f, &(offload_data->B2DS), offload_array);
	offload_data->type = B_field_type_2DS;
    offload_data->offload_array_length=offload_data->B2DS.offload_array_length;

	#if VERBOSE > 0
	    printf("\nLoaded 2D magnetic field (B_2D)\n");
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
    else if (strncmp(type, "B_3D",4) == 0) {
        hdf5_bfield_init_offload_3DS(f, &(offload_data->B3DS), offload_array);
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
		   offload_data->B3DS.r_min,offload_data->B3DS.r_max,offload_data->B3DS.n_r);
	    printf("- zmin, zmax, nz = %le, %le, %d\n",
		   offload_data->B3DS.z_min,offload_data->B3DS.z_max,offload_data->B3DS.n_z);
	#endif
        return 1;
    }
    else if (strncmp(type, "B_ST",4) == 0) {
        hdf5_bfield_init_offload_STS(f, &(offload_data->BSTS), offload_array);
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
    
    printf("\nFailed to load magnetic field\n");
    return -1;
}

void hdf5_bfield_init_offload_2D(hid_t f, B_2DS_offload_data* offload_data, real** offload_array) {
    herr_t err;
        
    err = H5LTread_dataset_int(f,"/bfield/B_2D/n_r",&(offload_data->n_r));
    err = H5LTread_dataset_int(f,"/bfield/B_2D/n_z",&(offload_data->n_z));
    err = H5LTread_dataset_double(f,"/bfield/B_2D/r_min",&(offload_data->r_min));
    err = H5LTread_dataset_double(f,"/bfield/B_2D/r_max",&(offload_data->r_max));
    err = H5LTread_dataset_double(f,"/bfield/B_2D/z_min",&(offload_data->z_min));
    err = H5LTread_dataset_double(f,"/bfield/B_2D/z_max",&(offload_data->z_max));

    offload_data->r_grid = (offload_data->r_max - offload_data->r_min)
                           / (offload_data->n_r - 1);
    offload_data->z_grid = (offload_data->z_max - offload_data->z_min)
                           / (offload_data->n_z - 1);

    /* Allocate offload_array; psi and each component (r,phi,z) is
     * size n_r*n_z */
    int B_size = offload_data->n_r * offload_data->n_z;
    *offload_array = (real*) malloc(4 * B_size * sizeof(real));
    offload_data->offload_array_length = 4 * B_size;

    err = H5LTread_dataset_double(f,"/bfield/B_2D/psi",&(*offload_array)[0]);
    err = H5LTread_dataset_double(f,"/bfield/B_2D/B_r",&(*offload_array)[B_size]);
    err = H5LTread_dataset_double(f,"/bfield/B_2D/B_phi",&(*offload_array)[2*B_size]);
    err = H5LTread_dataset_double(f,"/bfield/B_2D/B_z",&(*offload_array)[3*B_size]);

    /* Read the first two values; These are the poloidal flux (psi) values at
     * magnetic axis and at x point (that is, separatrix). */
    err = H5LTread_dataset_double(f,"/bfield/B_2D/psi0",&(offload_data->psi0));
    err = H5LTread_dataset_double(f,"/bfield/B_2D/psi1",&(offload_data->psi1));

    /* Read magnetic axis r and z coordinates */
    err = H5LTread_dataset_double(f,"/bfield/B_2D/axis_r",&(offload_data->axis_r));
    err = H5LTread_dataset_double(f,"/bfield/B_2D/axis_z",&(offload_data->axis_z));
}



void hdf5_bfield_init_offload_3D(hid_t f, B_3D_offload_data* offload_data, real** offload_array) {
    herr_t err;
    
    err = H5LTread_dataset_int(f,"/bfield/B_3D/n_r",&(offload_data->n_r));
    err = H5LTread_dataset_int(f,"/bfield/B_3D/n_phi",&(offload_data->n_phi));
    err = H5LTread_dataset_int(f,"/bfield/B_3D/n_z",&(offload_data->n_z));
    err = H5LTread_dataset_double(f,"/bfield/B_3D/r_min",&(offload_data->r_min));
    err = H5LTread_dataset_double(f,"/bfield/B_3D/r_max",&(offload_data->r_max));
    err = H5LTread_dataset_double(f,"/bfield/B_3D/z_min",&(offload_data->z_min));
    err = H5LTread_dataset_double(f,"/bfield/B_3D/z_max",&(offload_data->z_max));

    offload_data->r_grid = (offload_data->r_max - offload_data->r_min)
                           / (offload_data->n_r - 1);
    offload_data->z_grid = (offload_data->z_max - offload_data->z_min)
                           / (offload_data->n_z - 1);
    offload_data->phi_grid = 2*math_pi / (offload_data->n_phi);

    /* phi array starts from -0.5*phi_grid and ends at 0.5*phi_grid + 2pi! */
    offload_data->phi_min = -1.5 * offload_data->phi_grid;
    offload_data->phi_max = 1.5 * offload_data->phi_grid + 2*math_pi;

    /* Allocate offload_array */
    int psi_size = offload_data->n_r*offload_data->n_z;
    int B_size = offload_data->n_r*offload_data->n_z*(offload_data->n_phi+4);
    *offload_array = (real*) malloc((psi_size + 3 * B_size) * sizeof(real));
    offload_data->offload_array_length = psi_size + 3 * B_size;

    /* Read psi */
    err = H5LTread_dataset_double(f,"/bfield/B_3D/psi",&(*offload_array)[3*B_size]);

    /* Calculate 2D components of poloidal field */
    real* eq_B_r = (real*) malloc(psi_size*sizeof(real));
    real* eq_B_z = (real*) malloc(psi_size*sizeof(real));

    int i;
    for(i = 0; i < offload_data->n_r; i++) {
        int j;
        for(j = 0; j < offload_data->n_z; j++) {
            real psi[4];
            B_2D_bicubic_derivs(psi, 0, 0, i, j, offload_data->n_r, offload_data->r_grid, offload_data->z_grid, &(*offload_array)[3*B_size]);
            eq_B_r[j*offload_data->n_r + i] = 1/(2*math_pi) * -psi[3] / (offload_data->r_min + i*offload_data->r_grid);
            eq_B_z[j*offload_data->n_r + i] = 1/(2*math_pi) * psi[1] / (offload_data->r_min + i*offload_data->r_grid);
        }
    }

    /* Read the magnetic field */
    real* temp_B_r = (real*) malloc(psi_size*offload_data->n_phi*sizeof(real));
    real* temp_B_phi = (real*) malloc(psi_size*offload_data->n_phi*sizeof(real));
    real* temp_B_z = (real*) malloc(psi_size*offload_data->n_phi*sizeof(real));

    err = H5LTread_dataset_double(f,"/bfield/B_3D/B_r",temp_B_r);
    err = H5LTread_dataset_double(f,"/bfield/B_3D/B_phi",temp_B_phi);
    err = H5LTread_dataset_double(f,"/bfield/B_3D/B_z",temp_B_z);

    /* permute the phi and z dimensions */
    for(i = 0; i < offload_data->n_phi; i++) {
        int j;
        for(j = 0; j < offload_data->n_z; j++) {
            int k;
            for(k = 0; k < offload_data->n_r; k++) {
               (*offload_array)[2*psi_size+i*offload_data->n_z*offload_data->n_r
                          + j*offload_data->n_r + k] =
                        eq_B_r[j*offload_data->n_r + k] +
                        temp_B_r[j*offload_data->n_phi*offload_data->n_r
                          + i*offload_data->n_r + k];
               (*offload_array)[2*psi_size+i*offload_data->n_z*offload_data->n_r
                          + j*offload_data->n_r + k + B_size] =
                        temp_B_phi[j*offload_data->n_phi*offload_data->n_r
                          + i*offload_data->n_r + k];
               (*offload_array)[2*psi_size+i*offload_data->n_z*offload_data->n_r
                          + j*offload_data->n_r + k + 2*B_size] =
                        eq_B_z[j*offload_data->n_r + k] +
                        temp_B_z[j*offload_data->n_phi*offload_data->n_r
                          + i*offload_data->n_r + k];
            }
        }
    }
    
    free(eq_B_r);
    free(eq_B_z);
 
    /* Copy two phi slices into opposite ends for each field
     * component */
    memcpy(&(*offload_array)[0],
        &(*offload_array)[offload_data->n_phi*psi_size],
        psi_size * sizeof(real));
    memcpy(&(*offload_array)[B_size],
        &(*offload_array)[offload_data->n_phi*psi_size+B_size],
        psi_size * sizeof(real));
    memcpy(&(*offload_array)[2*B_size],
        &(*offload_array)[offload_data->n_phi*psi_size+2*B_size],
        psi_size * sizeof(real));

    memcpy(&(*offload_array)[psi_size],
        &(*offload_array)[(offload_data->n_phi+1)*psi_size],
        psi_size * sizeof(real));
    memcpy(&(*offload_array)[psi_size + B_size],
        &(*offload_array)[(offload_data->n_phi+1)*psi_size+B_size],
        psi_size * sizeof(real));
    memcpy(&(*offload_array)[psi_size + 2*B_size],
        &(*offload_array)[(offload_data->n_phi+1)*psi_size+2*B_size],
        psi_size * sizeof(real));

    memcpy(&(*offload_array)[(offload_data->n_phi+2)*psi_size],
        &(*offload_array)[2*psi_size],
        psi_size * sizeof(real));
    memcpy(&(*offload_array)[(offload_data->n_phi+2)*psi_size+B_size],
        &(*offload_array)[2*psi_size+B_size],
        psi_size * sizeof(real));
    memcpy(&(*offload_array)[(offload_data->n_phi+2)*psi_size+2*B_size],
        &(*offload_array)[2*psi_size+2*B_size],
        psi_size * sizeof(real));

    memcpy(&(*offload_array)[(offload_data->n_phi+3)*psi_size],
        &(*offload_array)[3*psi_size],
        psi_size * sizeof(real));
    memcpy(&(*offload_array)[(offload_data->n_phi+3)*psi_size+B_size],
        &(*offload_array)[3*psi_size+B_size],
        psi_size * sizeof(real));
    memcpy(&(*offload_array)[(offload_data->n_phi+3)*psi_size+2*B_size],
        &(*offload_array)[3*psi_size+2*B_size],
        psi_size * sizeof(real));

    /* Read the first two values; These are the poloidal flux (psi) values at
     * magnetic axis and at x point (that is, separatrix). */
    err = H5LTread_dataset_double(f,"/bfield/B_3D/psi0",&(offload_data->psi0));
    err = H5LTread_dataset_double(f,"/bfield/B_3D/psi1",&(offload_data->psi1));

    /* Read magnetic axis r and z coordinates */
    err = H5LTread_dataset_double(f,"/bfield/B_3D/axis_r",&(offload_data->axis_r));
    err = H5LTread_dataset_double(f,"/bfield/B_3D/axis_z",&(offload_data->axis_z));
}

void hdf5_bfield_init_offload_3DS(hid_t f, B_3DS_offload_data* offload_data, real** offload_array) {
    herr_t err;
    
    err = H5LTread_dataset_int(f,"/bfield/B_3D/n_r",&(offload_data->n_r));
    err = H5LTread_dataset_int(f,"/bfield/B_3D/n_phi",&(offload_data->n_phi));
    err = H5LTread_dataset_int(f,"/bfield/B_3D/n_z",&(offload_data->n_z));
    err = H5LTread_dataset_double(f,"/bfield/B_3D/r_min",&(offload_data->r_min));
    err = H5LTread_dataset_double(f,"/bfield/B_3D/r_max",&(offload_data->r_max));
    err = H5LTread_dataset_double(f,"/bfield/B_3D/z_min",&(offload_data->z_min));
    err = H5LTread_dataset_double(f,"/bfield/B_3D/z_max",&(offload_data->z_max));

    offload_data->r_grid = (offload_data->r_max - offload_data->r_min)
                           / (offload_data->n_r - 1);
    offload_data->z_grid = (offload_data->z_max - offload_data->z_min)
                           / (offload_data->n_z - 1);
    //offload_data->phi_grid = 2*math_pi / (offload_data->n_phi);

    /* phi array starts from -0.5*phi_grid and ends at 0.5*phi_grid + 2pi! */
    //offload_data->phi_min = -1.5 * offload_data->phi_grid;
    //offload_data->phi_max = 1.5 * offload_data->phi_grid + 2*math_pi;

    offload_data->phi_grid = 2*math_pi / (offload_data->n_phi);
    offload_data->phi_min = 0;
    offload_data->phi_max = 2*math_pi;

    /* Allocate offload_array */
    int psi_size = offload_data->n_r*offload_data->n_z;
    int B_size = offload_data->n_r*offload_data->n_z*offload_data->n_phi;
    *offload_array = (real*) malloc((psi_size + 3 * B_size) * sizeof(real));
    offload_data->offload_array_length = psi_size + 3 * B_size;

    /* Read psi */
    err = H5LTread_dataset_double(f,"/bfield/B_3D/psi",&(*offload_array)[3*B_size]);

    /* Read the magnetic field */
    real* temp_B_r = (real*) malloc(psi_size*offload_data->n_phi*sizeof(real));
    real* temp_B_phi = (real*) malloc(psi_size*offload_data->n_phi*sizeof(real));
    real* temp_B_z = (real*) malloc(psi_size*offload_data->n_phi*sizeof(real));

    err = H5LTread_dataset_double(f,"/bfield/B_3D/B_r",temp_B_r);
    err = H5LTread_dataset_double(f,"/bfield/B_3D/B_phi",temp_B_phi);
    err = H5LTread_dataset_double(f,"/bfield/B_3D/B_z",temp_B_z);

    /* permute the phi and z dimensions */
    int i;
    for(i = 0; i < offload_data->n_phi; i++) {
        int j;
        for(j = 0; j < offload_data->n_z; j++) {
            int k;
            for(k = 0; k < offload_data->n_r; k++) {
               (*offload_array)[i*offload_data->n_z*offload_data->n_r
                          + j*offload_data->n_r + k] =
                        temp_B_r[j*offload_data->n_phi*offload_data->n_r
                          + i*offload_data->n_r + k];
               (*offload_array)[i*offload_data->n_z*offload_data->n_r
                          + j*offload_data->n_r + k + B_size] =
                        temp_B_phi[j*offload_data->n_phi*offload_data->n_r
                          + i*offload_data->n_r + k];
               (*offload_array)[i*offload_data->n_z*offload_data->n_r
                          + j*offload_data->n_r + k + 2*B_size] =
                        temp_B_z[j*offload_data->n_phi*offload_data->n_r
                          + i*offload_data->n_r + k];
            }
        }
    }

    /* Read the first two values; These are the poloidal flux (psi) values at
     * magnetic axis and at x point (that is, separatrix). */
    err = H5LTread_dataset_double(f,"/bfield/B_3D/psi0",&(offload_data->psi0));
    err = H5LTread_dataset_double(f,"/bfield/B_3D/psi1",&(offload_data->psi1));

    /* Read magnetic axis r and z coordinates */
    err = H5LTread_dataset_double(f,"/bfield/B_3D/axis_r",&(offload_data->axis_r));
    err = H5LTread_dataset_double(f,"/bfield/B_3D/axis_z",&(offload_data->axis_z));

    free(temp_B_r);
    free(temp_B_phi);
    free(temp_B_z);
}

/**
 * @brief Load stellarator magnetic field data from h5 file
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
void hdf5_bfield_init_offload_ST(hid_t f, B_ST_offload_data* offload_data, real** offload_array) {
    herr_t err;
    int periods;
    hsize_t dims[3];

    /* Number of toroidal periods */
    err = H5LTread_dataset_int(f,"/bfield/B_ST/toroidalPeriods",&periods);
    offload_data->periods = periods;

    /* Read the coordinate data */
 
    err = H5LTread_dataset_int(f,"/bfield/B_ST/n_r",&(offload_data->n_r));
    err = H5LTread_dataset_int(f,"/bfield/B_ST/n_phi",&(offload_data->n_phi));
    err = H5LTread_dataset_int(f,"/bfield/B_ST/n_z",&(offload_data->n_z));
    err = H5LTread_dataset_double(f,"/bfield/B_ST/r_min",&(offload_data->r_min));
    err = H5LTread_dataset_double(f,"/bfield/B_ST/r_max",&(offload_data->r_max));
    err = H5LTread_dataset_double(f,"/bfield/B_ST/z_min",&(offload_data->z_min));
    err = H5LTread_dataset_double(f,"/bfield/B_ST/z_max",&(offload_data->z_max));
    err = H5LTread_dataset_double(f,"/bfield/B_ST/phi_min",&(offload_data->phi_min));
    err = H5LTread_dataset_double(f,"/bfield/B_ST/phi_max",&(offload_data->phi_max));

    offload_data->r_grid = (offload_data->r_max - offload_data->r_min)
                           / (offload_data->n_r - 1);
    /* Note! phi data in input.h5 in deg */
    offload_data->phi_min = offload_data->phi_min / 360.0 * 2.0 *math_pi;
    offload_data->phi_max = offload_data->phi_max / 360.0 * 2.0 *math_pi;

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

    err = H5LTread_dataset_double(f, "/bfield/B_ST/B_r", temp_B_r);
    err = H5LTread_dataset_double(f, "/bfield/B_ST/B_phi", temp_B_phi);
    err = H5LTread_dataset_double(f, "/bfield/B_ST/B_z", temp_B_z);
    err = H5LTread_dataset_double(f, "/bfield/B_ST/s", temp_B_s);

    /* We need to use stellarator symmetry here.
     * http://dx.doi.org/10.1016/S0167-2789(97)00216-9
     * The data is expected to include half a period.
     */
    int B_size = offload_data->n_r*offload_data->n_z*(2*(offload_data->n_phi - 1) + 1 + 3);
    int phi_size = offload_data->n_r*offload_data->n_z;
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
                off_ind = 2*phi_size + i_phi*offload_data->n_z*offload_data->n_r + i_z*offload_data->n_r + i_r;
                sym_ind = 2*phi_size + (2*(offload_data->n_phi - 1) - i_phi)*offload_data->n_z*offload_data->n_r
                    + (offload_data->n_z - i_z - 1)*offload_data->n_r + i_r;

                /* B_r */
                (*offload_array)[off_ind] =  temp_B_r[temp_ind];
                (*offload_array)[sym_ind] = -temp_B_r[temp_ind];

                // B_phi
                (*offload_array)[B_size + off_ind] = temp_B_phi[temp_ind];
                (*offload_array)[B_size + sym_ind] = temp_B_phi[temp_ind];

                // B_z
                (*offload_array)[2*B_size + off_ind] = temp_B_z[temp_ind];
                (*offload_array)[2*B_size + sym_ind] = temp_B_z[temp_ind];

                // B_s
                (*offload_array)[3*B_size + off_ind] = temp_B_s[temp_ind];
                (*offload_array)[3*B_size + sym_ind] = temp_B_s[temp_ind];
            }
        }
    }

    /* Phi data is now for one toroidal period */
    offload_data->n_phi = 2*(offload_data->n_phi - 1);
    offload_data->phi_max = offload_data->phi_min + (offload_data->n_phi - 1)*offload_data->phi_grid;

    /* Copy two phi slices into opposite ends for each field
     * component */
    memcpy(&(*offload_array)[0],
        &(*offload_array)[offload_data->n_phi*phi_size],
        phi_size * sizeof(real));
    memcpy(&(*offload_array)[B_size],
        &(*offload_array)[B_size + offload_data->n_phi*phi_size],
        phi_size * sizeof(real));
    memcpy(&(*offload_array)[2*B_size],
        &(*offload_array)[2*B_size + offload_data->n_phi*phi_size],
        phi_size * sizeof(real));
    memcpy(&(*offload_array)[3*B_size],
        &(*offload_array)[3*B_size + offload_data->n_phi*phi_size],
        phi_size * sizeof(real));

    memcpy(&(*offload_array)[phi_size],
        &(*offload_array)[(offload_data->n_phi+1)*phi_size],
        phi_size * sizeof(real));
    memcpy(&(*offload_array)[B_size + phi_size],
        &(*offload_array)[B_size + (offload_data->n_phi+1)*phi_size],
        phi_size * sizeof(real));
    memcpy(&(*offload_array)[2*B_size + phi_size],
        &(*offload_array)[2*B_size + (offload_data->n_phi+1)*phi_size],
        phi_size * sizeof(real));
    memcpy(&(*offload_array)[3*B_size + phi_size],
        &(*offload_array)[3*B_size + (offload_data->n_phi+1)*phi_size],
        phi_size * sizeof(real));

    memcpy(&(*offload_array)[(offload_data->n_phi+2)*phi_size],
        &(*offload_array)[2*phi_size],
        phi_size * sizeof(real));
    memcpy(&(*offload_array)[B_size + (offload_data->n_phi+2)*phi_size],
        &(*offload_array)[B_size + 2*phi_size],
        phi_size * sizeof(real));
    memcpy(&(*offload_array)[2*B_size + (offload_data->n_phi+2)*phi_size],
        &(*offload_array)[2*B_size + 2*phi_size],
        phi_size * sizeof(real));
    memcpy(&(*offload_array)[3*B_size + (offload_data->n_phi+2)*phi_size],
        &(*offload_array)[3*B_size + 2*phi_size],
        phi_size * sizeof(real));

    memcpy(&(*offload_array)[(offload_data->n_phi+3)*phi_size],
        &(*offload_array)[3*phi_size],
        phi_size * sizeof(real));
    memcpy(&(*offload_array)[B_size + (offload_data->n_phi+3)*phi_size],
        &(*offload_array)[B_size + 3*phi_size],
        phi_size * sizeof(real));
    memcpy(&(*offload_array)[2*B_size + (offload_data->n_phi+3)*phi_size],
        &(*offload_array)[2*B_size + 3*phi_size],
        phi_size * sizeof(real));
    memcpy(&(*offload_array)[3*B_size + (offload_data->n_phi+3)*phi_size],
        &(*offload_array)[3*B_size + 3*phi_size],
        phi_size * sizeof(real));

    /* Copying the segments changes the phi limits, but not n_phi */
    offload_data->phi_min = offload_data->phi_min - 2*offload_data->phi_grid;
    offload_data->phi_max = offload_data->phi_max + 2*offload_data->phi_grid;

    /* Dummy values for magnetic axis */
    offload_data->axis_r = 0;
    offload_data->axis_z = 0;

    free(temp_B_r);
    free(temp_B_phi);
    free(temp_B_z);
    free(temp_B_s);
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
void hdf5_bfield_init_offload_STS(hid_t f, B_STS_offload_data* offload_data, real** offload_array) {
    herr_t err;
    int periods;
    hsize_t dims[3];

    /* Number of toroidal periods */
    err = H5LTread_dataset_int(f,"/bfield/B_ST/toroidalPeriods",&periods);
    offload_data->periods = periods;

    /* Read the coordinate data */
 
    err = H5LTread_dataset_int(f,"/bfield/B_ST/n_r",&(offload_data->n_r));
    err = H5LTread_dataset_int(f,"/bfield/B_ST/n_phi",&(offload_data->n_phi));
    err = H5LTread_dataset_int(f,"/bfield/B_ST/n_z",&(offload_data->n_z));
    err = H5LTread_dataset_double(f,"/bfield/B_ST/r_min",&(offload_data->r_min));
    err = H5LTread_dataset_double(f,"/bfield/B_ST/r_max",&(offload_data->r_max));
    err = H5LTread_dataset_double(f,"/bfield/B_ST/z_min",&(offload_data->z_min));
    err = H5LTread_dataset_double(f,"/bfield/B_ST/z_max",&(offload_data->z_max));
    err = H5LTread_dataset_double(f,"/bfield/B_ST/phi_min",&(offload_data->phi_min));
    err = H5LTread_dataset_double(f,"/bfield/B_ST/phi_max",&(offload_data->phi_max));

    offload_data->r_grid = (offload_data->r_max - offload_data->r_min)
                           / (offload_data->n_r - 1);
    /* Note! phi data in input.h5 in deg */
    offload_data->phi_min = offload_data->phi_min / 360.0 * 2.0 *math_pi;
    offload_data->phi_max = offload_data->phi_max / 360.0 * 2.0 *math_pi;

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

    err = H5LTread_dataset_double(f, "/bfield/B_ST/B_r", temp_B_r);
    err = H5LTread_dataset_double(f, "/bfield/B_ST/B_phi", temp_B_phi);
    err = H5LTread_dataset_double(f, "/bfield/B_ST/B_z", temp_B_z);
    err = H5LTread_dataset_double(f, "/bfield/B_ST/s", temp_B_s);

    /* We need to use stellarator symmetry here.
     * http://dx.doi.org/10.1016/S0167-2789(97)00216-9
     * The data is expected to include half a period.
     */
    int B_size = offload_data->n_r*offload_data->n_z*(2*(offload_data->n_phi - 1) + 1 + 3);
    int phi_size = offload_data->n_r*offload_data->n_z;
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
                off_ind = 2*phi_size + i_phi*offload_data->n_z*offload_data->n_r + i_z*offload_data->n_r + i_r;
                sym_ind = 2*phi_size + (2*(offload_data->n_phi - 1) - i_phi)*offload_data->n_z*offload_data->n_r
                    + (offload_data->n_z - i_z - 1)*offload_data->n_r + i_r;

                /* B_r */
                (*offload_array)[off_ind] =  temp_B_r[temp_ind];
                (*offload_array)[sym_ind] = -temp_B_r[temp_ind];

                // B_phi
                (*offload_array)[B_size + off_ind] = temp_B_phi[temp_ind];
                (*offload_array)[B_size + sym_ind] = temp_B_phi[temp_ind];

                // B_z
                (*offload_array)[2*B_size + off_ind] = temp_B_z[temp_ind];
                (*offload_array)[2*B_size + sym_ind] = temp_B_z[temp_ind];

                // B_s
                (*offload_array)[3*B_size + off_ind] = temp_B_s[temp_ind];
                (*offload_array)[3*B_size + sym_ind] = temp_B_s[temp_ind];
            }
        }
    }

    /* Phi data is now for one toroidal period */
    offload_data->n_phi = 2*(offload_data->n_phi - 1);
    offload_data->phi_max = offload_data->phi_min + (offload_data->n_phi - 1)*offload_data->phi_grid;

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
void hdf5_bfield_init_offload_TC(hid_t f, B_TC_offload_data* offload_data, real** offload_array) {
    herr_t err;

    err = H5LTread_dataset_double(f,"/bfield/B_TC/axisr",&(offload_data->axisr));
    err = H5LTread_dataset_double(f,"/bfield/B_TC/axisz",&(offload_data->axisr));
    err = H5LTread_dataset_double(f,"/bfield/B_TC/psival",&(offload_data->psival));
    err = H5LTread_dataset_double(f,"/bfield/B_TC/rhoval",&(offload_data->rhoval));

    offload_data->offload_array_length = 12;

    *offload_array = (real*) malloc(offload_data->offload_array_length*sizeof(real));
    err = H5LTread_dataset_double(f,"/bfield/B_TC/Bxyz",&(*offload_array)[0]);
    err = H5LTread_dataset_double(f,"/bfield/B_TC/gradB",&(*offload_array)[3]);
    

}

/**
 * @brief Reads analytical magnetic field from hdf5 file
 * 
 * @param f hdf5 source file
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
void hdf5_bfield_init_offload_GS(hid_t f, B_GS_offload_data* offload_data, real** offload_array) {
    herr_t err;

    /* Equilibrium */
    err = H5LTread_dataset_double(f,"/bfield/B_GS/R0",&(offload_data->R0));
    err = H5LTread_dataset_double(f,"/bfield/B_GS/z0",&(offload_data->z0));
    err = H5LTread_dataset_double(f,"/bfield/B_GS/B_phi0",&(offload_data->B_phi0));
    err = H5LTread_dataset_double(f,"/bfield/B_GS/psi0",&(offload_data->psi0));
    err = H5LTread_dataset_double(f,"/bfield/B_GS/psi1",&(offload_data->psi1));
    err = H5LTread_dataset_double(f,"/bfield/B_GS/psi_mult",&(offload_data->psi_mult));

    /* Ripple */
    err = H5LTread_dataset_double(f,"/bfield/B_GS/delta0",&(offload_data->delta0));
    err = H5LTread_dataset_double(f,"/bfield/B_GS/alpha0",&(offload_data->alpha0));
    err = H5LTread_dataset_double(f,"/bfield/B_GS/a0",&(offload_data->a0));
    err = H5LTread_dataset_int(f,"/bfield/B_GS/Nripple",&(offload_data->Nripple));

    offload_data->offload_array_length = 13;

    *offload_array = (real*) malloc(offload_data->offload_array_length*sizeof(real));
    err = H5LTread_dataset_double(f,"/bfield/B_GS/psi_coeff",&(*offload_array)[0]);
    
}
