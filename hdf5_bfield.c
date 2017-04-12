/**
 * @file hdf5_bfield.c
 * @brief HDF5 format bfield reading
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "math.h"
#include "ascot5.h"
#include "B_ST.h"
#include "hdf5.h"
#include "hdf5_helpers.h"
#include "hdf5_hl.h"

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
void hdf5_bfield_init_offload_ST(B_ST_offload_data* offload_data, real** offload_array) {
    herr_t err;
    int version;
    int periods;
    hsize_t dims[3];

    hid_t f = hdf5_open("input.h5");

    err = H5LTfind_dataset(f, "/bfield/");
    err = H5LTfind_dataset(f, "/bfield/stellarator");

    /* */
    err = H5LTget_attribute_int(f, "/bfield/stellarator", "version", &version);

    /* Number of toroidal periods */
    err = H5LTread_dataset_int(f,"/bfield/stellarator/toroidalPeriods",&periods);
    offload_data->periods = periods;
    
    /* Read the coordinate data */
    /* Get dimensions of bfield data */
    err = H5LTget_dataset_info(f,"/bfield/stellarator/br", dims, NULL, NULL);

    /* Sizes of dimensions */
    int n_z   = dims[0];
    int n_phi = dims[1];
    int n_r   = dims[2];
    offload_data->n_z   = n_z;
    offload_data->n_phi = n_phi;
    offload_data->n_r   = n_r;

    real* temp_r   = (real*) malloc(offload_data->n_r*sizeof(real));
    real* temp_phi = (real*) malloc(offload_data->n_phi*sizeof(real));
    real* temp_z   = (real*) malloc(offload_data->n_z*sizeof(real));

    err = H5LTread_dataset_double(f, "/bfield/stellarator/r", temp_r);
    err = H5LTread_dataset_double(f, "/bfield/stellarator/phi", temp_phi);
    err = H5LTread_dataset_double(f, "/bfield/stellarator/z", temp_z);

    offload_data->r_min = temp_r[0];
    offload_data->r_max = temp_r[n_r - 1];
    offload_data->r_grid = (temp_r[1] - temp_r[0]);
    /* Note! phi data in input.h5 in deg */
    offload_data->phi_min = temp_phi[0] / 360 * 2*math_pi;
    offload_data->phi_max = temp_phi[n_phi - 1] / 360 * 2*math_pi;
    offload_data->phi_grid = (temp_phi[1] - temp_phi[0]) / 360 * 2*math_pi;

    offload_data->z_min = temp_z[0];
    offload_data->z_max = temp_z[n_z - 1];
    offload_data->z_grid = (temp_z[1] - temp_z[0]);

    free(temp_r);
    free(temp_phi);
    free(temp_z);

    /* Read the bfield data */
    /* Allocate enough space for half a period */
    int temp_B_size = n_r*n_z*n_phi;

    real* temp_B_r   = (real*) malloc(temp_B_size*sizeof(real));
    real* temp_B_phi = (real*) malloc(temp_B_size*sizeof(real));
    real* temp_B_z   = (real*) malloc(temp_B_size*sizeof(real));
    real* temp_B_s   = (real*) malloc(temp_B_size*sizeof(real));

    err = H5LTread_dataset_double(f, "/bfield/stellarator/br", temp_B_r);
    err = H5LTread_dataset_double(f, "/bfield/stellarator/bphi", temp_B_phi);
    err = H5LTread_dataset_double(f, "/bfield/stellarator/bz", temp_B_z);
    err = H5LTread_dataset_double(f, "/bfield/stellarator/s", temp_B_s);

    /* We need to use stellarator symmetry here.
     * http://dx.doi.org/10.1016/S0167-2789(97)00216-9
     * The data is expected to include half a period.
     */

    int B_size = n_r*n_z*(2*(n_phi - 1) + 1 + 3);
    int phi_size = n_r*n_z;
    *offload_array = (real*) malloc(4 * B_size * sizeof(real));
    offload_data->offload_array_length = 4 * B_size;

    int i_phi;
    int i_z;
    int i_r;
    int temp_ind, off_ind, sym_ind;
    for (i_phi = 0; i_phi < n_phi; i_phi++) {
        for (i_z = 0; i_z < n_z; i_z++) {
            for (i_r = 0; i_r < n_r; i_r++) {
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
                temp_ind = i_z*n_phi*n_r + i_phi*n_r + i_r;
                
                /* Index of data point in offload_array and corresponding stel.-symmetric index  */
                off_ind = 2*phi_size + i_phi*n_z*n_r + i_z*n_r + i_r;
                sym_ind = 2*phi_size + (2*(n_phi - 1) - i_phi)*n_z*n_r + (n_z - i_z - 1)*n_r + i_r;
                
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
    n_phi = 2*(n_phi - 1);
    offload_data->n_phi = n_phi;
    offload_data->phi_max = offload_data->phi_min + (n_phi - 1)*offload_data->phi_grid;

    /* Copy two phi slices into opposite ends for each field
     * component */
    memcpy(&(*offload_array)[0],
        &(*offload_array)[n_phi*phi_size],
        phi_size * sizeof(real));
    memcpy(&(*offload_array)[B_size],
        &(*offload_array)[B_size + n_phi*phi_size],
        phi_size * sizeof(real));
    memcpy(&(*offload_array)[2*B_size],
        &(*offload_array)[2*B_size + n_phi*phi_size],
        phi_size * sizeof(real));
    memcpy(&(*offload_array)[3*B_size],
        &(*offload_array)[3*B_size + n_phi*phi_size],
        phi_size * sizeof(real));

    memcpy(&(*offload_array)[phi_size],
        &(*offload_array)[(n_phi+1)*phi_size],
        phi_size * sizeof(real));
    memcpy(&(*offload_array)[B_size + phi_size],
        &(*offload_array)[B_size + (n_phi+1)*phi_size],
        phi_size * sizeof(real));
    memcpy(&(*offload_array)[2*B_size + phi_size],
        &(*offload_array)[2*B_size + (n_phi+1)*phi_size],
        phi_size * sizeof(real));
    memcpy(&(*offload_array)[3*B_size + phi_size],
        &(*offload_array)[3*B_size + (n_phi+1)*phi_size],
        phi_size * sizeof(real));

    memcpy(&(*offload_array)[(n_phi+2)*phi_size],
        &(*offload_array)[2*phi_size],
        phi_size * sizeof(real));
    memcpy(&(*offload_array)[B_size + (n_phi+2)*phi_size],
        &(*offload_array)[B_size + 2*phi_size],
        phi_size * sizeof(real));
    memcpy(&(*offload_array)[2*B_size + (n_phi+2)*phi_size],
        &(*offload_array)[2*B_size + 2*phi_size],
        phi_size * sizeof(real));
    memcpy(&(*offload_array)[3*B_size + (n_phi+2)*phi_size],
        &(*offload_array)[3*B_size + 2*phi_size],
        phi_size * sizeof(real));

    memcpy(&(*offload_array)[(n_phi+3)*phi_size],
        &(*offload_array)[3*phi_size],
        phi_size * sizeof(real));
    memcpy(&(*offload_array)[B_size + (n_phi+3)*phi_size],
        &(*offload_array)[B_size + 3*phi_size],
        phi_size * sizeof(real));
    memcpy(&(*offload_array)[2*B_size + (n_phi+3)*phi_size],
        &(*offload_array)[2*B_size + 3*phi_size],
        phi_size * sizeof(real));
    memcpy(&(*offload_array)[3*B_size + (n_phi+3)*phi_size],
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

    hdf5_close(f);
}
