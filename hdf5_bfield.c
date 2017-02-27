/**
 * @file B_ST.c
 * @brief Stellarator magnetic field with cubic interpolation
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "math.h"
#include "ascot5.h"
#include "B_ST.h"
#include "B_3D.h" /* for 3D interpolation routines */
#include "B_2D.h" /* for 2D interpolation routines */
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
 * @todo Constructing full 3D magnetic field from data segment
 * @todo Reading the magnetic axis
 * @todo Reading s and volume profiles
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
void hdf5_bfield_init_offload_ST(B_ST_offload_data* offload_data, real** offload_array) {
    herr_t err;
    int version;
    int periods;
    real period_length;
    hsize_t dims[3];

    hid_t f = hdf5_open("input.h5");

    err = H5LTfind_dataset(f, "/bfield/");
    err = H5LTfind_dataset(f, "/bfield/stellarator");

    /* */
    err = H5LTget_attribute_int(f, "/bfield/stellarator", "version", &version);

    /* Number of toroidal periods */
    err = H5LTread_dataset_int(f,"/bfield/stellarator/toroidalPeriods",&periods);

    /* Read the coordinate data */
    /* Get dimensions of bfield data */
    err = H5LTget_dataset_info(f,"/bfield/stellarator/br", dims, NULL, NULL);

    /* Sizes of dimensions */
    offload_data->n_r   = dims[0];
    offload_data->n_phi = dims[1];
    offload_data->n_z   = dims[2];

    real* temp_r   = (real*) malloc(offload_data->n_r*sizeof(real));
    real* temp_phi = (real*) malloc(offload_data->n_phi*sizeof(real));
    real* temp_z   = (real*) malloc(offload_data->n_z*sizeof(real));

    err = H5LTread_dataset_double(f, "/bfield/stellarator/br", &temp_r);
    err = H5LTread_dataset_double(f, "/bfield/stellarator/bphi", &temp_phi);
    err = H5LTread_dataset_double(f, "/bfield/stellarator/bz", &temp_z);

    /* Note! phi data includes both 0 deg and 36 deg! */

    offload_data->r_min = temp_r[0];
    offload_data->r_max = temp_r[offload_data->n_r - 1];
    offload_data->r_grid = (temp_r[1] - temp_r[0]);

    /* Switch to full toroidal phi. Data only contains 1/(periods * 2) of data */
    offload_data->n_phi = (offload_data->n_phi - 1) * periods * 2 + 1;
    offload_data->phi_grid = 2*math_pi / (offload_data->n_phi);
    /* phi array starts from -0.5*phi_grid and ends at 0.5*phi_grid + 2pi! */
    offload_data->phi_min = -1.5 * offload_data->phi_grid;
    offload_data->phi_max = 1.5 * offload_data->phi_grid + 2*math_pi;

    offload_data->z_min = temp_z[0];
    offload_data->z_max = temp_z[offload_data->n_z - 1];
    offload_data->z_grid = (temp_z[1] - temp_z[0]);

    free(temp_r);
    free(temp_phi);
    free(temp_z);

    /* Read the bfield data */

    int temp_B_size = offload_data->n_r*offload_data->n_z*offload_data->n_phi;

    real* temp_B_r   = (real*) malloc(temp_B_size*sizeof(real));
    real* temp_B_phi = (real*) malloc(temp_B_size*sizeof(real));
    real* temp_B_z   = (real*) malloc(temp_B_size*sizeof(real));
    real* temp_B_s   = (real*) malloc(temp_B_size*sizeof(real));

    err = H5LTread_dataset_double(f, "/bfield/stellarator/br", &temp_B_r);
    err = H5LTread_dataset_double(f, "/bfield/stellarator/bphi", &temp_B_phi);
    err = H5LTread_dataset_double(f, "/bfield/stellarator/bz", &temp_B_z);
    err = H5LTread_dataset_double(f, "/bfield/stellarator/s", &temp_B_s);

    /* Derive the 3D magnetic field from stellarator symmetry*/

    /* We need to use the stellarator symmetry here.
     * http://dx.doi.org/10.1016/S0167-2789(97)00216-9
     * The data is expected to include half a period.
     * First figure out the effective phi angle
     */
    /* phi = modulo( rpz(2), simu%periodLength) */
    /* if( 2*phi > simu%periodLength) then */
    /*     mirrored = .true. */
    /*     phi = simu%periodLength - phi */
    /*     z = -rpz(3) */
    /* else */
    /*     mirrored = .false. */
    /*     z = rpz(3) */
    /* end if */

    *offload_array = (real*) malloc(4*B_size*sizeof(real));

    // Loop over toroidal periods
    //     Copy magnetic field components stellarator-symmetrically
    //         B_r
    //         B_phi
    //         B_z

    // Possible solution: reverse each row & then reverse each column
    
    
    
    hdf5_close(fileid);


              /////////////////////////


    /* Read phi parameters */
    fscanf(f, "%*d %*d %d %*d %*d", &(offload_data->n_phi));

    /* Read r and z parameters */
    fscanf(f, "%lf %lf %d", &(offload_data->r_min), &(offload_data->r_max),
                            &(offload_data->n_r));
    fscanf(f, "%lf %lf %d", &(offload_data->z_min), &(offload_data->z_max),
                            &(offload_data->n_z));

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

    /* Skip phimaps */
    fscanf(f, "%*f %*f");

    /* Read psi */
    int i;
    for(i = 0; i < psi_size; i++) {
        fscanf(f, "%lf", &(*offload_array)[i + 3*B_size]);
    }

    /* Calculate 2D components of poloidal field */
    real* eq_B_r = (real*) malloc(psi_size*sizeof(real));
    real* eq_B_z = (real*) malloc(psi_size*sizeof(real));

    for(i = 0; i < offload_data->n_r; i++) {
        int j;
        for(j = 0; j < offload_data->n_z; j++) {
            real psi[4];
            B_2D_bicubic_derivs(psi, 0, 0, i, j, offload_data->n_r, offload_data->r_grid, offload_data->z_grid, &(*offload_array)[3*B_size]);
            eq_B_r[j*offload_data->n_r + i] = 1/(2*math_pi) * -psi[3] / (offload_data->r_min + i*offload_data->r_grid);
            eq_B_z[j*offload_data->n_r + i] = 1/(2*math_pi) * psi[1] / (offload_data->r_min + i*offload_data->r_grid);
        }
    }

    real* temp_B_r = (real*) malloc(psi_size*offload_data->n_phi*sizeof(real));
    real* temp_B_phi = (real*) malloc(psi_size*offload_data->n_phi*sizeof(real));
    real* temp_B_z = (real*) malloc(psi_size*offload_data->n_phi*sizeof(real));

    /* Read B_r */
    for(i = 0; i < psi_size*offload_data->n_phi; i++) {
        fscanf(f, "%lf", &temp_B_r[i]);
    }

    /* Read B_phi */
    for(i = 0; i < psi_size*offload_data->n_phi; i++) {
        fscanf(f, "%lf", &temp_B_phi[i]);
    }

    /* Read B_z */
    for(i = 0; i < psi_size*offload_data->n_phi; i++) {
        fscanf(f, "%lf", &temp_B_z[i]);
    }

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
    fclose(f);

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

}
