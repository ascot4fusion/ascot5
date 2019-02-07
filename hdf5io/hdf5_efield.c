/**
 * @file hdf5_efield.c
 * @brief Module for reading electric field data from HDF5 file
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
#include <hdf5_hl.h>
#include "../ascot5.h"
#include "../E_field.h"
#include "../Efield/E_TC.h"
#include "../Efield/E_1DS.h"
#include "../Efield/E_3D.h"
#include "../Efield/E_3DS.h"
#include "../print.h"
#include "hdf5_helpers.h"
#include "hdf5_efield.h"

int hdf5_efield_read_1DS(hid_t f, E_1DS_offload_data* offload_data,
                         real** offload_array, char* qid);
int hdf5_efield_read_TC(hid_t f, E_TC_offload_data* offload_data,
                        real** offload_array, char* qid);
int hdf5_efield_read_3D(hid_t f, E_3D_offload_data* offload_data,
                        real** offload_array, char* qid);
int hdf5_efield_read_3DS(hid_t f, E_3DS_offload_data* offload_data,
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

    hdf5_gen_path("/efield/E_3D-XXXXXXXXXX", qid, path);
    if(hdf5_find_group(f, path) == 0) {
      offload_data->type = E_field_type_3D;
      err = hdf5_efield_read_3D(f, &(offload_data->E3D),
                                offload_array, qid);
    }

    hdf5_gen_path("/efield/E_3DS-XXXXXXXXXX", qid, path);
    if(hdf5_find_group(f, path) == 0) {
      offload_data->type = E_field_type_3DS;
      err = hdf5_efield_read_3DS(f, &(offload_data->E3DS),
                                 offload_array, qid);
    }

    /* Initialize if data was read succesfully */
    if(!err) {
        err = E_field_init_offload(offload_data, offload_array);
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
}


/**
 * @brief Read electric field data of type E_3D
 *
 * The E_3D data is stored in HDF5 file under the group
 * /efield/E_3D-XXXXXXXXXX/ where X's mark the QID.
 *
 * This function assumes the group holds the following datasets:
 * (E data refers to \f$E_R\f$, \f$E_phi\f$, and \f$E_z\f$
 *
 * - int n_R Number of R grid points in the E data grid
 * - int n_phi Number of phi grid points in the E data grid
 * - int n_z Number of z grid points in the E data grid
 * - double R_min Minimum value in E data R grid [m]
 * - double R_max Maximum value in E data R grid [m]
 * - double phi_min Minimum value in E data phi grid [deg]
 * - double phi_max Maximum value in E data phi grid [deg]
 * - double z_min Minimum value in E data z grid [m]
 * - double z_max Maximum value in E data z grid [m]
 *
 * - double E_R   Electric field R component on the Rz-grid as
 *                a {nphi, nz, nR} matrix [T]
 * - double E_phi Electric field phi component on the Rz-grid as
 *                a {nphi, nz, nR} matrix [T]
 * - double E_z   Electric field z component on the Rz-grid as
 *                a {nphi, nz, nR} matrix [T]
 *
 * @param f HDF5 file identifier for a file which is opened and closed outside
 *          of this function
 * @param offload_data pointer to offload data struct which is allocated here
 * @param offload_array pointer to offload array which is allocated here and
 *                      used to store E_R, E_phi, and E_z values as
 *                      required by E_3D_init_offload()
 * @param qid QID of the E_3D field that is to be read
 *
 * @return zero if reading succeeded
 */
int hdf5_efield_read_3D(hid_t f, E_3D_offload_data* offload_data,
                         real** offload_array, char* qid) {
    #undef EPATH
    #define EPATH "/efield/E_3D-XXXXXXXXXX/"

  /* Read and initialize magnetic field Rpz-grid */
  if( hdf5_read_int(EPATH "n_R", &(offload_data->n_r),
                    f, qid, __FILE__, __LINE__) ) {return 1;}
  if( hdf5_read_int(EPATH "n_z", &(offload_data->n_z),
                    f, qid, __FILE__, __LINE__) ) {return 1;}
  if( hdf5_read_double(EPATH "R_min", &(offload_data->r_min),
                       f, qid, __FILE__, __LINE__) ) {return 1;}
  if( hdf5_read_double(EPATH "R_max", &(offload_data->r_max),
                       f, qid, __FILE__, __LINE__) ) {return 1;}
  if( hdf5_read_double(EPATH "z_min", &(offload_data->z_min),
                       f, qid, __FILE__, __LINE__) ) {return 1;}
  if( hdf5_read_double(EPATH "z_max", &(offload_data->z_max),
                       f, qid, __FILE__, __LINE__) ) {return 1;}

  if( hdf5_read_int(EPATH "n_phi", &(offload_data->n_phi),
                    f, qid, __FILE__, __LINE__) ) {return 1;}
  if( hdf5_read_double(EPATH "phi_min", &(offload_data->phi_min),
                       f, qid, __FILE__, __LINE__) ) {return 1;}
  if( hdf5_read_double(EPATH "phi_max", &(offload_data->phi_max),
                       f, qid, __FILE__, __LINE__) ) {return 1;}

  // Convert to radians
  offload_data->phi_min = math_deg2rad(offload_data->phi_min);
  offload_data->phi_max = math_deg2rad(offload_data->phi_max);

  /* Allocate offload_array storing the three components of E */
    int E_size = offload_data->n_r * offload_data->n_z
      * offload_data->n_phi;

    *offload_array = (real*) malloc((3 * E_size) * sizeof(real));
    offload_data->offload_array_length =  3 * E_size;

    /* Read the magnetic field */
    if( hdf5_read_double(EPATH "E_R", &(*offload_array)[0*E_size],
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(EPATH "E_phi", &(*offload_array)[1*E_size],
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(EPATH "E_z", &(*offload_array)[2*E_size],
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    return 0;
}

/**
 * @brief Read electric field data of type E_3DS
 *
 * The E_3DS data is stored in HDF5 file under the group
 * /efield/E_3DS-XXXXXXXXXX/ where X's mark the QID.
 *
 * This function assumes the group holds the following datasets:
 * (E data refers to \f$E_R\f$, \f$E_phi\f$, and \f$E_z\f$
 *
 * - int n_R Number of R grid points in the E data grid
 * - int n_phi Number of phi grid points in the E data grid
 * - int n_z Number of z grid points in the E data grid
 * - double R_min Minimum value in E data R grid [m]
 * - double R_max Maximum value in E data R grid [m]
 * - double phi_min Minimum value in E data phi grid [deg]
 * - double phi_max Maximum value in E data phi grid [deg]
 * - double z_min Minimum value in E data z grid [m]
 * - double z_max Maximum value in E data z grid [m]
 *
 * - double E_R   Electric field R component on the Rz-grid as
 *                a {nphi, nz, nR} matrix [T]
 * - double E_phi Electric field phi component on the Rz-grid as
 *                a {nphi, nz, nR} matrix [T]
 * - double E_z   Electric field z component on the Rz-grid as
 *                a {nphi, nz, nR} matrix [T]
 *
 * - double E_R   Electric field R component on the Rz-grid as
 *                a {nphi, nz, nR} matrix [T]
 * - double E_phi Electric field phi component on the Rz-grid as
 *                a {nphi, nz, nR} matrix [T]
 * - double E_z   Electric field z component on the Rz-grid as
 *                a {nphi, nz, nR} matrix [T]
 *
 * @param f HDF5 file identifier for a file which is opened and closed outside
 *          of this function
 * @param offload_data pointer to offload data struct which is allocated here
 * @param offload_array pointer to offload array which is allocated here and
 *                      used to store E_R, E_phi, and E_z values as
 *                      required by E_3DS_init_offload()
 * @param qid QID of the E_3DS field that is to be read
 *
 * @return zero if reading succeeded
 */
int hdf5_efield_read_3DS(hid_t f, E_3DS_offload_data* offload_data,
                         real** offload_array, char* qid) {
    #undef EPATH
    #define EPATH "/efield/E_3DS-XXXXXXXXXX/"

  /* Read and initialize magnetic field Rpz-grid */
  if( hdf5_read_int(EPATH "n_R", &(offload_data->n_r),
                    f, qid, __FILE__, __LINE__) ) {return 1;}
  if( hdf5_read_int(EPATH "n_z", &(offload_data->n_z),
                    f, qid, __FILE__, __LINE__) ) {return 1;}
  if( hdf5_read_double(EPATH "R_min", &(offload_data->r_min),
                       f, qid, __FILE__, __LINE__) ) {return 1;}
  if( hdf5_read_double(EPATH "R_max", &(offload_data->r_max),
                       f, qid, __FILE__, __LINE__) ) {return 1;}
  if( hdf5_read_double(EPATH "z_min", &(offload_data->z_min),
                       f, qid, __FILE__, __LINE__) ) {return 1;}
  if( hdf5_read_double(EPATH "z_max", &(offload_data->z_max),
                       f, qid, __FILE__, __LINE__) ) {return 1;}

  if( hdf5_read_int(EPATH "n_phi", &(offload_data->n_phi),
                    f, qid, __FILE__, __LINE__) ) {return 1;}
  if( hdf5_read_double(EPATH "phi_min", &(offload_data->phi_min),
                       f, qid, __FILE__, __LINE__) ) {return 1;}
  if( hdf5_read_double(EPATH "phi_max", &(offload_data->phi_max),
                       f, qid, __FILE__, __LINE__) ) {return 1;}


  // Convert to radians
  offload_data->phi_min = math_deg2rad(offload_data->phi_min);
  offload_data->phi_max = math_deg2rad(offload_data->phi_max);

  /* Allocate offload_array storing the three components of E */
    int E_size = offload_data->n_r * offload_data->n_z
      * offload_data->n_phi;

    *offload_array = (real*) malloc((3 * E_size) * sizeof(real));
    offload_data->offload_array_length =  3 * E_size;

    /* Read the magnetic field */
    if( hdf5_read_double(EPATH "E_R", &(*offload_array)[0*E_size],
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(EPATH "E_phi", &(*offload_array)[1*E_size],
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(EPATH "E_z", &(*offload_array)[2*E_size],
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    return 0;
}
