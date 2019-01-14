/**
 * @file hdf5_bfield.c
 * @brief Module for reading magnetic field data from HDF5 file
 *
 * This module contains routines to read magnetic field data from ASCOT5 HDF5
 * file. The reading is done by calling hdf5_bfield_init_offload() which reads
 * the data and uses that to initialize the offload data.
 *
 * Note: the routines within this module that read the data from HDF5 file
 * may use the offload data struct and offload arrays as (temporary) storage.
 * However, the actual initialization is done at the specific
 * B_field_init_offload() function the magnetic field data corresponds to. Check
 * from that function what the offload data and the offload array are expected
 * to contain. As a rule of thumb, the reading routines here should only read
 * the data and maybe do some trivial computations but nothing complicated.
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "../ascot5.h"
#include "../math.h"
#include "../print.h"
#include "../consts.h"
#include "../symmetry.h"
#include "../B_field.h"
#include "../Bfield/B_2DS.h"
#include "../Bfield/B_3DS.h"
#include "../Bfield/B_3DS_T.h"
#include "../Bfield/B_STS.h"
<<<<<<< HEAD
##include "hdf5_helpers.h"
=======
#include "../Bfield/B_TC.h"
#include "../Bfield/B_GS.h"
#include "hdf5_helpers.h"
>>>>>>> develop
#include "hdf5_bfield.h"

int hdf5_bfield_read_2DS(hid_t f, B_2DS_offload_data* offload_data,
                         real** offload_array, char* qid);
int hdf5_bfield_read_3DS(hid_t f, B_3DS_offload_data* offload_data,
                         real** offload_array, char* qid);
int hdf5_bfield_read_STS(hid_t f, B_STS_offload_data* offload_data,
                         real** offload_array, char* qid);
int hdf5_bfield_read_TC(hid_t f, B_TC_offload_data* offload_data,
                        real** offload_array, char* qid);
int hdf5_bfield_read_GS(hid_t f, B_GS_offload_data* offload_data,
                        real** offload_array, char* qid);

/**
 * @brief Initialize magnetic field offload data from HDF5 file
 *
 * This function determines to which type of magnetic field the given QID
 * corresponds to, and calls the corresponding reading routine. The read data is
 * then used to initialize the offload data by calling B_field_init_offload().
 *
 * The magnetic field data is stored under /bfield/ group in ASCOT5 HDF5 file.
 * Several magnetic fields of same or different type maybe stored there in
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
int hdf5_bfield_init_offload(hid_t f, B_field_offload_data* offload_data,
                             real** offload_array, char* qid) {

    char path[256]; // Storage array required for hdf5_gen_path() calls
    int err = 1;    // Error flag which is nullified if data is read succesfully

    /* Read data the QID corresponds to */

    hdf5_gen_path("/bfield/B_TC-XXXXXXXXXX", qid, path);
    if( !hdf5_find_group(f, path) ) {
        offload_data->type = B_field_type_TC;
        err = hdf5_bfield_read_TC(f, &(offload_data->BTC),
                                  offload_array, qid);
    }

    hdf5_gen_path("/bfield/B_GS-XXXXXXXXXX", qid, path);
    if( !hdf5_find_group(f, path) ) {
        offload_data->type = B_field_type_GS;
        err = hdf5_bfield_read_GS(f, &(offload_data->BGS),
                                  offload_array, qid);
    }

    hdf5_gen_path("/bfield/B_2DS-XXXXXXXXXX", qid, path);
    if( !hdf5_find_group(f, path) ) {
        offload_data->type = B_field_type_2DS;
        err = hdf5_bfield_read_2DS(f, &(offload_data->B2DS),
                                   offload_array, qid);
    }

    hdf5_gen_path("/bfield/B_3DS-XXXXXXXXXX", qid, path);
    if( !hdf5_find_group(f, path) ) {
        offload_data->type = B_field_type_3DS;
        err = hdf5_bfield_read_3DS(f, &(offload_data->B3DS),
                                   offload_array, qid);
    }

<<<<<<< HEAD
    hdf5_generate_qid_path("/bfield/B_3DS_T-XXXXXXXXXX", active, path);
    if(hdf5_find_group(f, path) == 0) {
      hdf5_bfield_init_offload_3DS_T(f, &(offload_data->B3DST), offload_array, active);
      offload_data->type = B_field_type_3DS_T;
      offload_data->offload_array_length=offload_data->B3DST.offload_array_length;

        #if VERBOSE > 0
      printf("\nLoaded 3D magnetic field (B_3DS_T) time interpolated\n");
      printf("with parameters:\n");
      printf("- %d time slices\n",
	     offload_data->B3DST.n_time);
      printf("- t = [");
      int i = 0;
      for (i=0;i<(offload_data->B3DST.n_time-1);i++){
	printf("%le, ",offload_data->B3DST.time[i]);
      }
      printf("%le]\n",offload_data->B3DST.time[offload_data->B3DST.n_time-1]);
      printf("- rmin, rmax, nr = %le, %le, %d\n",
	     offload_data->B3DST.psigrid_r_min,offload_data->B3DST.psigrid_r_max,offload_data->B3DST.psigrid_n_r);
      printf("- zmin, zmax, nz = %le, %le, %d\n",
	     offload_data->B3DST.psigrid_z_min,offload_data->B3DST.psigrid_z_max,offload_data->B3DST.psigrid_n_z);
      printf("with parameters:\n");
        #endif
      return 1;
    }

    hdf5_generate_qid_path("/bfield/B_STS-XXXXXXXXXX", active, path);
    if(hdf5_find_group(f, path) == 0) {
        hdf5_bfield_init_offload_STS(f, &(offload_data->BSTS), offload_array, active);
=======
    hdf5_gen_path("/bfield/B_STS-XXXXXXXXXX", qid, path);
    if( !hdf5_find_group(f, path) ) {
>>>>>>> develop
        offload_data->type = B_field_type_STS;
        err = hdf5_bfield_read_STS(f, &(offload_data->BSTS),
                                   offload_array, qid);
    }

    /* Initialize if data was read succesfully */
    if(!err) {
        err = B_field_init_offload(offload_data, offload_array);
    }

    return err;
}

/**
 * @brief Read magnetic field data of type B_2DS
 *
 * The B_2DS data is stored in HDF5 file under the group
 * /bfield/B_2DS-XXXXXXXXXX/ where X's mark the QID.
 *
 * This function assumes the group holds the following datasets:
 *
 * - int n_R Number of R grid points in the magnetic field grid
 * - int n_z Number of z grid points in the magnetic field grid
 * - double R_min Minimum value in R grid [m]
 * - double R_max Maximum value in R grid [m]
 * - double z_min Minimum value in z grid [m]
 * - double z_max Maximum value in z grid [m]
 * - double axis_R Magnetic axis R coordinate [m]
 * - double axis_z Magnetic axis z coordinate [m]
 * - double psi0 Poloidal magnetic flux value on magnetic axis [V*s*m^-1]
 * - double psi1 Poloidal magnetic flux value on separatrix [V*s*m^-1]
 * - double psi Poloidal magnetic flux on the Rz-grid as
 *              a {nz, nR} matrix [V*s*m^-1]
 * - double B_R   Magnetic field R component on the Rz-grid as
 *                a {nz, nR} matrix [T]
 * - double B_phi Magnetic field R component on the Rz-grid as
 *                a {nz, nR} matrix [T]
 * - double B_z   Magnetic field R component on the Rz-grid as
 *                a {nz, nR} matrix [T]
 *
 * @param f HDF5 file identifier for a file which is opened and closed outside
 *          of this function
 * @param offload_data pointer to offload data struct which is allocated here
 * @param offload_array pointer to offload array which is allocated here and
 *                      used to store psi, B_R, B_phi, and B_z values as
 *                      required by B_2DS_init_offload()
 * @param qid QID of the B_2DS field that is to be read
 *
 * @return zero if reading succeeded
 */
int hdf5_bfield_read_2DS(hid_t f, B_2DS_offload_data* offload_data,
                         real** offload_array, char* qid) {
    #undef BPATH
    #define BPATH "/bfield/B_2DS-XXXXXXXXXX/"

    /* Read and initialize psi and magnetic field Rz-grid */
    if( hdf5_read_int(BPATH "n_R", &(offload_data->n_r),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(BPATH "n_z", &(offload_data->n_z),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "R_min", &(offload_data->r_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "R_max", &(offload_data->r_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "z_min", &(offload_data->z_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "z_max", &(offload_data->z_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Allocate offload_array; psi and each component (B_R, B_phi, B_z) has
     * size = n_r*n_z */
    int B_size = offload_data->n_r * offload_data->n_z;
    *offload_array = (real*) malloc(4 * B_size * sizeof(real));

    /* Read psi and B values */
    if( hdf5_read_double(BPATH "psi", &(*offload_array)[0*B_size],
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "B_R", &(*offload_array)[1*B_size],
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "B_phi", &(*offload_array)[2*B_size],
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "B_z", &(*offload_array)[3*B_size],
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Read the poloidal flux (psi) values at magnetic axis and separatrix. */
    if( hdf5_read_double(BPATH "psi0", &(offload_data->psi0),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi1", &(offload_data->psi1),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Read magnetic axis R and z coordinates */
    if( hdf5_read_double(BPATH "axis_R", &(offload_data->axis_r),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "axis_z", &(offload_data->axis_z),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    return 0;
}

/**
 * @brief Read magnetic field data of type B_2DS
 *
 * The B_3DS data is stored in HDF5 file under the group
 * /bfield/B_3DS-XXXXXXXXXX/ where X's mark the QID.
 *
 * This function assumes the group holds the following datasets:
 * (B data refers to \f$B_R\f$, \f$B_phi\f$, and \f$B_z\f$ and psi data to
 *  \f$\psi\f$.)
 *
 * - int n_R Number of R grid points in the B data grid
 * - int n_phi Number of phi grid points in the B data grid
 * - int n_z Number of z grid points in the B data grid
 * - double R_min Minimum value in B data R grid [m]
 * - double R_max Maximum value in B data R grid [m]
 * - double phi_min Minimum value in B data phi grid [deg]
 * - double phi_max Maximum value in B data phi grid [deg]
 * - double z_min Minimum value in B data z grid [m]
 * - double z_max Maximum value in B data z grid [m]
 *
 * - int n_R Number of R grid points in the psi data grid
 * - int n_phi Number of phi grid points in the psi data grid
 * - int n_z Number of z grid points in the psi data grid
 * - double R_min Minimum value in psi data R grid [m]
 * - double R_max Maximum value in psi data R grid [m]
 * - double z_min Minimum value in psi data z grid [m]
 * - double z_max Maximum value in psi data z grid [m]
 *
 * - double axis_R Magnetic axis R coordinate [m]
 * - double axis_z Magnetic axis z coordinate [m]
 * - double psi0 Poloidal magnetic flux value on magnetic axis [V*s*m^-1]
 * - double psi1 Poloidal magnetic flux value on separatrix [V*s*m^-1]
 * - double psi Poloidal magnetic flux on the Rz-grid as
 *              a {nz, nR} matrix [V*s*m^-1]
 * - double B_R   Magnetic field R component on the Rz-grid as
 *                a {nphi, nz, nR} matrix [T]
 * - double B_phi Magnetic field R component on the Rz-grid as
 *                a {nphi, nz, nR} matrix [T]
 * - double B_z   Magnetic field R component on the Rz-grid as
 *                a {nphi, nz, nR} matrix [T]
 *
 * @param f HDF5 file identifier for a file which is opened and closed outside
 *          of this function
 * @param offload_data pointer to offload data struct which is allocated here
 * @param offload_array pointer to offload array which is allocated here and
 *                      used to store psi, B_R, B_phi, and B_z values as
 *                      required by B_3DS_init_offload()
 * @param qid QID of the B_3DS field that is to be read
 *
 * @return zero if reading succeeded
 */
int hdf5_bfield_read_3DS(hid_t f, B_3DS_offload_data* offload_data,
                         real** offload_array, char* qid) {
    #undef BPATH
    #define BPATH "/bfield/B_3DS-XXXXXXXXXX/"

    /* Read and initialize magnetic field Rpz-grid */
    if( hdf5_read_int(BPATH "n_R", &(offload_data->Bgrid_n_r),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(BPATH "n_z", &(offload_data->Bgrid_n_z),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "R_min", &(offload_data->Bgrid_r_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "R_max", &(offload_data->Bgrid_r_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "z_min", &(offload_data->Bgrid_z_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "z_max", &(offload_data->Bgrid_z_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    if( hdf5_read_int(BPATH "n_phi", &(offload_data->n_phi),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "phi_min", &(offload_data->phi_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "phi_max", &(offload_data->phi_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    // Convert to radians
    offload_data->phi_min = math_deg2rad(offload_data->phi_min);
    offload_data->phi_max = math_deg2rad(offload_data->phi_max);

    /* Read and initialize psi field Rz-grid */
    if( hdf5_read_int(BPATH "psigrid_n_R", &(offload_data->psigrid_n_r),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(BPATH "psigrid_n_z", &(offload_data->psigrid_n_z),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psigrid_R_min", &(offload_data->psigrid_r_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psigrid_R_max", &(offload_data->psigrid_r_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psigrid_z_min", &(offload_data->psigrid_z_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psigrid_z_max", &(offload_data->psigrid_z_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Allocate offload_array storing psi and the three components of B */
    int psi_size = offload_data->psigrid_n_r*offload_data->psigrid_n_z;
    int B_size = offload_data->Bgrid_n_r * offload_data->Bgrid_n_z
        * offload_data->n_phi;

    *offload_array = (real*) malloc((psi_size + 3 * B_size) * sizeof(real));
    offload_data->offload_array_length = psi_size + 3 * B_size;

    /* Read psi */
    if( hdf5_read_double(BPATH "psi", &(*offload_array)[3*B_size],
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Read the magnetic field */
    if( hdf5_read_double(BPATH "B_R", &(*offload_array)[0*B_size],
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "B_phi", &(*offload_array)[1*B_size],
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "B_z", &(*offload_array)[2*B_size],
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Read the poloidal flux (psi) values at magnetic axis and separatrix. */
    if( hdf5_read_double(BPATH "psi0", &(offload_data->psi0),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi1", &(offload_data->psi1),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Read magnetic axis R and z coordinates */
    if( hdf5_read_double(BPATH "axis_R", &(offload_data->axis_r),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "axis_z", &(offload_data->axis_z),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    return 0;
}





void hdf5_bfield_init_offload_3DS_T(hid_t f, B_3DS_T_offload_data* offload_data, real** offload_array, char* qid) {
  herr_t err;
  char path[256];

  
  /* Read number of time slices */

  err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/bfield/B_3DS_T-XXXXXXXXXX/n_time", qid, path), &(offload_data->n_time));
  if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;} 
  if(offload_data->n_time>N_MAX_TIME_SLICE){
    err = 46;
    printf("Error: n_time > N_MAX_TIME_SLICE");
    return;
  }

  err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_3DS_T-XXXXXXXXXX/time", qid, path), &(offload_data->time));
  if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}


  /* Read and initialize magnetic field Rz-grid */
  err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/bfield/B_3DS_T-XXXXXXXXXX/n_R", qid, path), &(offload_data->Bgrid_n_r));
  if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

  err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/bfield/B_3DS_T-XXXXXXXXXX/n_phi", qid, path), &(offload_data->n_phi));
  if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

  err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/bfield/B_3DS_T-XXXXXXXXXX/n_z", qid, path), &(offload_data->Bgrid_n_z));
  if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

  err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_3DS_T-XXXXXXXXXX/R_min", qid, path), &(offload_data->Bgrid_r_min));
  if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

  err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_3DS_T-XXXXXXXXXX/R_max", qid, path), &(offload_data->Bgrid_r_max));
  if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

  err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_3DS_T-XXXXXXXXXX/phi_min", qid, path), &(offload_data->phi_min));
  if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

  err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_3DS_T-XXXXXXXXXX/phi_max", qid, path), &(offload_data->phi_max));
  if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;

  err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_3DS_T-XXXXXXXXXX/z_min", qid, path), &(offload_data->Bgrid_z_min));
  if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

  err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_3DS_T-XXXXXXXXXX/z_max", qid, path), &(offload_data->Bgrid_z_max));
  if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}


  // Convert to radians
  offload_data->phi_min = offload_data->phi_min*(180.0/math_pi);
  offload_data->phi_max = offload_data->phi_max*(180.0/math_pi);

  /* Calculate grid size */
  offload_data->Bgrid_r_grid = (offload_data->Bgrid_r_max - offload_data->Bgrid_r_min)
    / (offload_data->Bgrid_n_r - 1);
  offload_data->Bgrid_z_grid = (offload_data->Bgrid_z_max - offload_data->Bgrid_z_min)
    / (offload_data->Bgrid_n_z - 1);
  offload_data->phi_grid = (offload_data->phi_max - offload_data->phi_min)
    / (offload_data->n_phi - 1);

  /* Read and initialize psi field Rz-grid */
  err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/bfield/B_3DS_T-XXXXXXXXXX/psigrid_n_R", qid, path), &(offload_data->psigrid_n_r));
  if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

  err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/bfield/B_3DS_T-XXXXXXXXXX/psigrid_n_z", qid, path), &(offload_data->psigrid_n_z));
  if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

  err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_3DS_T-XXXXXXXXXX/psigrid_R_min", qid, path), &(offload_data->psigrid_r_min));
  if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

  err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_3DS_T-XXXXXXXXXX/psigrid_R_max", qid, path), &(offload_data->psigrid_r_max));
  if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

  err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_3DS_T-XXXXXXXXXX/psigrid_z_min", qid, path), &(offload_data->psigrid_z_min));
  if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

  err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_3DS_T-XXXXXXXXXX/psigrid_z_max", qid, path), &(offload_data->psigrid_z_max));
  if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

  offload_data->psigrid_r_grid = (offload_data->psigrid_r_max - offload_data->psigrid_r_min)
    / (offload_data->psigrid_n_r - 1);
  offload_data->psigrid_z_grid = (offload_data->psigrid_z_max - offload_data->psigrid_z_min)
    / (offload_data->psigrid_n_z - 1);

  /* Allocate offload_array */
  int psi_size = offload_data->psigrid_n_r*offload_data->psigrid_n_z*offload_data->n_time;
  int B_size = offload_data->Bgrid_n_r*offload_data->Bgrid_n_z*offload_data->n_phi*offload_data->n_time;

  *offload_array = (real*) malloc((psi_size + 3 * B_size) * sizeof(real));
  offload_data->offload_array_length = psi_size + 3 * B_size;

  /* Read psi */
  err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_3DS_T-XXXXXXXXXX/psi", qid, path), &(*offload_array)[3*B_size]);
  if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

  /* Read the magnetic field */
  err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_3DS_T-XXXXXXXXXX/B_R", qid, path), &(*offload_array)[0*B_size]);
  if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

  err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_3DS_T-XXXXXXXXXX/B_phi", qid, path), &(*offload_array)[1*B_size]);
  if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

  err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_3DS_T-XXXXXXXXXX/B_z", qid, path), &(*offload_array)[2*B_size]);
  if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

  /* Read the first two values; These are the poloidal flux (psi) values at
   * magnetic axis and at x point (that is, separatrix). */
  err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_3DS_T-XXXXXXXXXX/psi0", qid, path), &(offload_data->psi0));
  if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

  err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_3DS_T-XXXXXXXXXX/psi1", qid, path), &(offload_data->psi1));
  if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

  /* Read magnetic axis r and z coordinates */
  err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_3DS_T-XXXXXXXXXX/axis_R", qid, path), &(offload_data->axis_r));
  if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}

  err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/bfield/B_3DS_T-XXXXXXXXXX/axis_z", qid, path), &(offload_data->axis_z));
  if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return;}
}





/**
 * @brief Read magnetic field data of type B_STS
 *
 * The B_STS data is stored in HDF5 file under the group
 * /bfield/B_STS-XXXXXXXXXX/ where X's mark the QID.
 *
 * This function assumes the group holds the following datasets:
 * (B data refers to \f$B_R\f$, \f$B_phi\f$, and \f$B_z\f$ and psi data to
 *  \f$\psi\f$.)
 *
 * - int n_R Number of R grid points in the B data grid
 * - int n_phi Number of phi grid points in the B data grid
 * - int n_z Number of z grid points in the B data grid
 * - double R_min Minimum value in B data R grid [m]
 * - double R_max Maximum value in B data R grid [m]
 * - double phi_min Minimum value in B data phi grid [deg]
 * - double phi_max Maximum value in B data phi grid [deg]
 * - double z_min Minimum value in B data z grid [m]
 * - double z_max Maximum value in B data z grid [m]
 *
 * - int psigrid_n_R Number of R grid points in the psi data grid
 * - int psigrid_n_phi Number of phi grid points in the psi data grid
 * - int psigrid_n_z Number of z grid points in the psi data grid
 * - double psigrid_R_min Minimum value in psi data R grid [m]
 * - double psigrid_R_max Maximum value in psi data R grid [m]
 * - double psigrid_phi_min Minimum value in psi data phi grid [deg]
 * - double psigrid_phi_max Maximum value in psi data phi grid [deg]
 * - double psigrid_z_min Minimum value in psi data z grid [m]
 * - double psigrid_z_max Maximum value in psi data z grid [m]
 *
 * - int n_axis Number of phi grid points in the magnetic axis data grid
 * - double axis_min Minimum value in magnetic axis data phi grid [deg]
 * - double axis_max Maximum value in magnetic axis data phi grid [deg]
 *
 * - double psi0 Poloidal magnetic flux value on magnetic axis [V*s*m^-1]
 * - double psi1 Poloidal magnetic flux value on separatrix [V*s*m^-1]
 * - double psi   Poloidal magnetic flux on the Rpz-grid as
 *                a {nphi, nz, nR} matrix [V*s*m^-1]
 * - double B_R   Magnetic field R component on the Rpz-grid as
 *                a {nphi, nz, nR} matrix [T]
 * - double B_phi Magnetic field R component on the Rpz-grid as
 *                a {nphi, nz, nR} matrix [T]
 * - double B_z   Magnetic field R component on the Rpz-grid as
 *                a {nphi, nz, nR} matrix [T]
 *
 * - double axis_R  Magnetic axis R location as a {n_axis} vector [m]
 * - double axis_z  Magnetic axis R location as a {n_axis} vector [m]
 *
 * - int toroidalPeriods  Fraction of the device the data represents.
 * - int symmetry_mode    Symmetry type of the data.
 *                (0 = stellarator symmetric, 1 = periodic)
 *
 * @param f HDF5 file identifier for a file which is opened and closed outside
 *          of this function
 * @param offload_data pointer to offload data struct which is allocated here
 * @param offload_array pointer to offload array which is allocated here and
 *                      used to store psi, B_R, B_phi, B_z, axis_r, and axis_z
 *                      values as required by B_STS_init_offload()
 * @param qid QID of the B_STS field that is to be read
 *
 * @return zero if reading succeeded
 */
int hdf5_bfield_read_STS(hid_t f, B_STS_offload_data* offload_data,
                                 real** offload_array, char* qid) {
    #undef BPATH
    #define BPATH "/bfield/B_STS-XXXXXXXXXX/"

    /* Read and initialize magnetic field Rpz-grid */
    if( hdf5_read_int(BPATH "n_R", &(offload_data->Bgrid_n_r),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(BPATH "n_z", &(offload_data->Bgrid_n_z),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "R_min", &(offload_data->Bgrid_r_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "R_max", &(offload_data->Bgrid_r_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "z_min", &(offload_data->Bgrid_z_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "z_max", &(offload_data->Bgrid_z_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    if( hdf5_read_int(BPATH "n_phi", &(offload_data->Bgrid_n_phi),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "phi_min", &(offload_data->Bgrid_phi_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "phi_max", &(offload_data->Bgrid_phi_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    // Convert to radians
    offload_data->Bgrid_phi_min = math_deg2rad(offload_data->Bgrid_phi_min);
    offload_data->Bgrid_phi_max = math_deg2rad(offload_data->Bgrid_phi_max);

    /* Read and initialize psi field Rpz-grid */
    if( hdf5_read_int(BPATH "psigrid_n_R", &(offload_data->psigrid_n_r),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(BPATH "psigrid_n_z", &(offload_data->psigrid_n_z),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psigrid_R_min", &(offload_data->psigrid_r_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psigrid_R_max", &(offload_data->psigrid_r_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psigrid_z_min", &(offload_data->psigrid_z_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psigrid_z_max", &(offload_data->psigrid_z_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    if( hdf5_read_int(BPATH "psigrid_n_phi", &(offload_data->psigrid_n_phi),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psigrid_phi_min",
                         &(offload_data->psigrid_phi_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psigrid_phi_max",
                         &(offload_data->psigrid_phi_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    // Convert to radians
    offload_data->psigrid_phi_min = math_deg2rad(offload_data->psigrid_phi_min);
    offload_data->psigrid_phi_max = math_deg2rad(offload_data->psigrid_phi_max);

    /* Read and initialize magnetic axis phi-grid */
    if( hdf5_read_int(BPATH "n_axis", &(offload_data->n_axis),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "axis_min", &(offload_data->axis_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "axis_max", &(offload_data->axis_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Note! phi data for magnetic axis in deg */
    offload_data->axis_min = math_deg2rad(offload_data->axis_min);
    offload_data->axis_max = math_deg2rad(offload_data->axis_max);

    /* Allocate offload_array storing psi and the three components of B */
    int psi_size = offload_data->psigrid_n_r*offload_data->psigrid_n_z
        * offload_data->psigrid_n_phi;
    int B_size = offload_data->Bgrid_n_r * offload_data->Bgrid_n_z
        * offload_data->Bgrid_n_phi;
    int axis_size = offload_data->n_axis;

    *offload_array = (real*) malloc((psi_size + 3 * B_size + 2 * axis_size)
                                    * sizeof(real));
    offload_data->offload_array_length = psi_size + 3 * B_size + 2 * axis_size;

    /* Read the magnetic field */
    if( hdf5_read_double(BPATH "B_R", &(*offload_array)[0*B_size],
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "B_phi", &(*offload_array)[1*B_size],
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "B_z", &(*offload_array)[2*B_size],
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Read psi */
    if( hdf5_read_double(BPATH "psi", &(*offload_array)[3*B_size],
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Read the magnetic axis */
    if( hdf5_read_double(BPATH "axis_R",
                         &(*offload_array)[3*B_size + psi_size],
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "axis_z",
                         &(*offload_array)[3*B_size + psi_size + axis_size],
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Read the poloidal flux (psi) values at magnetic axis and separatrix. */
    if( hdf5_read_double(BPATH "psi0", &(offload_data->psi0),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi1", &(offload_data->psi1),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Number of toroidal periods */
    if( hdf5_read_int(BPATH "toroidalPeriods", &(offload_data->n_periods),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    /* Symmetry mode */
    int symmetry_mode;
    if( hdf5_read_int(BPATH "symmetry_mode", &symmetry_mode,
                         f, qid, __FILE__, __LINE__) ) {
        /* If error, default to stellarator symmetry */
        symmetry_mode = 0;
    }
    switch(symmetry_mode) {
        case 0:
            offload_data->symmetry_mode = symmetry_type_stellarator;
            break;
        case 1:
            offload_data->symmetry_mode = symmetry_type_periodic;
            break;
        default:
            return 1;
            break;
    }

    return 0;
}

/**
 * @brief Read magnetic field data of type B_TC
 *
 * The B_TC data is stored in HDF5 file under the group
 * /bfield/B_TC-XXXXXXXXXX/ where X's mark the QID.
 *
 * This function assumes the group holds the following datasets:
 *
 * - double axisr A value that is returned when magnetic axis R coordinate is
 *                requested [m]
 * - double axisz A value that is returned when magnetic axis z coordinate is
 *                requested [m]
 * - double psival A value that is returned when poloidal magnetic flux value is
 *                 requested [V*s*m^-1]
 * - double rhoval A value that is returned when normalized poloidal flux value
 *                 is requested
 * - double Bxyz Magnetic field values at origo: [B_x, B_y, B_z] [T]
 * - double J Magnetic field Jacobian [dB_x/dx, dB_x/dy, dB_x/dz, dB_y/dx,
 *           dB_y/dy, dB_y/dz, dB_z/dx, dB_z/dy, dB_z/dz] [T/m]
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
int hdf5_bfield_read_TC(hid_t f, B_TC_offload_data* offload_data,
                        real** offload_array, char* qid) {
    #undef BPATH
    #define BPATH "/bfield/B_TC-XXXXXXXXXX/"

    if( hdf5_read_double(BPATH "axisr", &(offload_data->axisr),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "axisz", &(offload_data->axisz),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psival", &(offload_data->psival),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "rhoval", &(offload_data->rhoval),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "Bxyz", offload_data->B,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "J", offload_data->dB,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    *offload_array = NULL;
    return 0;
}

/**
 * @brief Read magnetic field data of type B_GS
 *
 * The B_GS data is stored in HDF5 file under the group
 * /bfield/B_GS-XXXXXXXXXX/ where X's mark the QID.
 *
 * This function assumes the group holds the following datasets:
 *
 * - double R0 Magnetic axis R coordinate [m]
 * - double z0 Magnetic axis z coordinate [m]
 * - double B_phi0 Toroidal magnetic field value at magnetic axis [T]
 * - double psi0 Poloidal flux value at magnetic axis [V*s*m^-1]
 * - double psi1 Poloidal flux value at separatrix [V*s*m^-1]
 * - double psi_mult Scaling factor for psi
 * - double psi_coeff Coefficients for evaluating psi [c_1, c_2, ..., c_12, A]
 * - double delta0 Ripple strength
 * - double alpha0 Ripple penetration
 * - double a0 Minor radius [m]
 * - int Nripple Number of toroidal field coils
 *
 * @param f HDF5 file identifier for a file which is opened and closed outside
 *          of this function
 * @param offload_data pointer to offload data struct which is allocated here
 * @param offload_array pointer to offload array but no data is stored there so
 *                      it is not allocated and NULL pointer is returned instead
 * @param qid QID of the B_GS field that is to be read
 *
 * @return zero if reading succeeded
 */
int hdf5_bfield_read_GS(hid_t f, B_GS_offload_data* offload_data,
                        real** offload_array, char* qid) {
    #undef BPATH
    #define BPATH "/bfield/B_GS-XXXXXXXXXX/"

    /* Equilibrium */
    if( hdf5_read_double(BPATH "R0", &(offload_data->R0),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "z0", &(offload_data->z0),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "B_phi0", &(offload_data->B_phi0),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi0", &(offload_data->psi0),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi1", &(offload_data->psi1),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi_mult", &(offload_data->psi_mult),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi_coeff", offload_data->psi_coeff,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Ripple */
    if( hdf5_read_double(BPATH "delta0", &(offload_data->delta0),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "alpha0", &(offload_data->alpha0),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "a0", &(offload_data->a0),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(BPATH "Nripple", &(offload_data->Nripple),
                      f, qid, __FILE__, __LINE__) ) {return 1;}

    *offload_array = NULL;

    return 0;
}
