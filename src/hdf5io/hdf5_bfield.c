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
#include "../B_field.h"
#include "../Bfield/B_2DS.h"
#include "../Bfield/B_3DS.h"
#include "../Bfield/B_STS.h"
#include "../Bfield/B_TC.h"
#include "../Bfield/B_GS.h"
#include "hdf5_helpers.h"
#include "hdf5_bfield.h"

#define BPATH /**< Macro that is used to store paths to data groups */

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

    hdf5_gen_path("/bfield/B_TC_XXXXXXXXXX", qid, path);
    if( !hdf5_find_group(f, path) ) {
        offload_data->type = B_field_type_TC;
        err = hdf5_bfield_read_TC(f, &(offload_data->BTC),
                                  offload_array, qid);
    }

    hdf5_gen_path("/bfield/B_GS_XXXXXXXXXX", qid, path);
    if( !hdf5_find_group(f, path) ) {
        offload_data->type = B_field_type_GS;
        err = hdf5_bfield_read_GS(f, &(offload_data->BGS),
                                  offload_array, qid);
    }

    hdf5_gen_path("/bfield/B_2DS_XXXXXXXXXX", qid, path);
    if( !hdf5_find_group(f, path) ) {
        offload_data->type = B_field_type_2DS;
        err = hdf5_bfield_read_2DS(f, &(offload_data->B2DS),
                                   offload_array, qid);
    }

    hdf5_gen_path("/bfield/B_3DS_XXXXXXXXXX", qid, path);
    if( !hdf5_find_group(f, path) ) {
        offload_data->type = B_field_type_3DS;
        err = hdf5_bfield_read_3DS(f, &(offload_data->B3DS),
                                   offload_array, qid);
    }

    hdf5_gen_path("/bfield/B_STS_XXXXXXXXXX", qid, path);
    if( !hdf5_find_group(f, path) ) {
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
 * /bfield/B_2DS_XXXXXXXXXX/ where X's mark the QID.
 *
 * This function assumes the group holds the following datasets:
 *
 * - int nr Number of R grid points in the magnetic field grid
 * - int nz Number of z grid points in the magnetic field grid
 * - double rmin Minimum value in R grid [m]
 * - double rmax Maximum value in R grid [m]
 * - double zmin Minimum value in z grid [m]
 * - double zmax Maximum value in z grid [m]
 * - double axisr Magnetic axis R coordinate [m]
 * - double axisz Magnetic axis z coordinate [m]
 * - double psi0 Poloidal magnetic flux value on magnetic axis [V*s*m^-1]
 * - double psi1 Poloidal magnetic flux value on separatrix [V*s*m^-1]
 * - double psi Poloidal magnetic flux on the Rz-grid as
 *              a {nz, nR} matrix [V*s*m^-1]
 * - double br   Magnetic field R component on the Rz-grid as
 *               a {nz, nR} matrix [T]
 * - double bphi Magnetic field R component on the Rz-grid as
 *               a {nz, nR} matrix [T]
 * - double bz   Magnetic field R component on the Rz-grid as
 *               a {nz, nR} matrix [T]
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
    #define BPATH "/bfield/B_2DS_XXXXXXXXXX/"

    /* Read and initialize psi and magnetic field Rz-grid */
    if( hdf5_read_int(BPATH "nr", &(offload_data->n_r),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(BPATH "nz", &(offload_data->n_z),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "rmin", &(offload_data->r_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "rmax", &(offload_data->r_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "zmin", &(offload_data->z_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "zmax", &(offload_data->z_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Allocate offload_array; psi and each component (B_R, B_phi, B_z) has
     * size = n_r*n_z */
    int B_size = offload_data->n_r * offload_data->n_z;
    *offload_array = (real*) malloc(4 * B_size * sizeof(real));

    /* Read psi and B values */
    if( hdf5_read_double(BPATH "psi", &(*offload_array)[0*B_size],
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "br", &(*offload_array)[1*B_size],
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "bphi", &(*offload_array)[2*B_size],
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "bz", &(*offload_array)[3*B_size],
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Read the poloidal flux (psi) values at magnetic axis and separatrix. */
    if( hdf5_read_double(BPATH "psi0", &(offload_data->psi0),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi1", &(offload_data->psi1),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Read magnetic axis R and z coordinates */
    if( hdf5_read_double(BPATH "axisr", &(offload_data->axis_r),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "axisz", &(offload_data->axis_z),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    return 0;
}

/**
 * @brief Read magnetic field data of type B_3DS
 *
 * The B_3DS data is stored in HDF5 file under the group
 * /bfield/B_3DS-XXXXXXXXXX/ where X's mark the QID.
 *
 * This function assumes the group holds the following datasets:
 * (B data refers to \f$B_R\f$, \f$B_phi\f$, and \f$B_z\f$ and psi data to
 *  \f$\psi\f$.)
 *
 * - int nr Number of R grid points in the B data grid
 * - int nphi Number of phi grid points in the B data grid
 * - int nz Number of z grid points in the B data grid
 * - double b_rmin Minimum value in B data R grid [m]
 * - double b_rmax Maximum value in B data R grid [m]
 * - double b_phimin Minimum value in B data phi grid [deg]
 * - double b_phimax Maximum value in B data phi grid [deg]
 * - double b_zmin Minimum value in B data z grid [m]
 * - double b_zmax Maximum value in B data z grid [m]
 *
 * - int b_nr Number of R grid points in the psi data grid
 * - int psi_nz Number of z grid points in the psi data grid
 * - double ps_rmin Minimum value in psi data R grid [m]
 * - double psi_rmax Maximum value in psi data R grid [m]
 * - double psi_zmin Minimum value in psi data z grid [m]
 * - double psi_zmax Maximum value in psi data z grid [m]
 *
 * - double axisr Magnetic axis R coordinate [m]
 * - double axisz Magnetic axis z coordinate [m]
 * - double psi0 Poloidal magnetic flux value on magnetic axis [V*s*m^-1]
 * - double psi1 Poloidal magnetic flux value on separatrix [V*s*m^-1]
 * - double psi Poloidal magnetic flux on the Rz-grid as
 *              a {psi_nz, psi_nr} matrix [V*s*m^-1]
 * - double br   Magnetic field R component on the Rz-grid as
 *               a {b_nz, b_nphi, b_nr} matrix [T]
 * - double bphi Magnetic field R component on the Rz-grid as
 *               a {b_nz, b_nphi, b_nr} matrix [T]
 * - double bz   Magnetic field R component on the Rz-grid as
 *               a {b_nz, b_nphi, b_nr} matrix [T]
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
    #define BPATH "/bfield/B_3DS_XXXXXXXXXX/"

    /* Read and initialize magnetic field Rpz-grid */
    if( hdf5_read_int(BPATH "b_nr", &(offload_data->Bgrid_n_r),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(BPATH "b_nz", &(offload_data->Bgrid_n_z),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "b_rmin", &(offload_data->Bgrid_r_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "b_rmax", &(offload_data->Bgrid_r_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "b_zmin", &(offload_data->Bgrid_z_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "b_zmax", &(offload_data->Bgrid_z_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    if( hdf5_read_int(BPATH "b_nphi", &(offload_data->Bgrid_n_phi),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "b_phimin", &(offload_data->Bgrid_phi_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "b_phimax", &(offload_data->Bgrid_phi_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    // Convert to radians
    offload_data->Bgrid_phi_min = math_deg2rad(offload_data->Bgrid_phi_min);
    offload_data->Bgrid_phi_max = math_deg2rad(offload_data->Bgrid_phi_max);

    /* Read and initialize psi field Rz-grid */
    if( hdf5_read_int(BPATH "psi_nr", &(offload_data->psigrid_n_r),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(BPATH "psi_nz", &(offload_data->psigrid_n_z),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi_rmin", &(offload_data->psigrid_r_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi_rmax", &(offload_data->psigrid_r_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi_zmin", &(offload_data->psigrid_z_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi_zmax", &(offload_data->psigrid_z_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Allocate offload_array storing psi and the three components of B */
    int psi_size = offload_data->psigrid_n_r*offload_data->psigrid_n_z;
    int B_size = offload_data->Bgrid_n_r * offload_data->Bgrid_n_z
        * offload_data->Bgrid_n_phi;

    *offload_array = (real*) malloc((psi_size + 3 * B_size) * sizeof(real));
    offload_data->offload_array_length = psi_size + 3 * B_size;

    /* Read psi */
    if( hdf5_read_double(BPATH "psi", &(*offload_array)[3*B_size],
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Read the magnetic field */
    if( hdf5_read_double(BPATH "br", &(*offload_array)[0*B_size],
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "bphi", &(*offload_array)[1*B_size],
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "bz", &(*offload_array)[2*B_size],
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Read the poloidal flux (psi) values at magnetic axis and separatrix. */
    if( hdf5_read_double(BPATH "psi0", &(offload_data->psi0),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi1", &(offload_data->psi1),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Read magnetic axis R and z coordinates */
    if( hdf5_read_double(BPATH "axisr", &(offload_data->axis_r),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "axisz", &(offload_data->axis_z),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    return 0;
}

/**
 * @brief Read magnetic field data of type B_STS
 *
 * The B_STS data is stored in HDF5 file under the group
 * /bfield/B_STS_XXXXXXXXXX/ where X's mark the QID.
 *
 * This function assumes the group holds the following datasets:
 * (B data refers to \f$B_R\f$, \f$B_phi\f$, and \f$B_z\f$ and psi data to
 *  \f$\psi\f$.)
 *
 * - int nr Number of R grid points in the B data grid
 * - int nphi Number of phi grid points in the B data grid
 * - int nz Number of z grid points in the B data grid
 * - double b_rmin Minimum value in B data R grid [m]
 * - double b_rmax Maximum value in B data R grid [m]
 * - double b_phimin Minimum value in B data phi grid [deg]
 * - double b_phimax Maximum value in B data phi grid [deg]
 * - double b_zmin Minimum value in B data z grid [m]
 * - double b_zmax Maximum value in B data z grid [m]
 *
 * - int psi_nr Number of R grid points in the psi data grid
 * - int psi_nphi Number of phi grid points in the psi data grid
 * - int psi_nz Number of z grid points in the psi data grid
 * - double psi_rmin Minimum value in psi data R grid [m]
 * - double psi_rmax Maximum value in psi data R grid [m]
 * - double psi_phimin Minimum value in psi data phi grid [deg]
 * - double psi_phimax Maximum value in psi data phi grid [deg]
 * - double psi_zmin Minimum value in psi data z grid [m]
 * - double psi_zmax Maximum value in psi data z grid [m]
 *
 * - int n_axis Number of phi grid points in the magnetic axis data grid
 * - double axis_min Minimum value in magnetic axis data phi grid [deg]
 * - double axis_max Maximum value in magnetic axis data phi grid [deg]
 *
 * - double psi0 Poloidal magnetic flux value on magnetic axis [V*s*m^-1]
 * - double psi1 Poloidal magnetic flux value on separatrix [V*s*m^-1]
 * - double psi  Poloidal magnetic flux on the Rpz-grid as
 *               a {psi_nphi, psi_nz, psi_nr} matrix [V*s*m^-1]
 * - double br   Magnetic field R component on the Rpz-grid as
 *               a {b_nphi, b_nz, b_nr} matrix [T]
 * - double bphi Magnetic field R component on the Rpz-grid as
 *               a {b_nphi, b_nz, b_nr} matrix [T]
 * - double bz   Magnetic field R component on the Rpz-grid as
 *               a {b_nphi, b_nz, b_nr} matrix [T]
 *
 * - double axis_R  Magnetic axis R location as a {n_axis} vector [m]
 * - double axis_z  Magnetic axis R location as a {n_axis} vector [m]
 *
 * - int toroidalPeriods  Fraction of the device the data represents.
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
    #define BPATH "/bfield/B_STS_XXXXXXXXXX/"

    /* Read and initialize magnetic field Rpz-grid */
    if( hdf5_read_int(BPATH "b_nr", &(offload_data->Bgrid_n_r),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(BPATH "b_nz", &(offload_data->Bgrid_n_z),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "b_rmin", &(offload_data->Bgrid_r_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "b_rmax", &(offload_data->Bgrid_r_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "b_zmin", &(offload_data->Bgrid_z_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "b_zmax", &(offload_data->Bgrid_z_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    if( hdf5_read_int(BPATH "b_nphi", &(offload_data->Bgrid_n_phi),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "b_phimin", &(offload_data->Bgrid_phi_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "b_phimax", &(offload_data->Bgrid_phi_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    // Convert to radians
    offload_data->Bgrid_phi_min = math_deg2rad(offload_data->Bgrid_phi_min);
    offload_data->Bgrid_phi_max = math_deg2rad(offload_data->Bgrid_phi_max);

    /* Read and initialize psi field Rpz-grid */
    if( hdf5_read_int(BPATH "psi_nr", &(offload_data->psigrid_n_r),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(BPATH "psi_nz", &(offload_data->psigrid_n_z),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi_rmin", &(offload_data->psigrid_r_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi_rmax", &(offload_data->psigrid_r_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi_zmin", &(offload_data->psigrid_z_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi_zmax", &(offload_data->psigrid_z_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    if( hdf5_read_int(BPATH "psi_nphi", &(offload_data->psigrid_n_phi),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi_phimin",
                         &(offload_data->psigrid_phi_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi_phimax",
                         &(offload_data->psigrid_phi_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    // Convert to radians
    offload_data->psigrid_phi_min = math_deg2rad(offload_data->psigrid_phi_min);
    offload_data->psigrid_phi_max = math_deg2rad(offload_data->psigrid_phi_max);

    /* Read and initialize magnetic axis phi-grid */
    if( hdf5_read_int(BPATH "axis_nphi", &(offload_data->n_axis),
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "axis_phimin", &(offload_data->axis_min),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "axis_phimax", &(offload_data->axis_max),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    // Convert to radians
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
    if( hdf5_read_double(BPATH "br", &(*offload_array)[0*B_size],
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "bphi", &(*offload_array)[1*B_size],
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "bz", &(*offload_array)[2*B_size],
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Read psi */
    if( hdf5_read_double(BPATH "psi", &(*offload_array)[3*B_size],
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Read the magnetic axis */
    if( hdf5_read_double(BPATH "axisr",
                         &(*offload_array)[3*B_size + psi_size],
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "axisz",
                         &(*offload_array)[3*B_size + psi_size + axis_size],
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Read the poloidal flux (psi) values at magnetic axis and separatrix. */
    if( hdf5_read_double(BPATH "psi0", &(offload_data->psi0),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi1", &(offload_data->psi1),
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    return 0;
}

/**
 * @brief Read magnetic field data of type B_TC
 *
 * The B_TC data is stored in HDF5 file under the group
 * /bfield/B_TC_XXXXXXXXXX/ where X's mark the QID.
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
 * - double jacobian Magnetic field Jacobian [dB_x/dx, dB_x/dy, dB_x/dz, dB_y/dx,
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
    #define BPATH "/bfield/B_TC_XXXXXXXXXX/"

    if( hdf5_read_double(BPATH "axisr", &(offload_data->axisr),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "axisz", &(offload_data->axisz),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psival", &(offload_data->psival),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "rhoval", &(offload_data->rhoval),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "bxyz", offload_data->B,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "jacobian", offload_data->dB,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    *offload_array = NULL;
    return 0;
}

/**
 * @brief Read magnetic field data of type B_GS
 *
 * The B_GS data is stored in HDF5 file under the group
 * /bfield/B_GS_XXXXXXXXXX/ where X's mark the QID.
 *
 * This function assumes the group holds the following datasets:
 *
 * - double r0 Magnetic axis R coordinate [m]
 * - double z0 Magnetic axis z coordinate [m]
 * - double bphi0 Toroidal magnetic field value at magnetic axis [T]
 * - double psi0 Poloidal flux value at magnetic axis [V*s*m^-1]
 * - double psi1 Poloidal flux value at separatrix [V*s*m^-1]
 * - double psi_mult Scaling factor for psi
 * - double coefficients Coefficients for evaluating psi [c_1, c_2, ..., c_12, A]
 * - double delta0 Ripple strength
 * - double alpha0 Ripple penetration
 * - double a0 Minor radius [m]
 * - int nripple Number of toroidal field coils
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
    #define BPATH "/bfield/B_GS_XXXXXXXXXX/"

    /* Equilibrium */
    if( hdf5_read_double(BPATH "r0", &(offload_data->R0),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "z0", &(offload_data->z0),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "raxis", &(offload_data->raxis),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "zaxis", &(offload_data->zaxis),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "bphi0", &(offload_data->B_phi0),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi0", &(offload_data->psi0),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi1", &(offload_data->psi1),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psimult", &(offload_data->psi_mult),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "coefficients", offload_data->psi_coeff,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Ripple */
    if( hdf5_read_double(BPATH "delta0", &(offload_data->delta0),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "alpha0", &(offload_data->alpha0),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "a0", &(offload_data->a0),
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(BPATH "nripple", &(offload_data->Nripple),
                      f, qid, __FILE__, __LINE__) ) {return 1;}

    *offload_array = NULL;

    return 0;
}
