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

int hdf5_bfield_read_2DS(hid_t f, B_2DS_data* data, char* qid);
int hdf5_bfield_read_3DS(hid_t f, B_3DS_data* data, char* qid);
int hdf5_bfield_read_STS(hid_t f, B_STS_data* data, char* qid);
int hdf5_bfield_read_TC(hid_t f, B_TC_data* data, char* qid);
int hdf5_bfield_read_GS(hid_t f, B_GS_data* data, char* qid);

/**
 * @brief Initialize magnetic field data from HDF5 file
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
 * @param data pointer to data struct which is initialized here
 * @param qid QID of the data that is to be read
 *
 * @return zero if reading and initialization succeeded
 */
int hdf5_bfield_init(hid_t f, B_field_data* data, char* qid) {

    char path[256]; // Storage array required for hdf5_gen_path() calls
    int err = 1;    // Error flag which is nullified if data is read succesfully

    /* Read data the QID corresponds to */

    hdf5_gen_path("/bfield/B_TC_XXXXXXXXXX", qid, path);
    if( !hdf5_find_group(f, path) ) {
        data->type = B_field_type_TC;
        err = hdf5_bfield_read_TC(f, &(data->BTC), qid);
    }

    hdf5_gen_path("/bfield/B_GS_XXXXXXXXXX", qid, path);
    if( !hdf5_find_group(f, path) ) {
        data->type = B_field_type_GS;
        err = hdf5_bfield_read_GS(f, &(data->BGS), qid);
    }

    hdf5_gen_path("/bfield/B_2DS_XXXXXXXXXX", qid, path);
    if( !hdf5_find_group(f, path) ) {
        data->type = B_field_type_2DS;
        err = hdf5_bfield_read_2DS(f, &(data->B2DS), qid);
    }

    hdf5_gen_path("/bfield/B_3DS_XXXXXXXXXX", qid, path);
    if( !hdf5_find_group(f, path) ) {
        data->type = B_field_type_3DS;
        err = hdf5_bfield_read_3DS(f, &(data->B3DS), qid);
    }

    hdf5_gen_path("/bfield/B_STS_XXXXXXXXXX", qid, path);
    if( !hdf5_find_group(f, path) ) {
        data->type = B_field_type_STS;
        err = hdf5_bfield_read_STS(f, &(data->BSTS), qid);
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
 * @param data pointer to data struct which is allocated here
 * @param qid QID of the B_2DS field that is to be read
 *
 * @return zero if reading succeeded
 */
int hdf5_bfield_read_2DS(hid_t f, B_2DS_data* data, char* qid) {
    #undef BPATH
    #define BPATH "/bfield/B_2DS_XXXXXXXXXX/"

    /* Read and initialize psi and magnetic field Rz-grid */
    int n_r, n_z;
    real r_min, r_max, z_min, z_max;
    if( hdf5_read_int(BPATH "nr", &n_r,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(BPATH "nz", &n_z,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "rmin", &r_min,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "rmax", &r_max,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "zmin", &z_min,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "zmax", &z_max,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Read psi and B values */
    real* psi = (real*)malloc(n_r*n_z*sizeof(real));
    real* br = (real*)malloc(n_r*n_z*sizeof(real));
    real* bphi = (real*)malloc(n_r*n_z*sizeof(real));
    real* bz = (real*)malloc(n_r*n_z*sizeof(real));
    if( hdf5_read_double(BPATH "psi", psi,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "br", br,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "bphi", bphi,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "bz", bz,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    real psi0, psi1, axisr, axisz;
    if( hdf5_read_double(BPATH "psi0", &psi0,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi1", &psi1,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "axisr", &axisr,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "axisz", &axisz,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    int err =  B_2DS_init(data, n_r, r_min, r_max, n_z, z_min, z_max,
                          axisr, axisz, psi0, psi1, psi, br, bphi, bz);
    free(psi);
    free(br);
    free(bphi);
    free(bz);
    return err;
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
 * @param data pointer to data struct which is allocated here
 * @param qid QID of the B_3DS field that is to be read
 *
 * @return zero if reading succeeded
 */
int hdf5_bfield_read_3DS(hid_t f, B_3DS_data* data, char* qid) {
    #undef BPATH
    #define BPATH "/bfield/B_3DS_XXXXXXXXXX/"

    /* Read and initialize magnetic field Rpz-grid */
    int b_n_r, b_n_phi, b_n_z;
    real b_r_min, b_r_max, b_phi_min, b_phi_max, b_z_min, b_z_max;
    if( hdf5_read_int(BPATH "b_nr", &b_n_r,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(BPATH "b_nz", &b_n_z,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "b_rmin", &b_r_min,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "b_rmax", &b_r_max,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "b_zmin", &b_z_min,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "b_zmax", &b_z_max,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    if( hdf5_read_int(BPATH "b_nphi", &b_n_phi,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "b_phimin", &b_phi_min,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "b_phimax", &b_phi_max,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    // Convert to radians
    b_phi_min = math_deg2rad(b_phi_min);
    b_phi_max = math_deg2rad(b_phi_max);

    /* Read and initialize psi field Rz-grid */
    int p_n_r, p_n_z;
    real p_r_min, p_r_max, p_z_min, p_z_max;
    if( hdf5_read_int(BPATH "psi_nr", &p_n_r,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(BPATH "psi_nz", &p_n_z,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi_rmin", &p_r_min,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi_rmax", &p_r_max,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi_zmin", &p_z_min,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi_zmax", &p_z_max,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Read psi and B values */
    real* psi = (real*)malloc(p_n_r*p_n_z*sizeof(real));
    real* br = (real*)malloc(b_n_r*b_n_phi*b_n_z*sizeof(real));
    real* bphi = (real*)malloc(b_n_r*b_n_phi*b_n_z*sizeof(real));
    real* bz = (real*)malloc(b_n_r*b_n_phi*b_n_z*sizeof(real));
    if( hdf5_read_double(BPATH "psi", psi,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "br", br,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "bphi", bphi,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "bz", bz,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Read the poloidal flux (psi) values at magnetic axis and separatrix. */
    real psi0, psi1, axisr, axisz;
    if( hdf5_read_double(BPATH "psi0", &psi0,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi1", &psi1,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Read magnetic axis R and z coordinates */
    if( hdf5_read_double(BPATH "axisr", &axisr,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "axisz", &axisz,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    int err = B_3DS_init(data, p_n_r, p_r_min, p_r_max, p_n_z, p_z_min, p_z_max,
                         b_n_r, b_r_min, b_r_max, b_n_phi, b_phi_min, b_phi_max,
                         b_n_z, b_z_min, b_z_max, axisr, axisz, psi0, psi1,
                         psi, br, bphi, bz);
    free(psi);
    free(br);
    free(bphi);
    free(bz);
    return err;
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
 * @param data pointer to data struct which is allocated here
 * @param qid QID of the B_STS field that is to be read
 *
 * @return zero if reading succeeded
 */
int hdf5_bfield_read_STS(hid_t f, B_STS_data* data, char* qid) {
    #undef BPATH
    #define BPATH "/bfield/B_STS_XXXXXXXXXX/"

    /* Read and initialize magnetic field Rpz-grid */
    int b_n_r, b_n_phi, b_n_z;
    real b_r_min, b_r_max, b_z_min, b_z_max, b_phi_min, b_phi_max;
    if( hdf5_read_int(BPATH "b_nr", &b_n_r,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(BPATH "b_nz", &b_n_z,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "b_rmin", &b_r_min,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "b_rmax", &b_r_max,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "b_zmin", &b_z_min,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "b_zmax", &b_z_max,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    if( hdf5_read_int(BPATH "b_nphi", &b_n_phi,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "b_phimin", &b_phi_min,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "b_phimax", &b_phi_max,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    // Convert to radians
    b_phi_min = math_deg2rad(b_phi_min);
    b_phi_max = math_deg2rad(b_phi_max);

    /* Read and initialize psi field Rpz-grid */
    int p_n_r, p_n_phi, p_n_z;
    real p_r_min, p_r_max, p_z_min, p_z_max, p_phi_min, p_phi_max;
    if( hdf5_read_int(BPATH "psi_nr", &p_n_r,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(BPATH "psi_nz", &p_n_z,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi_rmin", &p_r_min,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi_rmax", &p_r_max,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi_zmin", &p_z_min,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi_zmax", &p_z_max,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    if( hdf5_read_int(BPATH "psi_nphi", &p_n_phi,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi_phimin", &p_phi_min,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi_phimax", &p_phi_max,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    // Convert to radians
    p_phi_min = math_deg2rad(p_phi_min);
    p_phi_max = math_deg2rad(p_phi_max);

    /* Read and initialize magnetic axis phi-grid */
    int naxis;
    real axis_min, axis_max;
    if( hdf5_read_int(BPATH "axis_nphi", &naxis,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "axis_phimin", &axis_min,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "axis_phimax", &axis_max,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    // Convert to radians
    axis_min = math_deg2rad(axis_min);
    axis_max = math_deg2rad(axis_max);

    /* Read the poloidal flux (psi) values at magnetic axis and separatrix. */
    real psi0, psi1;
    if( hdf5_read_double(BPATH "psi0", &psi0,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi1", &psi1,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Read the magnetic field and psi */
    real* psi = (real*)malloc(p_n_r*p_n_phi*p_n_z*sizeof(real));
    real* br = (real*)malloc(b_n_r*b_n_phi*b_n_z*sizeof(real));
    real* bphi = (real*)malloc(b_n_r*b_n_phi*b_n_z*sizeof(real));
    real* bz = (real*)malloc(b_n_r*b_n_phi*b_n_z*sizeof(real));
    if( hdf5_read_double(BPATH "br", br,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "bphi", bphi,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "bz", bz,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi", psi,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Read the magnetic axis */
    real* axisr = (real*)malloc(naxis*sizeof(real));
    real* axisz = (real*)malloc(naxis*sizeof(real));
    if( hdf5_read_double(BPATH "axisr", axisr,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "axisz", axisz,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    int err = B_STS_init(data, p_n_r, p_r_min, p_r_max,
                         p_n_phi, p_phi_min, p_phi_max, p_n_z, p_z_min, p_z_max,
                         b_n_r, b_r_min, b_r_max, b_n_phi, b_phi_min, b_phi_max,
                         b_n_z, b_z_min, b_z_max, naxis, axis_min, axis_max,
                         axisr, axisz, psi0, psi1, psi, br, bphi, bz);
    free(psi);
    free(br);
    free(bphi);
    free(bz);
    free(axisr);
    free(axisz);

    return err;
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
 * @param data pointer to offload data struct which is allocated here
 * @param qid QID of the B_TC field that is to be read
 *
 * @return zero if reading succeeded
 */
int hdf5_bfield_read_TC(hid_t f, B_TC_data* data, char* qid) {
    #undef BPATH
    #define BPATH "/bfield/B_TC_XXXXXXXXXX/"

    real axisr, axisz, psival, rhoval, B[3], dB[9];
    if( hdf5_read_double(BPATH "axisr", &axisr,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "axisz", &axisz,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psival", &psival,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "rhoval", &rhoval,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "bxyz", B,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "jacobian", dB,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    int err = B_TC_init(data, axisr, axisz, psival, rhoval, B, dB);
    return err;
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
 * @param data pointer to data struct which is allocated here
 * @param qid QID of the B_GS field that is to be read
 *
 * @return zero if reading succeeded
 */
int hdf5_bfield_read_GS(hid_t f, B_GS_data* data, char* qid) {
    #undef BPATH
    #define BPATH "/bfield/B_GS_XXXXXXXXXX/"

    /* Equilibrium */
    real R0, z0, raxis, zaxis, B_phi0, psi0, psi1, psi_mult, psi_coeff[14];
    if( hdf5_read_double(BPATH "r0", &R0,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "z0", &z0,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "raxis", &raxis,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "zaxis", &zaxis,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "bphi0", &B_phi0,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi0", &psi0,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psi1", &psi1,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "psimult", &psi_mult,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "coefficients", psi_coeff,
                         f, qid, __FILE__, __LINE__) ) {return 1;}

    /* Ripple */
    int Nripple;
    real delta0, alpha0, a0;
    if( hdf5_read_double(BPATH "delta0", &delta0,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "alpha0", &alpha0,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(BPATH "a0", &a0,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(BPATH "nripple", &Nripple,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    int err = B_GS_init(data, R0, z0, raxis, zaxis, B_phi0, psi0, psi1,
                        psi_mult, psi_coeff, Nripple, a0, alpha0, delta0);
    return err;
}
