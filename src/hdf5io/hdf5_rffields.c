/**
 * @file hdf5_rffields.c
 * @brief Module for reading RF wave electric and magnetic fields from the HDF5.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "../ascot5.h"
#include "../print.h"
#include "hdf5_helpers.h"
#include "hdf5_rffields.h"
#include "../rf_fields_fo.h"

int hdf5_efield_init(hid_t f, RF2D_fields* data, char* qid) {
    char path[256];
    int err = 1;

    #ifdef RFPATH
    #undef RFPATH
    #endif
    #define RFPATH "/rffields/RF_2D_XXXXXXXXXX"

    // Internal parameters to store temporarily what 
    // we read from the HDF5 file
    int n_r, n_z; // Number of grid points in r and z directions
    real rmin, rmax, zmin, zmax; // Min and max values of r and z
    real omega; // Frequency of the RF field
    int n_tor; // Toroidal mode number

    /* Read data the QID corresponds to */
    hdf5_gen_path(RFPATH, qid, path);
    if(hdf5_find_group(f, path) == 0) {
        print_err("Error: RF_2D group not found in HDF5 file\n");
        return 1;
    }

    /* Read data */
    if( hdf5_read_int(RFPATH "n_r", &n_r,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(RFPATH "n_z", &n_z,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(RFPATH "rmin", &rmin,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(RFPATH "rmax", &rmax,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(RFPATH "zmin", &zmin,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(RFPATH "zmax", &zmax,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    
    int nsize = n_r*n_z;
    real* buffer = (real*) malloc(12*nsize*sizeof(real));
    real* fields[12]; // Pointers to the fields in the buffer
    const char** fieldnames = {"Er_real", "Er_imag", "Ephi_real", "Ephi_imag", "Ez_real", "Ez_imag",
                               "Br_real", "Br_imag", "Bphi_real", "Bphi_imag", "Bz_real", "Bz_imag"
                             };
    for(int i = 0; i < 12; i++) {
        fields[i] = buffer + i*nsize;
        if( hdf5_read_double(RFPATH "field", fieldnames[i],
                             f, qid, __FILE__, __LINE__) ) {err = 1; break;}
    }
    if(err) {
        free(buffer);
        return 1;
    }

    /* Initialize the RF2D_fields struct */
    err = RF2D_fields_init(data, rmin, rmax, n_r, zmin, zmax, n_z, n_tor, omega, \
                           fields[0], fields[1], fields[2], fields[3], fields[4], fields[5], \
                           fields[6], fields[7], fields[8], fields[9], fields[10], fields[11]);
    free(buffer);

    return err;

}