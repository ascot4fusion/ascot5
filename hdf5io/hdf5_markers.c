/**
 * @file hdf5_markers.c
 * @brief Read markers from HDF5 file
 *
 * Markers must be read by calling hdf5_markers_init() contained in this module
 * This module contains reading routines for all marker types.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "hdf5_helpers.h"
#include "../particle.h"
#include "../math.h"
#include "../print.h"
#include "../consts.h"
#include "hdf5_markers.h"

#define MRKPATH /**< Macro that is used to store paths to data groups */

int hdf5_markers_read_particle(hid_t f, int* nmrk, input_particle** p,
                               char* qid);
int hdf5_markers_read_guiding_center(hid_t f, int* nmrk, input_particle** p,
                                     char* qid);
int hdf5_markers_read_field_line(hid_t f, int* nmrk, input_particle** p,
                                 char* qid);

/**
 * @brief Read marker input.
 *
 * Reads all marker types and places them on input_particle array.
 *
 * @param f HDF5 file containing marker input
 * @param n pointer where the number of markers that was read will be stored
 * @param p pointer to array where markers are stored
 * @param qid QID of the marker data
 *
 * @return zero on success
 */
int hdf5_markers_read(hid_t f, int *n, input_particle** p, char* qid) {

    char path[256];
    int err = 1;

    /* Read data the QID corresponds to */

    hdf5_gen_path("/marker/prt-XXXXXXXXXX", qid, path);
    if(hdf5_find_group(f, path) == 0) {
        if( !hdf5_markers_read_particle(f, n, p, qid) ) {
            print_out(VERBOSE_IO,"\nLoaded %d particles.\n", *n);
            err = 0;
        }
    }

    hdf5_generate_qid_path("/marker/gc-XXXXXXXXXX", qid, path);
    if(hdf5_find_group(f, path) == 0) {
        if( !hdf5_markers_read_guiding_center(f, n, p, qid) ) {
            print_out(VERBOSE_IO,"\nLoaded %d guiding centers.\n", *n);
            err = 0;
        }
    }

    hdf5_generate_qid_path("/marker/fl-XXXXXXXXXX", qid, path);
    if(hdf5_find_group(f, path) == 0) {
        if( !hdf5_markers_read_field_line(f, n, p, qid) ) {
            print_out(VERBOSE_IO,"\nLoaded %d field lines.\n", *n);
            err = 0;
        }
    }

    return err;
}

/**
 * @brief Read particle input.
 *
 * Reads particles and places them on input_particle array.
 *
 * @param f HDF5 file containing marker input
 * @param nmrk pointer where the number of markers that was read will be stored
 * @param mrk pointer to array where markers are stored
 * @param qid QID of the marker data
 *
 * @return zero on success
 */
int hdf5_markers_read_particle(hid_t f, int* nmrk, input_particle** mrk,
                               char* qid) {
    #undef MRKPATH
    #define MRKPATH "/marker/prt-XXXXXXXXXX/"

    integer n;
    if( hdf5_read_long(MRKPATH "n", &n,
                       f, qid, __FILE__, __LINE__) ) {return 1;}
    *nmrk = n;

    real* r      = malloc(n * sizeof(real));
    real* phi    = malloc(n * sizeof(real));
    real* z      = malloc(n * sizeof(real));
    real* v_r    = malloc(n * sizeof(real));
    real* v_phi  = malloc(n * sizeof(real));
    real* v_z    = malloc(n * sizeof(real));
    real* mass   = malloc(n * sizeof(real));
    int* charge  = malloc(n * sizeof(int));
    real* weight = malloc(n * sizeof(real));
    real* time   = malloc(n * sizeof(real));
    integer* id  = malloc(n * sizeof(integer));

    if( hdf5_read_double(MRKPATH "r", r,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MRKPATH "phi", phi,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MRKPATH "z", z,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MRKPATH "vr", v_r,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MRKPATH "vphi", v_phi,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MRKPATH "vz", v_z,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MRKPATH "mass", mass,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(MRKPATH "charge", charge,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MRKPATH "weight", weight,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MRKPATH "time", time,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_long(MRKPATH "id", id,
                       f, qid, __FILE__, __LINE__) ) {return 1;}

    *mrk = (input_particle*) malloc(n * sizeof(input_particle));
    input_particle* p = *mrk;

    for(integer i = 0; i < n; i++) {
        p[i].p.r      = r[i];
        p[i].p.phi    = phi[i] * CONST_PI / 180;
        p[i].p.z      = z[i];
        p[i].p.v_r    = v_r[i];
        p[i].p.v_phi  = v_phi[i];
        p[i].p.v_z    = v_z[i];
        p[i].p.mass   = mass[i] * CONST_U;
        p[i].p.charge = charge[i] * CONST_E;
        p[i].p.weight = weight[i];
        p[i].p.time   = time[i];
        p[i].p.id     = (integer) id[i];
        p[i].type     = input_particle_type_p;
    }

    free(r);
    free(phi);
    free(z);
    free(v_r);
    free(v_phi);
    free(v_z);
    free(mass);
    free(charge);
    free(weight);
    free(time);
    free(id);

    return 0;
}

/**
 * @brief Read guiding center input.
 *
 * Reads guiding centers and places them on input_particle array.
 *
 * @param f HDF5 file containing marker input
 * @param nmrk pointer where the number of markers that was read will be stored
 * @param mrk pointer to array where markers are stored
 * @param qid QID of the marker data
 *
 * @return zero on success
 */
int hdf5_markers_read_guiding_center(hid_t f, int* nmrk, input_particle** mrk,
                                     char* qid) {
    #undef MRKPATH
    #define MRKPATH "/marker/gc-XXXXXXXXXX/"

    integer n;
    if( hdf5_read_long(MRKPATH "n", &n,
                       f, qid, __FILE__, __LINE__) ) {return 1;}
    *nmrk = n;

    real* r      = malloc(n * sizeof(real));
    real* phi    = malloc(n * sizeof(real));
    real* z      = malloc(n * sizeof(real));
    real* energy = malloc(n * sizeof(real));
    real* pitch  = malloc(n * sizeof(real));
    real* zeta   = malloc(n * sizeof(real));
    real* mass   = malloc(n * sizeof(real));
    int* charge  = malloc(n * sizeof(int));
    real* weight = malloc(n * sizeof(real));
    real* time   = malloc(n * sizeof(real));
    integer* id  = malloc(n * sizeof(integer));

    if( hdf5_read_double(MRKPATH "r", r,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MRKPATH "phi", phi,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MRKPATH "z", z,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MRKPATH "energy", energy,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MRKPATH "pitch", pitch,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MRKPATH "zeta", zeta,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MRKPATH "mass", mass,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(MRKPATH "charge", charge,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MRKPATH "weight", weight,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MRKPATH "time", time,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_long(MRKPATH "id", id,
                       f, qid, __FILE__, __LINE__) ) {return 1;}

    *mrk = (input_particle*) malloc(n * sizeof(input_particle));
    input_particle* p = *mrk;

    for(integer i = 0; i < n; i++) {
        p[i].p_gc.r      = r[i];
        p[i].p_gc.phi    = phi[i] * CONST_PI / 180;
        p[i].p_gc.z      = z[i];
        p[i].p_gc.energy = energy[i] * CONST_E;
        p[i].p_gc.pitch  = pitch[i];
        p[i].p_gc.zeta   = zeta[i];
        p[i].p_gc.mass   = mass[i] * CONST_U;
        p[i].p_gc.charge = charge[i] * CONST_E;
        p[i].p_gc.weight = weight[i];
        p[i].p_gc.time   = time[i];
        p[i].p_gc.id     = (integer) id[i];
        p[i].type        = input_particle_type_gc;
    }

    free(r);
    free(phi);
    free(z);
    free(energy);
    free(pitch);
    free(zeta);
    free(mass);
    free(charge);
    free(weight);
    free(time);
    free(id);

    return 0;
}

/**
 * @brief Read field line input.
 *
 * Reads field lines and places them on input_particle array.
 *
 * @param f HDF5 file containing marker input
 * @param nmrk pointer where the number of markers that was read will be stored
 * @param mrk pointer to array where markers are stored
 * @param qid QID of the marker data
 *
 * @return zero on success
 */
int hdf5_markers_read_field_line(hid_t f, int* nmrk, input_particle** mrk,
                                 char* qid) {
    #undef MRKPATH
    #define MRKPATH "/marker/fl-XXXXXXXXXX/"

    integer n;
    if( hdf5_read_long(MRKPATH "n", &n,
                       f, qid, __FILE__, __LINE__) ) {return 1;}
    *nmrk = n;

    real* r      = malloc(n * sizeof(real));
    real* phi    = malloc(n * sizeof(real));
    real* z      = malloc(n * sizeof(real));
    real* pitch  = malloc(n * sizeof(real));
    real* weight = malloc(n * sizeof(real));
    real* time   = malloc(n * sizeof(real));
    integer* id  = malloc(n * sizeof(integer));

    if( hdf5_read_double(MRKPATH "r", r,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MRKPATH "phi", phi,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MRKPATH "z", z,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MRKPATH "pitch", pitch,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MRKPATH "weight", weight,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(MRKPATH "time", time,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_long(MRKPATH "id", id,
                       f, qid, __FILE__, __LINE__) ) {return 1;}

    *mrk = (input_particle*) malloc(n * sizeof(input_particle));
    input_particle* p = *mrk;

    for(integer i = 0; i < n; i++) {
        p[i].p_ml.r      = r[i];
        p[i].p_ml.phi    = phi[i] * CONST_PI / 180;
        p[i].p_ml.z      = z[i];
        p[i].p_ml.pitch  = pitch[i];
        p[i].p_ml.weight = weight[i];
        p[i].p_ml.time   = time[i];
        p[i].p_ml.id     = (integer)id[i];
        p[i].type        = input_particle_type_ml;
    }

    free(r);
    free(phi);
    free(z);
    free(pitch);
    free(weight);
    free(time);
    free(id);

    return 0;
}
