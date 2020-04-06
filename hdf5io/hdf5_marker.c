/**
 * @file hdf5_markers.c
 * @brief Read markers from HDF5 file
 *
 * Markers must be read by calling hdf5_marker_init() contained in this module
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
#include "hdf5_marker.h"

#define MRKPATH /**< Macro that is used to store paths to data groups */

int hdf5_marker_read_particle(hid_t f, int* nmrk, input_particle** p,
                              char* qid);
int hdf5_marker_read_guiding_center(hid_t f, int* nmrk, input_particle** p,
                                    char* qid);
int hdf5_marker_read_field_line(hid_t f, int* nmrk, input_particle** p,
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
int hdf5_marker_read(hid_t f, int *n, input_particle** p, char* qid) {

    char path[256];
    int err = 1;

    /* Read data the QID corresponds to */

    hdf5_gen_path("/marker/prt_XXXXXXXXXX", qid, path);
    if(hdf5_find_group(f, path) == 0) {
        if( !hdf5_marker_read_particle(f, n, p, qid) ) {
            print_out(VERBOSE_IO,"\nLoaded %d particles.\n", *n);
            err = 0;
        }
    }

    hdf5_generate_qid_path("/marker/gc_XXXXXXXXXX", qid, path);
    if(hdf5_find_group(f, path) == 0) {
        if( !hdf5_marker_read_guiding_center(f, n, p, qid) ) {
            print_out(VERBOSE_IO,"\nLoaded %d guiding centers.\n", *n);
            err = 0;
        }
    }

    hdf5_generate_qid_path("/marker/fl_XXXXXXXXXX", qid, path);
    if(hdf5_find_group(f, path) == 0) {
        if( !hdf5_marker_read_field_line(f, n, p, qid) ) {
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
int hdf5_marker_read_particle(hid_t f, int* nmrk, input_particle** mrk,
                              char* qid) {
    #undef MRKPATH
    #define MRKPATH "/marker/prt_XXXXXXXXXX/"

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
    int* anum    = malloc(n * sizeof(int));
    int* znum    = malloc(n * sizeof(int));
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
    if( hdf5_read_int(MRKPATH "anum", anum,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(MRKPATH "znum", znum,
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
        p[i].p.anum   = anum[i];
        p[i].p.znum   = znum[i];
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
    free(anum);
    free(znum);
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
int hdf5_marker_read_guiding_center(hid_t f, int* nmrk, input_particle** mrk,
                                    char* qid) {
    #undef MRKPATH
    #define MRKPATH "/marker/gc_XXXXXXXXXX/"

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
    int* anum    = malloc(n * sizeof(int));
    int* znum    = malloc(n * sizeof(int));
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
    if( hdf5_read_int(MRKPATH "anum", anum,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(MRKPATH "znum", znum,
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
        p[i].p_gc.anum   = anum[i];
        p[i].p_gc.znum   = znum[i];
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
    free(anum);
    free(znum);
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
int hdf5_marker_read_field_line(hid_t f, int* nmrk, input_particle** mrk,
                                char* qid) {
    #undef MRKPATH
    #define MRKPATH "/marker/fl_XXXXXXXXXX/"

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

/**
 * @brief Write particle input.
 *
 * Write particles ín the input_particle array to hdf5 file.
 *
 * @param f HDF5 file to be written
 * @param n number of markers to be stored
 * @param p pointer to array where markers are stored
 * @param qid QID for the marker data
 *
 * @return zero on success
 */
int hdf5_marker_write_particle(hid_t f, int n, input_particle* p, char* qid) {
    #undef MRKPATH
    #define MRKPATH "/marker/prt_XXXXXXXXXX/"

    char path[256];
    hdf5_gen_path(MRKPATH, qid, path);

    hdf5_create_group(f, path);
    hid_t grp = H5Gopen(f, path, H5P_DEFAULT);

    real* r      = malloc(n * sizeof(real));
    real* phi    = malloc(n * sizeof(real));
    real* z      = malloc(n * sizeof(real));
    real* v_r    = malloc(n * sizeof(real));
    real* v_phi  = malloc(n * sizeof(real));
    real* v_z    = malloc(n * sizeof(real));
    real* mass   = malloc(n * sizeof(real));
    int* charge  = malloc(n * sizeof(real));
    int* anum    = malloc(n * sizeof(int));
    int* znum    = malloc(n * sizeof(int));
    real* weight = malloc(n * sizeof(real));
    real* time   = malloc(n * sizeof(real));
    integer* id  = malloc(n * sizeof(integer));

    for(int i = 0; i < n; i++) {
        r[i] = p[i].p.r;
        phi[i] = p[i].p.phi;
        z[i] = p[i].p.z;
        v_r[i] = p[i].p.v_r;
        v_phi[i] = p[i].p.v_phi;
        v_z[i] = p[i].p.v_z;
        mass[i] = p[i].p.mass / CONST_U;
        charge[i] = round(p[i].p.charge / CONST_E);
        anum[i] = p[i].p.anum;
        znum[i] = p[i].p.znum;
        weight[i] = p[i].p.weight;
        time[i] = p[i].p.time;
        id[i] = p[i].p.id;
    }

    hdf5_write_extendible_dataset_int(grp, "n", 1, &n);
    hdf5_write_extendible_dataset_double(grp, "r", n, r);
    hdf5_write_extendible_dataset_double(grp, "phi", n, phi);
    hdf5_write_extendible_dataset_double(grp, "z", n, z);
    hdf5_write_extendible_dataset_double(grp, "vr", n, v_r);
    hdf5_write_extendible_dataset_double(grp, "vphi", n, v_phi);
    hdf5_write_extendible_dataset_double(grp, "vz", n, v_z);
    hdf5_write_extendible_dataset_double(grp, "mass", n, mass);
    hdf5_write_extendible_dataset_int(grp, "charge", n, charge);
    hdf5_write_extendible_dataset_int(grp, "anum", n, anum);
    hdf5_write_extendible_dataset_int(grp, "znum", n, znum);
    hdf5_write_extendible_dataset_double(grp, "weight", n, weight);
    hdf5_write_extendible_dataset_double(grp, "time", n, time);
    hdf5_write_extendible_dataset_long(grp, "id", n, id);

    hdf5_write_string_attribute(f, path, "description",  "");
    hdf5_write_string_attribute(f, path, "date",  "");

    /* Set this run as active. */
    hdf5_write_string_attribute(f, "/marker", "active",  qid);

    H5Gclose(grp);

    free(r);
    free(phi);
    free(z);
    free(v_r);
    free(v_phi);
    free(v_z);
    free(mass);
    free(charge);
    free(anum);
    free(znum);
    free(weight);
    free(time);
    free(id);

    return 0;
}
