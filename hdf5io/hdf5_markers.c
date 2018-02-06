/**
 * @file hdf5_markers.c
 * @brief HDF5 format simulation marker input
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "hdf5_helpers.h"
#include "../particle.h"
#include "../math.h"
#include "../consts.h"
#include "hdf5_markers.h"

/**
 * @brief Read marker input.
 */
int hdf5_markers_init(hid_t f, int *n, input_particle** p) {
    herr_t err;

    #if VERBOSE > 0
        printf("\nReading marker input from the HDF5 file...\n");
    #endif
    
    err = hdf5_find_group(f, "/marker/");
    if(err < 0) {
        return -1;
    }
    
    char active[11];
    err = H5LTget_attribute_string(f, "/marker/", "active", active);
    if(err < 0) {
        return -1;
    }
    active[10] = '\0';
    
    #if VERBOSE > 0
        printf("Active qid is %s\n", active);
    #endif

    /* Go through all different input types and see which one the active qid corresponds to.
     * Then read this input. */
    char path[256];
	
    hdf5_generate_qid_path("/marker/particle-XXXXXXXXXX", active, path);
    if(hdf5_find_group(f, path) == 0) {
	hdf5_markers_init_particle(f, n, p, active);
	#if VERBOSE > 0
	    printf("\nLoaded %d particles.\n", *n);
	#endif
	return 1;
    }
    
    hdf5_generate_qid_path("/marker/guiding_center-XXXXXXXXXX", active, path);
    if(hdf5_find_group(f, path) == 0) {
	hdf5_markers_init_guiding_center(f, n, p, active);
	#if VERBOSE > 0
	    printf("\nLoaded %d guiding centers.\n", *n);
	#endif
	return 1;
    }
    
    hdf5_generate_qid_path("/marker/field_line-XXXXXXXXXX", active, path);
    if(hdf5_find_group(f, path) == 0) {
	hdf5_markers_init_field_line(f, n, p, active);
	#if VERBOSE > 0
	    printf("\nLoaded %d field lines.\n", *n);
	#endif
	return 1;
    }

    return -1;
}


void hdf5_markers_init_particle(hid_t f, int* nmrk, input_particle** mrk, char* qid) {
    herr_t err;
    char path[256];
    int i;

    int n;
    err = H5LTread_dataset_long(f, hdf5_generate_qid_path("/marker/particle-XXXXXXXXXX/n", qid, path), &n);
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
     
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/marker/particle-XXXXXXXXXX/r", qid, path), r);
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/marker/particle-XXXXXXXXXX/phi", qid, path), phi);    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/marker/particle-XXXXXXXXXX/z", qid, path), z);
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/marker/particle-XXXXXXXXXX/v_r", qid, path), v_r);
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/marker/particle-XXXXXXXXXX/v_phi", qid, path), v_phi);
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/marker/particle-XXXXXXXXXX/v_z", qid, path), v_z);
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/marker/particle-XXXXXXXXXX/mass", qid, path), mass);
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/marker/particle-XXXXXXXXXX/charge", qid, path), charge);
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/marker/particle-XXXXXXXXXX/weight", qid, path), weight);   
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/marker/particle-XXXXXXXXXX/time", qid, path), time);
    err = H5LTread_dataset_long(f, hdf5_generate_qid_path("/marker/particle-XXXXXXXXXX/id", qid, path), id);

    *mrk = (input_particle*) malloc(n * sizeof(input_particle));
    input_particle* p = *mrk;
    
    for(i = 0; i < n; i++) {
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
}


void hdf5_markers_init_guiding_center(hid_t f, int* nmrk, input_particle** mrk, char* qid) {
    herr_t err;
    char path[256];
    int i;

    int n;
    err = H5LTread_dataset_long(f, hdf5_generate_qid_path("/marker/guiding_center-XXXXXXXXXX/n", qid, path), &n);
    *nmrk = n;
    
    real* r      = malloc(n * sizeof(real));
    real* phi    = malloc(n * sizeof(real));
    real* z      = malloc(n * sizeof(real));
    real* energy = malloc(n * sizeof(real));
    real* pitch  = malloc(n * sizeof(real));
    real* theta  = malloc(n * sizeof(real));
    real* mass   = malloc(n * sizeof(real));
    int* charge  = malloc(n * sizeof(int));
    real* weight = malloc(n * sizeof(real));
    real* time   = malloc(n * sizeof(real));
    integer* id  = malloc(n * sizeof(integer));
     
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/marker/guiding_center-XXXXXXXXXX/r", qid, path), r);
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/marker/guiding_center-XXXXXXXXXX/phi", qid, path), phi);    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/marker/guiding_center-XXXXXXXXXX/z", qid, path), z);
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/marker/guiding_center-XXXXXXXXXX/energy", qid, path), energy);
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/marker/guiding_center-XXXXXXXXXX/pitch", qid, path), pitch); 
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/marker/guiding_center-XXXXXXXXXX/theta", qid, path), theta);
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/marker/guiding_center-XXXXXXXXXX/mass", qid, path), mass);
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/marker/guiding_center-XXXXXXXXXX/charge", qid, path), charge);
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/marker/guiding_center-XXXXXXXXXX/weight", qid, path), weight);   
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/marker/guiding_center-XXXXXXXXXX/time", qid, path), time);
    err = H5LTread_dataset_long(f, hdf5_generate_qid_path("/marker/guiding_center-XXXXXXXXXX/id", qid, path), id);

    *mrk = (input_particle*) malloc(n * sizeof(input_particle));
    input_particle* p = *mrk;
    
    for(i = 0; i < n; i++) {
        p[i].p_gc.r      = r[i];
        p[i].p_gc.phi    = phi[i] * CONST_PI / 180;
        p[i].p_gc.z      = z[i];
        p[i].p_gc.energy = energy[i] * CONST_E;
        p[i].p_gc.pitch  = pitch[i];
	p[i].p_gc.theta  = theta[i];
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
    free(theta);
    free(mass);
    free(charge);
    free(weight);
    free(time);
    free(id);
}

void hdf5_markers_init_field_line(hid_t f, int* nmrk, input_particle** mrk, char* qid) {
    herr_t err;
    char path[256];
    int i;

    int n;
    err = H5LTread_dataset_long(f, hdf5_generate_qid_path("/marker/field_line-XXXXXXXXXX/n", qid, path), &n);
    *nmrk = n;
    
    real* r      = malloc(n * sizeof(real));
    real* phi    = malloc(n * sizeof(real));
    real* z      = malloc(n * sizeof(real));
    real* pitch  = malloc(n * sizeof(real));
    real* weight = malloc(n * sizeof(real));
    real* time   = malloc(n * sizeof(real));
    integer* id  = malloc(n * sizeof(integer));
    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/marker/field_line-XXXXXXXXXX/r", qid, path), r);
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/marker/field_line-XXXXXXXXXX/phi", qid, path), phi);    
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/marker/field_line-XXXXXXXXXX/z", qid, path), z);
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/marker/field_line-XXXXXXXXXX/pitch", qid, path), pitch);
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/marker/field_line-XXXXXXXXXX/weight", qid, path), weight);
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/marker/field_line-XXXXXXXXXX/time", qid, path), time);
    err = H5LTread_dataset_long(f, hdf5_generate_qid_path("/marker/field_line-XXXXXXXXXX/id", qid, path), id);

    *mrk = (input_particle*) malloc(n * sizeof(input_particle));
    input_particle* p = *mrk;
    
    for(i = 0; i < n; i++) {
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
}
