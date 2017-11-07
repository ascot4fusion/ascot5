/**
 * @file hdf5_markers.c
 * @brief HDF5 format simulation marker input
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include "hdf5_hl.h"
#include "../particle.h"
#include "../math.h"
#include "../consts.h"
#include "hdf5_markers.h"

/**
 * @brief Write magnetic field component
 *
 * This function writes the magnetic field components in a hardcoded grid to
 * BR.out, Bphi.out and Bz.out for further processing with prepare_magn_bkg.m.
 *
 * @param Bdata pointer to magnetic field data struct
 */
void hdf5_markers_init(hid_t f, int *n, input_particle** p) {
    herr_t err;

    err = H5LTfind_dataset(f, "/markers/");
    if(err < 0) {
        return;
    }

    int n_particle;
    int n_guiding_center;
    int n_field_line;    
    err = H5LTget_attribute_int(f, "/markers/", "n_particle", &n_particle);
    err = H5LTget_attribute_int(f, "/markers/", "n_guiding_center", &n_guiding_center);
    err = H5LTget_attribute_int(f, "/markers/", "n_field_line", &n_field_line);
    
    *p = (input_particle*) malloc((n_particle + n_guiding_center + n_field_line) * sizeof(input_particle));

    /* Pointers to beginning of different data series to make code more
     * readable */
    input_particle* particle = &(*p)[0];
    input_particle* guiding_center = &(*p)[n_particle];
    input_particle* field_line = &(*p)[n_particle + n_guiding_center];

    hdf5_markers_init_particle(f, n_particle, particle);
    hdf5_markers_init_guiding_center(f, n_guiding_center, guiding_center);
    hdf5_markers_init_field_line(f, n_field_line, field_line);

    *n = n_particle + n_guiding_center + n_field_line;

}


void hdf5_markers_init_particle(hid_t f, int n, input_particle* p) {
    herr_t err;
    int i;
    
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
     
    err = H5LTread_dataset_double(f,"/markers/particle/r", r);
    err = H5LTread_dataset_double(f,"/markers/particle/phi", phi);    
    err = H5LTread_dataset_double(f,"/markers/particle/z", z);
    err = H5LTread_dataset_double(f,"/markers/particle/v_r", v_r);
    err = H5LTread_dataset_double(f,"/markers/particle/v_phi", v_phi);
    err = H5LTread_dataset_double(f,"/markers/particle/v_z", v_z);
    err = H5LTread_dataset_double(f,"/markers/particle/mass", mass);
    err = H5LTread_dataset_int(f,"/markers/particle/charge", charge);
    err = H5LTread_dataset_double(f,"/markers/particle/weight", weight);   
    err = H5LTread_dataset_double(f,"/markers/particle/time", time);
    err = H5LTread_dataset_long(f,"/markers/particle/id", id);

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


void hdf5_markers_init_guiding_center(hid_t f, int n, input_particle* p) {
    herr_t err;
    int i;
    
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
     
    err = H5LTread_dataset_double(f,"/markers/guiding_center/r", r);
    err = H5LTread_dataset_double(f,"/markers/guiding_center/phi", phi);    
    err = H5LTread_dataset_double(f,"/markers/guiding_center/z", z);
    err = H5LTread_dataset_double(f,"/markers/guiding_center/energy", energy);
    err = H5LTread_dataset_double(f,"/markers/guiding_center/pitch", pitch); 
    err = H5LTread_dataset_double(f,"/markers/guiding_center/theta", theta);
    err = H5LTread_dataset_double(f,"/markers/guiding_center/mass", mass);
    err = H5LTread_dataset_int(f,"/markers/guiding_center/charge", charge);
    err = H5LTread_dataset_double(f,"/markers/guiding_center/weight", weight);   
    err = H5LTread_dataset_double(f,"/markers/guiding_center/time", time);
    err = H5LTread_dataset_long(f,"/markers/guiding_center/id", id);
        
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

void hdf5_markers_init_field_line(hid_t f, int n, input_particle* p) {
    herr_t err;
    int i;
    
    real* r      = malloc(n * sizeof(real));
    real* phi    = malloc(n * sizeof(real));
    real* z      = malloc(n * sizeof(real));
    real* pitch  = malloc(n * sizeof(real));
    real* weight = malloc(n * sizeof(real));
    real* time   = malloc(n * sizeof(real));
    integer* id  = malloc(n * sizeof(integer));
     
    err = H5LTread_dataset_double(f,"/markers/field_line/r", r);
    err = H5LTread_dataset_double(f,"/markers/field_line/phi", phi);    
    err = H5LTread_dataset_double(f,"/markers/field_line/z", z);
    err = H5LTread_dataset_double(f,"/markers/field_line/pitch", pitch);
    err = H5LTread_dataset_double(f,"/markers/field_line/weight", weight);
    err = H5LTread_dataset_double(f,"/markers/field_line/time", time);
    err = H5LTread_dataset_long(f,"/markers/field_line/id", id);
        
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
