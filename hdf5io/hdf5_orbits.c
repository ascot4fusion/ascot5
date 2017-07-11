#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "hdf5_helpers.h"
#include "../simulate.h"
#include "../particle.h"
#include "../ascot5.h"
#include "../diag_orb.h"

void hdf5_orbits_write(sim_data* sim) {
    diag_orb_data* diag = &sim->diag_data.orbits;
    int size = diag->size;
    if(!size) {
	// No data here!
	return;
    }
    hid_t file = hdf5_open("ascot.h5");
    hid_t group = hdf5_create_group(file, "orbits");
    hsize_t dims[1];
    dims[0] = size;

    if(diag->type == diag_orb_type_fo) {
	group = hdf5_create_group(group, "fo");
	if(diag->writelist) {
	    
	    int i;
	    diag_orb_dat* top = diag->writelist;
	    real* data = malloc(size*sizeof(real));
	    diag_orb_dat* list;
	    
	    list = top;
	    for(i = 0; i < size; i++) {
		data[i] = list->fo.r;
		list = list->prev;
	    }
	    H5LTmake_dataset(group, "R", 1, dims, H5T_IEEE_F64LE, data);
	    H5LTset_attribute_string(group, "R", "unit", "deg");

	    list = top;
	    for(i = 0; i < size; i++) {
		data[i] = list->fo.phi;
		list = list->prev;
	    }
	    H5LTmake_dataset(group, "phi", 1, dims, H5T_IEEE_F64LE, data);
	    H5LTset_attribute_string(group, "phi", "unit", "deg");

	    list = top;
	    for(i = 0; i < size; i++) {
		data[i] = list->fo.z;
		list = list->prev;
	    }
	    H5LTmake_dataset(group, "z", 1, dims, H5T_IEEE_F64LE, data);
	    H5LTset_attribute_string(group, "z", "unit", "deg");

	    list = top;
	    for(i = 0; i < size; i++) {
		data[i] = list->fo.rdot;
		list = list->prev;
	    }
	    H5LTmake_dataset(group, "vR", 1, dims, H5T_IEEE_F64LE, data);
	    H5LTset_attribute_string(group, "vR", "unit", "deg");

	    list = top;
	    for(i = 0; i < size; i++) {
		data[i] = list->fo.r * list->fo.phidot;
		list = list->prev;
	    }
	    H5LTmake_dataset(group, "vphi", 1, dims, H5T_IEEE_F64LE, data);
	    H5LTset_attribute_string(group, "vpar", "unit", "deg");

	    list = top;
	    for(i = 0; i < size; i++) {
		data[i] = list->fo.zdot;
		list = list->prev;
	    }
	    H5LTmake_dataset(group, "vz", 1, dims, H5T_IEEE_F64LE, data);
	    H5LTset_attribute_string(group, "vz", "unit", "deg");

	    free(data);

	    real* intdata = malloc(size*sizeof(integer));
	    list = top;
	    for(i = 0; i < size; i++) {
		intdata[i] = list->fo.id;
		list = list->prev;
	    }
	    H5LTmake_dataset(group, "id", 1, dims, H5T_STD_I64LE, intdata);
	    H5LTset_attribute_string(group, "id", "unit", "deg");

	    free(intdata);
	}
    }
    else if(diag->type == diag_orb_type_gc) {
	group = hdf5_create_group(group, "gc");
	if(diag->writelist) {
	    
	    int i;
	    diag_orb_dat* top = diag->writelist;
	    real* data = malloc(size*sizeof(real));
	    diag_orb_dat* list;
	    
	    list = top;
	    for(i = 0; i < size; i++) {
		data[i] = list->gc.r;
		list = list->prev;
	    }
	    H5LTmake_dataset(group, "R", 1, dims, H5T_IEEE_F64LE, data);
	    H5LTset_attribute_string(group, "R", "unit", "deg");

	    list = top;
	    for(i = 0; i < size; i++) {
		data[i] = list->gc.phi;
		list = list->prev;
	    }
	    H5LTmake_dataset(group, "phi", 1, dims, H5T_IEEE_F64LE, data);
	    H5LTset_attribute_string(group, "phi", "unit", "deg");

	    list = top;
	    for(i = 0; i < size; i++) {
		data[i] = list->gc.z;
		list = list->prev;
	    }
	    H5LTmake_dataset(group, "z", 1, dims, H5T_IEEE_F64LE, data);
	    H5LTset_attribute_string(group, "z", "unit", "deg");

	    list = top;
	    for(i = 0; i < size; i++) {
		data[i] = list->gc.mu;
		list = list->prev;
	    }
	    H5LTmake_dataset(group, "mu", 1, dims, H5T_IEEE_F64LE, data);
	    H5LTset_attribute_string(group, "mu", "unit", "deg");

	    list = top;
	    for(i = 0; i < size; i++) {
		data[i] = list->gc.vpar;
		list = list->prev;
	    }
	    H5LTmake_dataset(group, "vpar", 1, dims, H5T_IEEE_F64LE, data);
	    H5LTset_attribute_string(group, "vpar", "unit", "deg");

	    free(data);

	    real* intdata = malloc(size*sizeof(integer));
	    list = top;
	    for(i = 0; i < size; i++) {
		intdata[i] = list->gc.id;
		list = list->prev;
	    }
	    H5LTmake_dataset(group, "id", 1, dims, H5T_STD_I64LE, intdata);
	    H5LTset_attribute_string(group, "id", "unit", "deg");

	    free(intdata);
	}
    }
    else if(diag->type == diag_orb_type_ml) {
	group = hdf5_create_group(group, "ml");
	
	if(diag->writelist) {
	    
	    int i;
	    diag_orb_dat* top = diag->writelist;
	    real* data = malloc(size*sizeof(real));
	    diag_orb_dat* list;
	    
	    list = top;
	    for(i = 0; i < size; i++) {
		data[i] = list->ml.r;
		list = list->prev;
	    }
	    H5LTmake_dataset(group, "R", 1, dims, H5T_IEEE_F64LE, data);
	    H5LTset_attribute_string(group, "R", "unit", "deg");

	    list = top;
	    for(i = 0; i < size; i++) {
		data[i] = list->ml.phi;
		list = list->prev;
	    }
	    H5LTmake_dataset(group, "phi", 1, dims, H5T_IEEE_F64LE, data);
	    H5LTset_attribute_string(group, "phi", "unit", "deg");

	    list = top;
	    for(i = 0; i < size; i++) {
		data[i] = list->ml.z;
		list = list->prev;
	    }
	    H5LTmake_dataset(group, "z", 1, dims, H5T_IEEE_F64LE, data);
	    H5LTset_attribute_string(group, "z", "unit", "deg");

	    free(data);

	    real* intdata = malloc(size*sizeof(integer));
	    list = top;
	    for(i = 0; i < size; i++) {
		intdata[i] = list->ml.id;
		list = list->prev;
	    }
	    H5LTmake_dataset(group, "id", 1, dims, H5T_STD_I64LE, intdata);
	    H5LTset_attribute_string(group, "id", "unit", "deg");

	    free(intdata);
	}
    }

    hdf5_close(file);
}
