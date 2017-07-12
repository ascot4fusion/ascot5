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
#include "../consts.h"

void hdf5_orbits_writeset(hid_t group,  diag_orb_dat_type type, diag_orb_dat* list, int size, int* mask, const char* dataset);
void hdf5_orbits_writeset_fo(hid_t group, diag_orb_dat* list, int size, int* mask);
void hdf5_orbits_writeset_gc(hid_t group, diag_orb_dat* list, int size, int* mask);
void hdf5_orbits_writeset_ml(hid_t group, diag_orb_dat* list, int size, int* mask);

void hdf5_orbits_write(sim_data* sim, char* out) {
    diag_orb_data* diag = &sim->diag_data.orbits;
    diag_orb_dat* top = diag->writelist;

    int size = diag->size;
    int i;
    if(!size) {
	// No data here!
	return;
    }
    hid_t file = hdf5_open(out);
    hid_t group = hdf5_create_group(file, "orbits");
    hid_t grp;
	    
    if(diag->mode == DIAG_ORB_ORBIT) {
	
	int mask[size];
	for(i=0; i<size; i++) {
	    mask[i] = 1;
	}

	if(diag->type == diag_orb_type_fo) {
	    grp = hdf5_create_group(group, "fo"); 
	    hdf5_orbits_writeset_fo(grp, top, size, mask);
	}
	else if(diag->type == diag_orb_type_gc) {
	    grp = hdf5_create_group(group, "gc"); 
	    hdf5_orbits_writeset_gc(grp, top, size, mask);
	}
	else if(diag->type == diag_orb_type_ml) {
	    grp = hdf5_create_group(group, "ml"); 
	    hdf5_orbits_writeset_ml(grp, top, size, mask);
	}
    }
    else if(diag->mode == DIAG_ORB_POINCARE) {
	
	for(int ip=0; ip<diag->npoloidalplots; ip++) {
	    
	    diag_orb_dat* list = top;
	    int mask[size];
	    for(i=0; i<size; i++) {
		mask[i] = 0;
		if(list->poincareId == ip) {
		    mask[i] = 1;
		}
		list = list->prev;
	    }

	    char groupname[16];
	    if(diag->type == diag_orb_type_fo) {
		sprintf(groupname, "fopol%d", (int)diag->poloidalangles[ip]);
		grp = hdf5_create_group(group, groupname); 
		hdf5_orbits_writeset_fo(grp, top, size, mask);
	    }
	    else if(diag->type == diag_orb_type_gc) {
		sprintf(groupname, "gcpol%d", (int)diag->poloidalangles[ip]);
		grp = hdf5_create_group(group, groupname); 
		hdf5_orbits_writeset_gc(grp, top, size, mask);
	    }
	    else if(diag->type == diag_orb_type_ml) {
		sprintf(groupname, "mlpol%d", (int)diag->poloidalangles[ip]);
		grp = hdf5_create_group(group, groupname);
		hdf5_orbits_writeset_ml(grp, top, size, mask);
	    }
	}

	for(int ip=0; ip<diag->ntoroidalplots; ip++) {
	    
	    diag_orb_dat* list = top;
	    int mask[size];
	    for(i=0; i<size; i++) {
		mask[i] = 0;
		if((list->poincareId - DIAG_ORB_MAXPOINCARES) == ip) {
		    mask[i] = 1;
		}
		list = list->prev;
	    }

	    char groupname[16];
	    if(diag->type == diag_orb_type_fo) {
		sprintf(groupname, "fotor%d", (int)diag->poloidalangles[ip]);
	        grp = hdf5_create_group(group, groupname); 
		hdf5_orbits_writeset_fo(grp, top, size, mask);
	    }
	    else if(diag->type == diag_orb_type_gc) {
		sprintf(groupname, "gctor%d", (int)diag->poloidalangles[ip]);
		grp = hdf5_create_group(group, groupname);
		hdf5_orbits_writeset_gc(grp, top, size, mask);
	    }
	    else if(diag->type == diag_orb_type_ml) {
		sprintf(groupname, "mltor%d", (int)diag->poloidalangles[ip]);
		grp = hdf5_create_group(group, groupname);
		hdf5_orbits_writeset_ml(grp, top, size, mask);
	    }
	}
    }


    hdf5_close(file);
}

void hdf5_orbits_writeset(hid_t group,  diag_orb_dat_type type, diag_orb_dat* list, int size, int* mask, const char* dataset) {
   
    int i;      // List index
    int id = 0; // Datavec index
    integer datasize = 0;
    for(i = 0; i < size; i++) {
	datasize += mask[i];
    }
    hsize_t dims[1];
    dims[0] = datasize;
    
    if(!datasize) {
	// Empty set
	return;
    }

    /* Common fields (R,phi,z,time,id)*/
    if(strcmp(dataset,"R")==0) {
	real* datavec = (real*)malloc(datasize*sizeof(real));
	
	for(i = 0; i < size; i++) {
	    if(!mask[i]){continue;}
	    
	    if(type == diag_orb_type_fo) {
		datavec[id] = list->fo.r;
	    }
	    else if(type == diag_orb_type_gc) {
		datavec[id] = list->gc.r;
	    }
	    else if(type == diag_orb_type_ml) {
		datavec[id] = list->ml.r;
	    }
	    list = list->prev;
	    id++;
	}

	H5LTmake_dataset(group, "R", 1, dims, H5T_IEEE_F64LE, datavec);
	H5LTset_attribute_string(group, "R", "unit", "m");

	free(datavec);
	return;
    }

    if(strcmp(dataset,"phi")==0) {
	real* datavec = (real*)malloc(datasize*sizeof(real));

	for(i = 0; i < size; i++) {
	    if(!mask[i]){continue;}

	    if(type == diag_orb_type_fo) {
		datavec[id] = list->fo.phi*(CONST_PI/180);
	    }
	    else if(type == diag_orb_type_gc) {
		datavec[id] = list->gc.phi*(CONST_PI/180);
	    }
	    else if(type == diag_orb_type_ml) {
		datavec[id] = list->ml.phi*(CONST_PI/180);
	    }
	    list = list->prev;
	    id++;
	}

	H5LTmake_dataset(group, "phi", 1, dims, H5T_IEEE_F64LE, datavec);
	H5LTset_attribute_string(group, "phi", "unit", "deg");

	free(datavec);
	return;
    }

    if(strcmp(dataset,"z")==0) {
	real* datavec = (real*)malloc(datasize*sizeof(real));

	for(i = 0; i < size; i++) {
	    if(!mask[i]){continue;}

	    if(type == diag_orb_type_fo) {
		datavec[id] = list->fo.z;
	    }
	    else if(type == diag_orb_type_gc) {
		datavec[id] = list->gc.z;
	    }
	    else if(type == diag_orb_type_ml) {
		datavec[id] = list->ml.z;
	    }
	    list = list->prev;
	    id++;
	}

	H5LTmake_dataset(group, "z", 1, dims, H5T_IEEE_F64LE, datavec);
	H5LTset_attribute_string(group, "z", "unit", "m");

	free(datavec);
	return;
    }

    if(strcmp(dataset,"time")==0) {
	real* datavec = (real*)malloc(datasize*sizeof(real));

	for(i = 0; i < size; i++) {
	    if(!mask[i]){continue;}

	    if(type == diag_orb_type_fo) {
		datavec[id] = list->fo.time;
	    }
	    else if(type == diag_orb_type_gc) {
		datavec[id] = list->gc.time;
	    }
	    else if(type == diag_orb_type_ml) {
		datavec[id] = list->ml.time;
	    }
	    list = list->prev;
	    id++;
	}

	H5LTmake_dataset(group, "time", 1, dims, H5T_IEEE_F64LE, datavec);
	H5LTset_attribute_string(group, "time", "unit", "s");

	free(datavec);
	return;
    }

    if(strcmp(dataset,"id")==0) {
	integer* datavec = (integer*)malloc(datasize*sizeof(integer));

	for(i = 0; i < size; i++) {
	    if(!mask[i]){continue;}

	    if(type == diag_orb_type_fo) {
		datavec[id] = list->fo.id;
	    }
	    else if(type == diag_orb_type_gc) {
		datavec[id] = list->gc.id;
	    }
	    else if(type == diag_orb_type_ml) {
		datavec[id] = list->ml.id;
	    }
	    list = list->prev;
	    id++;
	}

	H5LTmake_dataset(group, "id", 1, dims, H5T_STD_I64LE, datavec);
	H5LTset_attribute_string(group, "id", "unit", "1");

	free(datavec);
	return;
    }

    /* Magnetic field components */
     if(strcmp(dataset,"B_r")==0) {
	real* datavec = (real*)malloc(datasize*sizeof(real));

	for(i = 0; i < size; i++) {
	    if(!mask[i]){continue;}

	    if(type == diag_orb_type_fo) {
		datavec[id] = list->fo.B_r;
	    }
	    else if(type == diag_orb_type_gc) {
		datavec[id] = list->gc.B_r;
	    }
	    else if(type == diag_orb_type_ml) {
		datavec[id] = list->ml.B_r;
	    }
	    list = list->prev;
	    id++;
	}

	H5LTmake_dataset(group, "B_R", 1, dims, H5T_IEEE_F64LE, datavec);
	H5LTset_attribute_string(group, "B_R", "unit", "T");

	free(datavec);
	return;
    }

    if(strcmp(dataset,"B_phi")==0) {
	real* datavec = (real*)malloc(datasize*sizeof(real));

	for(i = 0; i < size; i++) {
	    if(!mask[i]){continue;}

	    if(type == diag_orb_type_fo) {
		datavec[id] = list->fo.B_phi;
	    }
	    else if(type == diag_orb_type_gc) {
		datavec[id] = list->gc.B_phi;
	    }
	    else if(type == diag_orb_type_ml) {
		datavec[id] = list->ml.B_phi;
	    }
	    list = list->prev;
	    id++;
	}

	H5LTmake_dataset(group, "B_phi", 1, dims, H5T_IEEE_F64LE, datavec);
	H5LTset_attribute_string(group, "B_phi", "unit", "T");

	free(datavec);
	return;
    }

    if(strcmp(dataset,"B_z")==0) {
	real* datavec = (real*)malloc(datasize*sizeof(real));

	for(i = 0; i < size; i++) {
	    if(!mask[i]){continue;}

	    if(type == diag_orb_type_fo) {
		datavec[id] = list->fo.B_z;
	    }
	    else if(type == diag_orb_type_gc) {
		datavec[id] = list->gc.B_z;
	    }
	    else if(type == diag_orb_type_ml) {
		datavec[id] = list->ml.B_z;
	    }
	    list = list->prev;
	    id++;
	}

	H5LTmake_dataset(group, "B_z", 1, dims, H5T_IEEE_F64LE, datavec);
	H5LTset_attribute_string(group, "B_z", "unit", "T");

	free(datavec);
	return;
    }   

    /* Particle/guiding center specific (mass, charge, weight) */
    if(strcmp(dataset,"mass")==0) {
	real* datavec = (real*)malloc(datasize*sizeof(real));

	for(i = 0; i < size; i++) {
	    if(!mask[i]){continue;}

	    if(type == diag_orb_type_fo) {
		datavec[id] = list->fo.mass;
	    }
	    else if(type == diag_orb_type_gc) {
		datavec[id] = list->gc.mass;
	    }
	    list = list->prev;
	    id++;
	}

	H5LTmake_dataset(group, "mass", 1, dims, H5T_IEEE_F64LE, datavec);
	H5LTset_attribute_string(group, "mass", "unit", "kg");

	free(datavec);
	return;
    }
    
    if(strcmp(dataset,"charge")==0) {
	int* datavec = (int*)malloc(datasize*sizeof(int));

	for(i = 0; i < size; i++) {
	    if(!mask[i]){continue;}

	    if(type == diag_orb_type_fo) {
		datavec[id] = list->fo.charge;
	    }
	    else if(type == diag_orb_type_gc) {
		datavec[id] = list->gc.charge;
	    }
	    list = list->prev;
	    id++;
	}

	H5LTmake_dataset(group, "charge", 1, dims, H5T_STD_I64LE, datavec);
	H5LTset_attribute_string(group, "charge", "unit", "e");

	free(datavec);
	return;
    }

    if(strcmp(dataset,"weight")==0) {
	real* datavec = (real*)malloc(datasize*sizeof(real));

	for(i = 0; i < size; i++) {
	    if(!mask[i]){continue;}

	    if(type == diag_orb_type_fo) {
		datavec[id] = list->fo.weight;
	    }
	    else if(type == diag_orb_type_gc) {
		datavec[id] = list->gc.weight;
	    }
	    list = list->prev;
	    id++;
	}

	H5LTmake_dataset(group, "weight", 1, dims, H5T_IEEE_F64LE, datavec);
	H5LTset_attribute_string(group, "weight", "unit", "markers/s");

	free(datavec);
	return;
    }

    /* Particle specific */
    if(strcmp(dataset,"v_R")==0) {
	real* datavec = (real*)malloc(datasize*sizeof(real));

	for(i = 0; i < size; i++) {
	    if(!mask[i]){continue;}

	    if(type == diag_orb_type_fo) {
		datavec[id] = list->fo.rdot;
	    }
	    list = list->prev;
	    id++;
	}

	H5LTmake_dataset(group, "v_R", 1, dims, H5T_IEEE_F64LE, datavec);
	H5LTset_attribute_string(group, "v_R", "unit", "m/s");

	free(datavec);
	return;
    }
    if(strcmp(dataset,"v_phi")==0) {
	real* datavec = (real*)malloc(datasize*sizeof(real));

	for(i = 0; i < size; i++) {
	    if(!mask[i]){continue;}

	    if(type == diag_orb_type_fo) {
		datavec[id] = list->fo.phidot * list->fo.r;
	    }
	    list = list->prev;
	    id++;
	}

	H5LTmake_dataset(group, "v_phi", 1, dims, H5T_IEEE_F64LE, datavec);
	H5LTset_attribute_string(group, "v_phi", "unit", "m/s");

	free(datavec);
	return;
    }
    if(strcmp(dataset,"v_z")==0) {
	real* datavec = (real*)malloc(datasize*sizeof(real));

	for(i = 0; i < size; i++) {
	    if(!mask[i]){continue;}

	    if(type == diag_orb_type_fo) {
		datavec[id] = list->fo.zdot;
	    }
	    list = list->prev;
	    id++;
	}

	H5LTmake_dataset(group, "v_z", 1, dims, H5T_IEEE_F64LE, datavec);
	H5LTset_attribute_string(group, "v_z", "unit", "m/s");

	free(datavec);
	return;
    }

    /* Guiding center specific */
    if(strcmp(dataset,"vpar")==0) {
	real* datavec = (real*)malloc(datasize*sizeof(real));

	for(i = 0; i < size; i++) {
	    if(!mask[i]){continue;}

	    if(type == diag_orb_type_gc) {
		datavec[id] = list->gc.vpar;
	    }
	    list = list->prev;
	    id++;
	}

	H5LTmake_dataset(group, "vpar", 1, dims, H5T_IEEE_F64LE, datavec);
	H5LTset_attribute_string(group, "vpar", "unit", "m/s");

	free(datavec);	
	return;
    }
    if(strcmp(dataset,"mu")==0) {
	real* datavec = (real*)malloc(datasize*sizeof(real));

	for(i = 0; i < size; i++) {
	    if(!mask[i]){continue;}

	    if(type == diag_orb_type_gc) {
		datavec[id] = list->gc.mu/CONST_E;
	    }
	    list = list->prev;
	    id++;
	}

	H5LTmake_dataset(group, "mu", 1, dims, H5T_IEEE_F64LE, datavec);
	H5LTset_attribute_string(group, "mu", "unit", "eV/T");

	free(datavec);
	return;
    }
    if(strcmp(dataset,"theta")==0) {
	real* datavec = (real*)malloc(datasize*sizeof(real));

	for(i = 0; i < size; i++) {
	    if(!mask[i]){continue;}

	    if(type == diag_orb_type_gc) {
		datavec[id] = list->gc.theta;
	    }
	    list = list->prev;
	    id++;
	}

	H5LTmake_dataset(group, "theta", 1, dims, H5T_IEEE_F64LE, datavec);
	H5LTset_attribute_string(group, "theta", "unit", "rad");

	free(datavec);
	return;
    }
}

void hdf5_orbits_writeset_fo(hid_t group, diag_orb_dat* list, int size, int* mask) {
    hdf5_orbits_writeset(group, diag_orb_type_fo, list, size, mask, "R");
    hdf5_orbits_writeset(group, diag_orb_type_fo, list, size, mask, "phi");
    hdf5_orbits_writeset(group, diag_orb_type_fo, list, size, mask, "z");
    hdf5_orbits_writeset(group, diag_orb_type_fo, list, size, mask, "v_R");
    hdf5_orbits_writeset(group, diag_orb_type_fo, list, size, mask, "v_phi");
    hdf5_orbits_writeset(group, diag_orb_type_fo, list, size, mask, "v_z");
    hdf5_orbits_writeset(group, diag_orb_type_fo, list, size, mask, "time");
    hdf5_orbits_writeset(group, diag_orb_type_fo, list, size, mask, "id");
    hdf5_orbits_writeset(group, diag_orb_type_fo, list, size, mask, "weight");
    hdf5_orbits_writeset(group, diag_orb_type_fo, list, size, mask, "mass");
    hdf5_orbits_writeset(group, diag_orb_type_fo, list, size, mask, "charge");
    hdf5_orbits_writeset(group, diag_orb_type_fo, list, size, mask, "B_r");
    hdf5_orbits_writeset(group, diag_orb_type_fo, list, size, mask, "B_phi");
    hdf5_orbits_writeset(group, diag_orb_type_fo, list, size, mask, "B_z");
}

void hdf5_orbits_writeset_gc(hid_t group, diag_orb_dat* list, int size, int* mask) {
    hdf5_orbits_writeset(group, diag_orb_type_gc, list, size, mask, "R");
    hdf5_orbits_writeset(group, diag_orb_type_gc, list, size, mask, "phi");
    hdf5_orbits_writeset(group, diag_orb_type_gc, list, size, mask, "z");
    hdf5_orbits_writeset(group, diag_orb_type_gc, list, size, mask, "vpar");
    hdf5_orbits_writeset(group, diag_orb_type_gc, list, size, mask, "mu");
    hdf5_orbits_writeset(group, diag_orb_type_gc, list, size, mask, "theta");
    hdf5_orbits_writeset(group, diag_orb_type_gc, list, size, mask, "time");
    hdf5_orbits_writeset(group, diag_orb_type_gc, list, size, mask, "id");
    hdf5_orbits_writeset(group, diag_orb_type_gc, list, size, mask, "weight");
    hdf5_orbits_writeset(group, diag_orb_type_gc, list, size, mask, "mass");
    hdf5_orbits_writeset(group, diag_orb_type_gc, list, size, mask, "charge");
    hdf5_orbits_writeset(group, diag_orb_type_gc, list, size, mask, "B_r");
    hdf5_orbits_writeset(group, diag_orb_type_gc, list, size, mask, "B_phi");
    hdf5_orbits_writeset(group, diag_orb_type_gc, list, size, mask, "B_z");
}

void hdf5_orbits_writeset_ml(hid_t group, diag_orb_dat* list, int size, int* mask) {
    hdf5_orbits_writeset(group, diag_orb_type_ml, list, size, mask, "R");
    hdf5_orbits_writeset(group, diag_orb_type_ml, list, size, mask, "phi");
    hdf5_orbits_writeset(group, diag_orb_type_ml, list, size, mask, "z");
    hdf5_orbits_writeset(group, diag_orb_type_ml, list, size, mask, "time");
    hdf5_orbits_writeset(group, diag_orb_type_ml, list, size, mask, "id");
    hdf5_orbits_writeset(group, diag_orb_type_ml, list, size, mask, "B_r");
    hdf5_orbits_writeset(group, diag_orb_type_ml, list, size, mask, "B_phi");
    hdf5_orbits_writeset(group, diag_orb_type_ml, list, size, mask, "B_z");
}
