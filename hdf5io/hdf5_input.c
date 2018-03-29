#include <stdio.h>
#include <string.h>
#include <time.h>
#include "../simulate.h"
#include <hdf5.h>
#include "hdf5_helpers.h"
#include "hdf5_hl.h"
#include "hdf5_input.h"
#include "hdf5_simulate.h"
#include "hdf5_bfield.h"
#include "hdf5_plasma.h"
#include "hdf5_neutral.h"
#include "hdf5_efield.h"
#include "hdf5_wall.h"
#include "hdf5_markers.h"

int hdf5_input(sim_offload_data* sim, 
	       real** B_offload_array,
	       real** E_offload_array,
	       real** plasma_offload_array,
	       real** neutral_offload_array,
	       real** wall_offload_array,
               input_particle** p,
               int* n_markers){
    
    /* This init disables automatic error messages.
     * We want to generate our own that are more informative.*/
    hdf5_init();
    herr_t err;


    /* Check if requested hdf5 exists and open it */
    hid_t f = hdf5_open(sim->hdf5_in);
    if(f < 0) {
	printf("\nError: Attempted to access input file %s.\n",sim->hdf5_in);
	printf("File not found.\n");
	return -1;
    }

    /* Check that input contains all relevant groups */
    err = hdf5_find_group(f, "/options/");
    if(err < 0) {
	printf("\nError: Options not found within %s.\n",sim->hdf5_in);
	return -1;
    }

    err = hdf5_find_group(f, "/bfield/");
    if(err < 0) {
	printf("\nError: Magnetic field not found within %s.\n",sim->hdf5_in);
	return -1;
    }

    err = hdf5_find_group(f, "/efield/");
    if(err < 0) {
	printf("\nError: Electric field not found within %s.\n",sim->hdf5_in);
	return -1;
    }

    err = hdf5_find_group(f, "/plasma/");
    if(err < 0) {
	printf("\nError: Plasma not found within %s.\n",sim->hdf5_in);
	return -1;
    }

    err = hdf5_find_group(f, "/wall/");
    if(err < 0) {
	printf("\nError: Wall not found within %s.\n",sim->hdf5_in);
	return -1;
    }

    err = hdf5_find_group(f, "/marker/");
    if(err < 0) {
	printf("\nError: Marker input not found within %s.\n",sim->hdf5_in);
	return -1;
    }

    /* Read input from hdf5 and initialize */
    hdf5_simulate(f, sim);
    if(err < 0) {
	printf("\nError: Failed to initialize options from %s.\n",sim->hdf5_in);
	return -1;
    }

    hdf5_bfield_init_offload(f,&(sim->B_offload_data), B_offload_array);
    if(err < 0) {
	printf("\nError: Failed to initialize magnetic field from %s.\n",sim->hdf5_in);
	return -1;
    }

    hdf5_efield_init_offload(f,&(sim->E_offload_data), E_offload_array);
    if(err < 0) {
	printf("\nError: Failed to initialize electric field from %s.\n",sim->hdf5_in);
	return -1;
    }

    hdf5_plasma_init_offload(f,&(sim->plasma_offload_data), plasma_offload_array);
    if(err < 0) {
	printf("\nError: Failed to initialize plasma from %s.\n",sim->hdf5_in);
	return -1;
    }

    hdf5_neutral_init_offload(f,&(sim->neutral_offload_data), neutral_offload_array);
    if(err < 0) {
	printf("\nError: Failed to initialize neutral density from %s.\n",sim->hdf5_in);
	return -1;
    }

    hdf5_wall_init_offload(f,&(sim->wall_offload_data), wall_offload_array);
    if(err < 0) {
	printf("\nError: Failed to initialize wall from %s.\n",sim->hdf5_in);
	return -1;
    }

    hdf5_markers_init(f, n_markers, p);
    if(err < 0) {
	printf("\nError: Failed to initialize markers from %s.\n",sim->hdf5_in);
	return -1;
    }

    /* Close the hdf5 file */
    err = hdf5_close(f);
    if(err < 0) {
	printf("\nError: Could not close hdf5 file %s.\n",sim->hdf5_in);
	return -1;
    }

    return 0;
}

int hdf5_initoutput(sim_offload_data* sim, char* qid) {
    
    /* Create new file for the output if one does not yet exist. */
    hid_t fout = hdf5_create(sim->hdf5_out);
    if(fout < 0) {
	printf("Note: Output file %s is already present.\n", sim->hdf5_out);
    }
    hdf5_close(fout);
    
    /* Open output file and create results section if one does not yet exist. */
    fout = hdf5_open(sim->hdf5_out);
    
    herr_t  err = hdf5_find_group(fout, "/results/");
    if(err) {
	/* No results group exists so we make one */
	hdf5_create_group(fout, "/results/");
    }

    /* Create a group for this specific run. */
    char path[256];
    hdf5_generate_qid_path("/results/run-XXXXXXXXXX", qid, path);
    hid_t newgroup = H5Gcreate2(fout, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose (newgroup);
    #if VERBOSE > 0
        printf("\nThe qid of this run is %s\n", qid);
    #endif
	
    /* If a run with identical qid exists, abort. */
    if(newgroup < 0) {
	printf("Error: A run with qid %s already exists.\n", qid);
	hdf5_close(fout);
	return -1;
    }

    /* Set this run as the active result. */
    H5LTset_attribute_string(fout, "/results", "active", qid);
    
    /* Read input data qids and store them here. */
    char inputqid[11];
    inputqid[10] = '\0';

    hid_t fin = hdf5_open(sim->hdf5_in);
    
    H5LTget_attribute_string(fin, "/options/", "active", inputqid);
    H5LTset_attribute_string(fout, path, "qid_options", inputqid);

    H5LTget_attribute_string(fin, "/bfield/", "active", inputqid);
    H5LTset_attribute_string(fout, path, "qid_bfield", inputqid);

    H5LTget_attribute_string(fin, "/efield/", "active", inputqid);
    H5LTset_attribute_string(fout, path, "qid_efield", inputqid);

    H5LTget_attribute_string(fin, "/plasma/", "active", inputqid);
    H5LTset_attribute_string(fout, path, "qid_plasma", inputqid);

    H5LTget_attribute_string(fin, "/wall/", "active", inputqid);
    H5LTset_attribute_string(fout, path, "qid_wall", inputqid);

    H5LTget_attribute_string(fin, "/marker/", "active", inputqid);
    H5LTset_attribute_string(fout, path, "qid_marker", inputqid);

    hdf5_close(fin);
    
    /* Finally we set a description and date, and close the file. */
    H5LTset_attribute_string(fout, path, "description", "");
    
    time_t t = time(NULL);
    struct tm tm = *localtime(&t);
    char date[21];
    sprintf(date, "%04d-%02d-%02d %02d:%02d:%02d.", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    H5LTset_attribute_string(fout, path, "date", date);
    hdf5_close(fout);
    
    return 0;
}
