/**
 * @file hdf5_interface.c
 * @brief HDF5 operations are accessed through here
 *
 * This module handles IO operations to HDF5 file. Accessing HDF5 files
 * from the main program should be done using this module.
 */
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "ascot5.h"
#include "simulate.h"
#include "print.h"
#include "gitver.h"
#include "hdf5_interface.h"
#include "hdf5io/hdf5_helpers.h"
#include "hdf5io/hdf5_options.h"
#include "hdf5io/hdf5_bfield.h"
#include "hdf5io/hdf5_plasma.h"
#include "hdf5io/hdf5_neutral.h"
#include "hdf5io/hdf5_efield.h"
#include "hdf5io/hdf5_wall.h"
#include "hdf5io/hdf5_markers.h"
#include "hdf5io/hdf5_particlestate.h"
#include "hdf5io/hdf5_dist.h"
#include "hdf5io/hdf5_orbits.h"

/**
 * @brief Read and initialize input data
 *
 * This function reads input from HDF5 file, initializes it, and
 * allocates offload arrays and returns the pointers to them.
 *
 * @param sim pointer to simulation offload struct
 * @param B_offload_array pointer to magnetic field offload array
 * @param E_offload_array pointer to electric field offload array
 * @param plasma_offload_array pointer to plasma data offload array
 * @param neutral_offload_array pointer to neutral data offload array
 * @param wall_offload_array pointer to wall offload array
 * @param p pointer to marker offload data
 * @param n_markers pointer to integer notating how many markers were read
 *
 * @return zero if reading and initialization succeeded
 */
int hdf5_interface_read_input(sim_offload_data* sim,
                              real** B_offload_array,
                              real** E_offload_array,
                              real** plasma_offload_array,
                              real** neutral_offload_array,
                              real** wall_offload_array,
                              input_particle** p,
                              int* n_markers){

    print_out(VERBOSE_IO, "\nReading and initializing input.\n");

    /* This init disables automatic error messages.
     * We want to generate our own that are more informative.*/
    hdf5_init();

    /* Check if requested HDF5 file exists and open it */
    print_out(VERBOSE_IO, "\nInput file is %s.\n", sim->hdf5_in);
    hid_t f = hdf5_open(sim->hdf5_in);
    if(f < 0) {
        print_err("Error: File not found.");
        return 1;
    }

    /* Check that input contains all relevant groups */

    if(hdf5_find_group(f, "/options/")) {
        print_err("Error: No options in input file.");
        return 1;
    }

    if(hdf5_find_group(f, "/bfield/")) {
        print_err("Error: No magnetic field in input file.");
        return 1;
    }

    if(hdf5_find_group(f, "/efield/")) {
        print_err("Error: No electric field in input file.");
        return 1;
    }

    if(hdf5_find_group(f, "/plasma/")) {
        print_err("Error: No plasma data in input file.");
        return 1;
    }

    if(hdf5_find_group(f, "/neutral/")) {
        print_err("Error: No neutral data in input file.");
        return 1;
    }

    if(hdf5_find_group(f, "/wall/")) {
        print_err("Error: No wall data in input file.");
        return 1;
    }

    if(hdf5_find_group(f, "/marker/")) {
        print_err("Error: No marker data in input file.");
        return 1;
    }

    /* Read input from hdf5 and initialize */
    char qid[11];


    print_out(VERBOSE_IO, "\nReading options input.\n");
    if( hdf5_get_active_qid(f, "/options/", qid) ) {
        print_err("Error: Active QID not declared.");
        return 1;
    }
    print_out(VERBOSE_IO, "Active QID is %s\n", qid);
    if( hdf5_options_read(f, sim, qid) ) {
        print_err("Error: Failed to initialize options.\n");
        return 1;
    }
    print_out(VERBOSE_IO, "Options read and initialized.\n");


    print_out(VERBOSE_IO, "\nReading magnetic field input.\n");
    if( hdf5_get_active_qid(f, "/bfield/", qid) ) {
        print_err("Error: Active QID not declared.");
        return 1;
    }
    print_out(VERBOSE_IO, "Active QID is %s\n", qid);
    if( hdf5_bfield_init_offload(f, &(sim->B_offload_data),
                                 B_offload_array, qid) ) {
        print_err("Error: Failed to initialize magnetic field.\n");
        return 1;
    }
    print_out(VERBOSE_IO, "Magnetic field read and initialized.\n");


    print_out(VERBOSE_IO, "\nReading electric field input.\n");
    if( hdf5_get_active_qid(f, "/efield/", qid) ) {
        print_err("Error: Active QID not declared.");
        return 1;
    }
    print_out(VERBOSE_IO, "Active QID is %s\n", qid);
    if( hdf5_efield_init_offload(f, &(sim->E_offload_data),
                                 E_offload_array, qid) ) {
        print_err("Error: Failed to initialize electric field.\n");
        return 1;
    }
    print_out(VERBOSE_IO, "Electric field read and initialized.\n");


    print_out(VERBOSE_IO, "\nReading plasma input.\n");
    if( hdf5_get_active_qid(f, "/plasma/", qid) ) {
        print_err("Error: Active QID not declared.");
        return 1;
    }
    print_out(VERBOSE_IO, "Active QID is %s\n", qid);
    if( hdf5_plasma_init_offload(f, &(sim->plasma_offload_data),
                                 plasma_offload_array, qid) ) {
        print_err("Error: Failed to initialize plasma data.\n");
        return 1;
    }
    print_out(VERBOSE_IO, "Plasma data read and initialized.\n");


    print_out(VERBOSE_IO, "\nReading neutral input.\n");
    if( hdf5_get_active_qid(f, "/neutral/", qid) ) {
        print_err("Error: Active QID not declared.");
        return 1;
    }
    print_out(VERBOSE_IO, "Active QID is %s\n", qid);
    if( hdf5_neutral_init_offload(f, &(sim->neutral_offload_data),
                                  neutral_offload_array, qid) ) {
        print_err("Error: Failed to initialize neutral data.\n");
        return 1;
    }
    print_out(VERBOSE_IO, "Neutral data read and initialized.\n");


    print_out(VERBOSE_IO, "\nReading wall input.\n");
    if( hdf5_get_active_qid(f, "/wall/", qid) ) {
        print_err("Error: Active QID not declared.");
        return 1;
    }
    print_out(VERBOSE_IO, "Active QID is %s\n", qid);
    if( hdf5_wall_init_offload(f, &(sim->wall_offload_data),
                               wall_offload_array, qid) ) {
        print_err("Error: Failed to initialize wall.\n");
        return 1;
    }
    print_out(VERBOSE_IO, "Wall data read and initialized.\n");


    print_out(VERBOSE_IO, "\nReading marker input.\n");
    if( hdf5_get_active_qid(f, "/marker/", qid) ) {
        print_err("Error: Active QID not declared.");
        return 1;
    }
    print_out(VERBOSE_IO, "Active QID is %s\n", qid);
    if( hdf5_markers_read(f, n_markers, p, qid) ) {
        print_err("Error: Failed to read markers.\n");
        return 1;
    }
    print_out(VERBOSE_IO, "Marker data read and initialized.\n");


    /* Close the hdf5 file */
    if( hdf5_close(f) ) {
        print_err("Error: Could not close the file.\n");
        return 1;
    }

    print_out(VERBOSE_IO, "\nAll input read and initialized.\n");
    return 0;
}

/**
 * @brief Initialize run group
 *
 * This functions creates results group (if one does not already exist) and
 * creates run group corresponding to this run. Run group is named as
 * /results/run-XXXXXXXXXX/ where X's are the qid of current run.
 *
 * The group is initialized by writing qids of all used inputs as string
 * attributes in the run group. Also the date and empty "details" fields
 * are written.
 *
 * @param sim pointer to simulation offload struct
 * @param qid qid of this run
 *
 * @return Zero if initialization succeeded
 */
int hdf5_interface_init_results(sim_offload_data* sim, char* qid) {

    /* Create new file for the output if one does not yet exist. */
    hid_t fout = hdf5_create(sim->hdf5_out);
    if(fout < 0) {
        print_out(VERBOSE_IO, "Note: Output file %s is already present.\n",
                  sim->hdf5_out);
    }
    hdf5_close(fout);

    /* Open output file and create results section if one does not yet exist. */
    fout = hdf5_open(sim->hdf5_out);

    if( hdf5_find_group(fout, "/results/") ) {
        hdf5_create_group(fout, "/results/");
    }

    /* Create a run group for this specific run. */
    char path[256];
    hdf5_gen_path("/results/run-XXXXXXXXXX", qid, path);
    hid_t newgroup = H5Gcreate2(fout, path,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose (newgroup);
    print_out(VERBOSE_IO, "\nThe qid of this run is %s\n", qid);

    /* If a run with identical qid exists, abort. */
    if(newgroup < 0) {
        print_err("Error: A run with qid %s already exists.\n", qid);
        hdf5_close(fout);
        return 1;
    }

    /* Set this run as the active run. */
    hdf5_write_string_attribute(fout, "/results", "active",  qid);

    /* Read input data qids and store them here. */
    char inputqid[11];
    inputqid[10] = '\0';

    hid_t fin = hdf5_open(sim->hdf5_in);

    H5LTget_attribute_string(fin, "/options/", "active", inputqid);
    hdf5_write_string_attribute(fout, path, "qid_options",  inputqid);

    H5LTget_attribute_string(fin, "/bfield/", "active", inputqid);
    hdf5_write_string_attribute(fout, path, "qid_bfield",  inputqid);

    H5LTget_attribute_string(fin, "/efield/", "active", inputqid);
    hdf5_write_string_attribute(fout, path, "qid_efield",  inputqid);

    H5LTget_attribute_string(fin, "/plasma/", "active", inputqid);
    hdf5_write_string_attribute(fout, path, "qid_plasma",  inputqid);

    H5LTget_attribute_string(fin, "/neutral/", "active", inputqid);
    hdf5_write_string_attribute(fout, path, "qid_neutral",  inputqid);

    H5LTget_attribute_string(fin, "/wall/", "active", inputqid);
    hdf5_write_string_attribute(fout, path, "qid_wall",  inputqid);

    H5LTget_attribute_string(fin, "/marker/", "active", inputqid);
    hdf5_write_string_attribute(fout, path, "qid_marker",  inputqid);

    hdf5_close(fin);

    /* Set a description, repository status, and date; close the file. */
    hdf5_write_string_attribute(fout, path, "description",  sim->description);

#ifdef GIT_VERSION
    char git[256];
    sprintf(git, "Tag %s Branch %s", GIT_VERSION, GIT_BRANCH);
    hdf5_write_string_attribute(fout, path, "repository", git);
#else
    hdf5_write_string_attribute(fout, path, "repository",
                                "Not under version control");
#endif

    time_t t = time(NULL);
    struct tm tm = *localtime(&t);
    char date[21];
    sprintf(date, "%04d-%02d-%02d %02d:%02d:%02d.", tm.tm_year + 1900,
            tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    hdf5_write_string_attribute(fout, path, "date",  date);
    hdf5_close(fout);

    return 0;
}

/**
 * @brief Write marker state to HDF5 output
 *
 * @param fn HDF5 output filename
 * @param state name of the state to be written
 * @param n number of markers in marker array
 * @param p array of markers to be written
 *
 * @results Zero if state was written succesfully
 */
int hdf5_interface_write_state(char* fn, char* state, integer n,
                             particle_state* p) {
    hid_t f = hdf5_open(fn);
    if(f < 0) {
        print_err("Error: File not found.\n");
        return 1;
    }

    char qid[11];
    if( hdf5_get_active_qid(f, "/results/", qid) ) {
        print_err("Error: Active QID was not written to results group.\n");
        hdf5_close(f);
        return 1;
    }

    if( hdf5_particlestate_write(f, qid, state, n, p) ) {
        print_err("Error: State could not be written.\n");
        hdf5_close(f);
        return 1;
    }

    hdf5_close(f);
    return 0;
}

/**
 * @brief Write diagnostics to HDF5 output
 *
 * @param sim pointer to simulation offload data
 * @param diag_offload_array diagnostics offload array
 * @param out HDF5 output filename
 *
 * @return Zero if diagnostics were written succesfully
 */
int hdf5_interface_write_diagnostics(sim_offload_data* sim,
                                     real* diag_offload_array, char* out) {

    print_out(VERBOSE_IO, "\nWriting diagnostics output.\n");

    hid_t f = hdf5_open(out);
    if(f < 0) {
        print_err("Error: File not found.\n");
        return 1;
    }

    char qid[11];
    if( hdf5_get_active_qid(f, "/results/", qid) ) {
        print_err("Error: Active QID was not written to results group.\n");
        hdf5_close(f);
        return 1;
    }
    hdf5_close(f);

    if(sim->diag_offload_data.dist5D_collect) {
        hdf5_dist_write_5D(
            &sim->diag_offload_data.dist5D,
            &diag_offload_array[sim->diag_offload_data.offload_dist5D_index],
            out, qid);
    }

    if(sim->diag_offload_data.dist6D_collect) {
        hdf5_dist_write_6D(
            &sim->diag_offload_data.dist6D,
            &diag_offload_array[sim->diag_offload_data.offload_dist6D_index],
            out, qid);
    }
    if(sim->diag_offload_data.distrho5D_collect) {
        hdf5_dist_write_rho5D(
            &sim->diag_offload_data.distrho5D,
            &diag_offload_array[sim->diag_offload_data.offload_distrho5D_index],
            out, qid);
    }

    if(sim->diag_offload_data.distrho6D_collect) {
        hdf5_dist_write_rho6D(
            &sim->diag_offload_data.distrho6D,
            &diag_offload_array[sim->diag_offload_data.offload_distrho6D_index],
            out, qid);
    }

    if(sim->diag_offload_data.diagorb_collect) {
        hid_t f = hdf5_open(out);
        print_out(VERBOSE_IO, "Writing orbit diagnostics.\n");

        int idx = sim->diag_offload_data.offload_diagorb_index;
        if( hdf5_orbits_write(f, qid, &sim->diag_offload_data.diagorb,
                              &diag_offload_array[idx]) ) {
            print_err("Warning: Orbit diagnostics could not be written.\n");
        }
        hdf5_close(f);
    }

    print_out(VERBOSE_IO, "\nDiagnostics output written.\n");

    return 0;
}

/**
 * @brief Fetch active qid within the given group
 *
 * Each input group (/bfield/, /options/, etc.) has qid string indicating which
 * of the subgroups is active, i.e., meant to be used within this simulation.
 *
 * This function fetches the active qid assuming the file is opened and closed
 * outside of this function.
 *
 * @param f HDF5 file identifier
 * @param group group string including the "/"s on both sides e.g. /bfield/
 * @param qid array where qid will be stored
 *
 * @return Zero on success
 */
int hdf5_get_active_qid(hid_t f, const char* group, char qid[11]) {
    if( H5LTget_attribute_string(f, group, "active", qid) ) {
        return 1;
    }
    qid[10] = '\0';

    return 0;
}
