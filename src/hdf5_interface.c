/**
 * @file hdf5_interface.c
 * @brief HDF5 operations are accessed through here
 *
 * This module handles IO operations to HDF5 file. Accessing HDF5 files
 * from the main program should be done using this module.
 */
#define _XOPEN_SOURCE 700
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "ascot5.h"
#include "simulate.h"
#include "print.h"
#include "gitver.h"
#include "compiler_flags.h"
#include "hdf5_interface.h"
#include "hdf5io/hdf5_helpers.h"
#include "hdf5io/hdf5_options.h"
#include "hdf5io/hdf5_bfield.h"
#include "hdf5io/hdf5_plasma.h"
#include "hdf5io/hdf5_neutral.h"
#include "hdf5io/hdf5_efield.h"
#include "hdf5io/hdf5_wall.h"
#include "hdf5io/hdf5_boozer.h"
#include "hdf5io/hdf5_mhd.h"
#include "hdf5io/hdf5_marker.h"
#include "hdf5io/hdf5_state.h"
#include "hdf5io/hdf5_dist.h"
#include "hdf5io/hdf5_orbit.h"
#include "hdf5io/hdf5_transcoef.h"
#include "hdf5io/hdf5_asigma.h"
#include "hdf5io/hdf5_nbi.h"

int hdf5_get_active_qid(hid_t f, const char* group, char qid[11]);

/**
 * @brief Read and initialize input data
 *
 * This function reads input from HDF5 file, initializes it, and
 * allocates offload arrays and returns the pointers to them.
 *
 * @param sim pointer to simulation offload struct
 * @param input_active bitflags for input types to read
 * @param B_offload_array pointer to magnetic field offload array
 * @param E_offload_array pointer to electric field offload array
 * @param plasma_offload_array pointer to plasma data offload array
 * @param neutral_offload_array pointer to neutral data offload array
 * @param wall_offload_array pointer to wall offload array
 * @param wall_int_offload_array pointer to wall integer offload array
 * @param boozer_offload_array pointer to boozer offload array
 * @param mhd_offload_array pointer to mhd offload array
 * @param asigma_offload_array pointer to atomic offload array
 * @param nbi_offload_array pointer to neutral beam injector data offload array
 * @param p pointer to marker offload data
 * @param n_markers pointer to integer notating how many markers were read
 *
 * @return zero if reading and initialization succeeded
 */
int hdf5_interface_read_input(sim_offload_data* sim,
                              int input_active,
                              real** B_offload_array,
                              real** E_offload_array,
                              real** plasma_offload_array,
                              real** neutral_offload_array,
                              real** wall_offload_array,
                              int**  wall_int_offload_array,
                              real** boozer_offload_array,
                              real** mhd_offload_array,
                              real** asigma_offload_array,
                              real** nbi_offload_array,
                              input_particle** p,
                              int* n_markers){

    print_out(VERBOSE_IO, "\nReading and initializing input.\n");

    /* This init disables automatic error messages.
     * We want to generate our own that are more informative.*/
    hdf5_init();

    /* Check if requested HDF5 file exists and open it */
    print_out(VERBOSE_IO, "\nInput file is %s.\n", sim->hdf5_in);
    hid_t f = hdf5_open_ro(sim->hdf5_in);
    if(f < 0) {
        print_err("Error: File not found.");
        return 1;
    }

    /* Read active input from hdf5 and initialize */
    char qid[11];

    if(input_active & hdf5_input_options) {
        if(hdf5_find_group(f, "/options/")) {
            print_err("Error: No options in input file.");
            return 1;
        }
        print_out(VERBOSE_IO, "\nReading options input.\n");
        if(sim->qid_options[0] != '\0') {
            strcpy(qid, sim->qid_options);
        }
        else if( hdf5_get_active_qid(f, "/options/", qid) ) {
            print_err("Error: Active QID not declared.");
            return 1;
        }
        strcpy(sim->qid_options, qid);
        print_out(VERBOSE_IO, "Active QID is %s\n", qid);
        if( hdf5_options_read(f, sim, qid) ) {
            print_err("Error: Failed to initialize options.\n");
            return 1;
        }
        print_out(VERBOSE_IO, "Options read and initialized.\n");
    }

    if(input_active & hdf5_input_bfield) {
        if(hdf5_find_group(f, "/bfield/")) {
            print_err("Error: No magnetic field in input file.");
            return 1;
        }
        print_out(VERBOSE_IO, "\nReading magnetic field input.\n");
        if(sim->qid_bfield[0] != '\0') {
            strcpy(qid, sim->qid_bfield);
        }
        else if( hdf5_get_active_qid(f, "/bfield/", qid) ) {
            print_err("Error: Active QID not declared.");
            return 1;
        }
        strcpy(sim->qid_bfield, qid);
        print_out(VERBOSE_IO, "Active QID is %s\n", qid);
        if( hdf5_bfield_init_offload(f, &(sim->B_offload_data),
                                     B_offload_array, qid) ) {
            print_err("Error: Failed to initialize magnetic field.\n");
            return 1;
        }
        print_out(VERBOSE_IO, "Magnetic field read and initialized.\n");
    }

    if(input_active & hdf5_input_efield) {
        if(hdf5_find_group(f, "/efield/")) {
            print_err("Error: No electric field in input file.");
            return 1;
        }
        print_out(VERBOSE_IO, "\nReading electric field input.\n");
        if(sim->qid_efield[0] != '\0') {
            strcpy(qid, sim->qid_efield);
        }
        else if( hdf5_get_active_qid(f, "/efield/", qid) ) {
            print_err("Error: Active QID not declared.");
            return 1;
        }
        strcpy(sim->qid_efield, qid);
        print_out(VERBOSE_IO, "Active QID is %s\n", qid);
        if( hdf5_efield_init_offload(f, &(sim->E_offload_data),
                                     E_offload_array, qid) ) {
            print_err("Error: Failed to initialize electric field.\n");
            return 1;
        }
        print_out(VERBOSE_IO, "Electric field read and initialized.\n");
    }

    if(input_active & hdf5_input_plasma) {
        if(hdf5_find_group(f, "/plasma/")) {
            print_err("Error: No plasma data in input file.");
            return 1;
        }
        print_out(VERBOSE_IO, "\nReading plasma input.\n");
        if(sim->qid_plasma[0] != '\0') {
            strcpy(qid, sim->qid_plasma);
        }
        else if( hdf5_get_active_qid(f, "/plasma/", qid) ) {
            print_err("Error: Active QID not declared.");
            return 1;
        }
        strcpy(sim->qid_plasma, qid);
        print_out(VERBOSE_IO, "Active QID is %s\n", qid);
        if( hdf5_plasma_init_offload(f, &(sim->plasma_offload_data),
                                     plasma_offload_array, qid) ) {
            print_err("Error: Failed to initialize plasma data.\n");
            return 1;
        }
        print_out(VERBOSE_IO, "Plasma data read and initialized.\n");
    }

    if(input_active & hdf5_input_neutral) {
        if(hdf5_find_group(f, "/neutral/")) {
            print_err("Error: No neutral data in input file.");
            return 1;
        }
        print_out(VERBOSE_IO, "\nReading neutral input.\n");
        if(sim->qid_neutral[0] != '\0') {
            strcpy(qid, sim->qid_neutral);
        }
        else if( hdf5_get_active_qid(f, "/neutral/", qid) ) {
            print_err("Error: Active QID not declared.");
            return 1;
        }
        strcpy(sim->qid_neutral, qid);
        print_out(VERBOSE_IO, "Active QID is %s\n", qid);
        if( hdf5_neutral_init_offload(f, &(sim->neutral_offload_data),
                                      neutral_offload_array, qid) ) {
            print_err("Error: Failed to initialize neutral data.\n");
            return 1;
        }
        print_out(VERBOSE_IO, "Neutral data read and initialized.\n");
    }

    if(input_active & hdf5_input_wall) {
        if(hdf5_find_group(f, "/wall/")) {
            print_err("Error: No wall data in input file.");
            return 1;
        }
        print_out(VERBOSE_IO, "\nReading wall input.\n");
        if(sim->qid_wall[0] != '\0') {
            strcpy(qid, sim->qid_wall);
        }
        else if( hdf5_get_active_qid(f, "/wall/", qid) ) {
            print_err("Error: Active QID not declared.");
            return 1;
        }
        strcpy(sim->qid_wall, qid);
        print_out(VERBOSE_IO, "Active QID is %s\n", qid);
        if( hdf5_wall_init_offload(f, &(sim->wall_offload_data),
                                   wall_offload_array, wall_int_offload_array,
                                   qid) ) {
            print_err("Error: Failed to initialize wall.\n");
            return 1;
        }
        print_out(VERBOSE_IO, "Wall data read and initialized.\n");
    }

    if(input_active & hdf5_input_boozer) {
        if(hdf5_find_group(f, "/boozer/")) {
            print_err("Error: No boozer data in input file.");
            return 1;
        }
        print_out(VERBOSE_IO, "\nReading boozer input.\n");
        if(sim->qid_boozer[0] != '\0') {
            strcpy(qid, sim->qid_boozer);
        }
        else if( hdf5_get_active_qid(f, "/boozer/", qid) ) {
            print_err("Error: Active QID not declared.");
            return 1;
        }
        strcpy(sim->qid_boozer, qid);
        print_out(VERBOSE_IO, "Active QID is %s\n", qid);
        if( hdf5_boozer_init_offload(f, &(sim->boozer_offload_data),
                                     boozer_offload_array, qid) ) {
            print_err("Error: Failed to read boozer input.\n");
            return 1;
        }
        print_out(VERBOSE_IO, "Boozer data read and initialized.\n");
    }

    if(input_active & hdf5_input_mhd) {
        if(hdf5_find_group(f, "/mhd/")) {
            print_err("Error: No MHD data in input file.");
            return 1;
        }
        print_out(VERBOSE_IO, "\nReading MHD input.\n");
        if(sim->qid_mhd[0] != '\0') {
            strcpy(qid, sim->qid_mhd);
        }
        else if( hdf5_get_active_qid(f, "/mhd/", qid) ) {
            print_err("Error: Active QID not declared.");
            return 1;
        }
        strcpy(sim->qid_mhd, qid);
        print_out(VERBOSE_IO, "Active QID is %s\n", qid);
        if( hdf5_mhd_init_offload(f, &(sim->mhd_offload_data),
                                  mhd_offload_array, qid) ) {
            print_err("Error: Failed to read MHD input.\n");
            return 1;
        }
        print_out(VERBOSE_IO, "MHD data read and initialized.\n");
    }

    if(input_active & hdf5_input_asigma) {
        if(hdf5_find_group(f, "/asigma/")) {
            print_err("Error: No atomic reaction data in input file.");
            return 1;
        }
        print_out(VERBOSE_IO, "\nReading atomic reaction input.\n");
        if(sim->qid_asigma[0] != '\0') {
            strcpy(qid, sim->qid_asigma);
        }
        else if( hdf5_get_active_qid(f, "/asigma/", qid) ) {
            print_err("Error: Active QID not declared.");
            return 1;
        }
        strcpy(sim->qid_asigma, qid);
        print_out(VERBOSE_IO, "Active QID is %s\n", qid);
        if( hdf5_asigma_init_offload(f, &(sim->asigma_offload_data),
                                     asigma_offload_array, qid) ) {
            print_err("Error: Failed to initialize atomic reaction data.\n");
            return 1;
        }
        print_out(VERBOSE_IO, "Atomic reaction data read and initialized.\n");
    }

    if(input_active & hdf5_input_nbi) {
        if(hdf5_find_group(f, "/nbi/")) {
            print_err("Error: No NBI data in input file.");
            return 1;
        }
        print_out(VERBOSE_IO, "\nReading NBI input.\n");
        if(sim->qid_nbi[0] != '\0') {
            strcpy(qid, sim->qid_nbi);
        }
        else if( hdf5_get_active_qid(f, "/nbi/", qid) ) {
            print_err("Error: Active QID not declared.");
            return 1;
        }
        strcpy(sim->qid_nbi, qid);
        print_out(VERBOSE_IO, "Active QID is %s\n", qid);
        if( hdf5_nbi_init_offload(f, &(sim->nbi_offload_data),
                                  nbi_offload_array, qid) ) {
            print_err("Error: Failed to initialize NBI data.\n");
            return 1;
        }
        print_out(VERBOSE_IO, "NBI data read and initialized.\n");
    }

    if(input_active & hdf5_input_marker) {
        if(hdf5_find_group(f, "/marker/")) {
            print_err("Error: No marker data in input file.");
            return 1;
        }
        print_out(VERBOSE_IO, "\nReading marker input.\n");
        if(sim->qid_marker[0] != '\0') {
            strcpy(qid, sim->qid_marker);
        }
        else if( hdf5_get_active_qid(f, "/marker/", qid) ) {
            print_err("Error: Active QID not declared.");
            return 1;
        }
        strcpy(sim->qid_marker, qid);
        print_out(VERBOSE_IO, "Active QID is %s\n", qid);
        if( hdf5_marker_read(f, n_markers, p, qid) ) {
            print_err("Error: Failed to read markers.\n");
            return 1;
        }
        print_out(VERBOSE_IO, "Marker data read and initialized.\n");
    }

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
 * /results/\<run\>_XXXXXXXXXX/ where X's are the qid of current run and \<run\>
 * is the type of the run: "run" for ascot5_main and "bbnbi" for bbnbi5.
 *
 * The group is initialized by writing qids of all used inputs as string
 * attributes in the run group. Also the date and empty "details" fields
 * are written.
 *
 * @param sim pointer to simulation offload struct
 * @param qid qid of this run
 * @param run type of this run
 *
 * @return Zero if initialization succeeded
 */
int hdf5_interface_init_results(sim_offload_data* sim, char* qid, char* run) {

    /* Create new file for the output if one does not yet exist. */
    hid_t fout = hdf5_create(sim->hdf5_out);
    if(fout < 0) {
        print_out(VERBOSE_IO, "Note: Output file %s is already present.\n",
                  sim->hdf5_out);
    }
    else {
        hdf5_close(fout);
    }

    /* Open output file. */
    fout = hdf5_open(sim->hdf5_out);
    if(fout < 0) {
        print_err("Error: Output file %s doesn't exists.\n", sim->hdf5_out);
        return 1;
    }

    /* Create a run group for this specific run (and results group if one */
    /* doesn't exist already.                                             */
    char runpath[200], path[256];
    sprintf(runpath, "/results/%s_XXXXXXXXXX", run);
    hdf5_gen_path(runpath, qid, path);
    hid_t newgroup = hdf5_create_group(fout, path);

    /* If a run with identical qid exists, abort. */
    print_out(VERBOSE_IO, "\nThe qid of this run is %s\n", qid);
    if(newgroup < 0) {
        print_err("Error: A run with qid %s already exists.\n", qid);
        hdf5_close(fout);
        return 1;
    }
    H5Gclose(newgroup);

    /* Set this run as the active run. */
    hdf5_write_string_attribute(fout, "/results", "active",  qid);

    /* Open input file (if different file than the output) */
    hid_t fin = fout;
    if( strcmp(sim->hdf5_in, sim->hdf5_out) != 0 ) {
        fin = hdf5_open(sim->hdf5_in);
    }

    /* Record input QIDs from sim struct to output group. */
    if(sim->qid_options[0] != '\0') {
        hdf5_write_string_attribute(
            fout, path, "qid_options", sim->qid_options);
    }

    if(sim->qid_bfield[0] != '\0') {
        hdf5_write_string_attribute(
            fout, path, "qid_bfield", sim->qid_bfield);
    }

    if(sim->qid_efield[0] != '\0') {
        hdf5_write_string_attribute(
            fout, path, "qid_efield", sim->qid_efield);
    }

    if(sim->qid_plasma[0] != '\0') {
        hdf5_write_string_attribute(
            fout, path, "qid_plasma", sim->qid_plasma);
    }

    if(sim->qid_neutral[0] != '\0') {
        hdf5_write_string_attribute(
            fout, path, "qid_neutral", sim->qid_neutral);
    }

    if(sim->qid_wall[0] != '\0') {
        hdf5_write_string_attribute(
            fout, path, "qid_wall", sim->qid_wall);
    }

    if(sim->qid_marker[0] != '\0') {
        hdf5_write_string_attribute(
            fout, path, "qid_marker", sim->qid_marker);
    }

    if(sim->qid_boozer[0] != '\0') {
        hdf5_write_string_attribute(
            fout, path, "qid_boozer", sim->qid_boozer);
    }

    if(sim->qid_mhd[0] != '\0') {
        hdf5_write_string_attribute(
            fout, path, "qid_mhd", sim->qid_mhd);
    }

    if(sim->qid_asigma[0] != '\0') {
        hdf5_write_string_attribute(
            fout, path, "qid_asigma", sim->qid_asigma);
    }

    if(sim->qid_nbi[0] != '\0') {
        hdf5_write_string_attribute(
            fout, path, "qid_nbi", sim->qid_nbi);
    }

    /* If input and output are different files, close input */
    if( strcmp(sim->hdf5_in, sim->hdf5_out) != 0 ) {
        hdf5_close(fin);
    }

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
    sprintf(date, "%04hu-%02hu-%02hu %02hu:%02hu:%02hu.", tm.tm_year + 1900,
            tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    hdf5_write_string_attribute(fout, path, "date",  date);

    /* Write compiler flags from stringified macros */
    hdf5_write_string_attribute(fout, path, "CFLAGS", str_macro(CFLAGS));
    hdf5_write_string_attribute(fout, path, "CC", str_macro(CC));

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
 * @return Zero if state was written succesfully
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
    char run[256];
    run[0] = '\0';
    if(run[0] == '\0') {
        /* Check if this an ascot5_main run */
        sprintf(run, "/results/run_%s/", qid);
        if( hdf5_find_group(f, run) < 0 ) {
            run[0] = '\0';
        }
    }
    if(run[0] == '\0') {
        /* Check if this an bbnbi5 run */
        sprintf(run, "/results/bbnbi_%s/", qid);
        if( hdf5_find_group(f, run) < 0 ) {
            run[0] = '\0';
        }
    }
    if(run[0] == '\0') {
        /* Check if this an afsi5 run */
        sprintf(run, "/results/afsi_%s/", qid);
        if( hdf5_find_group(f, run) < 0 ) {
            run[0] = '\0';
        }
    }
    if(run[0] == '\0') {
        print_err("Error: Run group not found.\n");
    }

    if( hdf5_state_write(f, run, state, n, p) ) {
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
    char path[256]; /* For storing dataset names */
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
    char run[256];
    run[0] = '\0';
    if(run[0] == '\0') {
        /* Check if this an ascot5_main run */
        sprintf(run, "/results/run_%s/", qid);
        if( hdf5_find_group(f, run) < 0 ) {
            run[0] = '\0';
        }
    }
    if(run[0] == '\0') {
        /* Check if this an bbnbi5 run */
        sprintf(run, "/results/bbnbi_%s/", qid);
        if( hdf5_find_group(f, run) < 0 ) {
            run[0] = '\0';
        }
    }
    if(run[0] == '\0') {
        /* Check if this an afsi5 run */
        sprintf(run, "/results/afsi_%s/", qid);
        if( hdf5_find_group(f, run) < 0 ) {
            run[0] = '\0';
        }
    }
    if(run[0] == '\0') {
        print_err("Error: Run group not found.\n");
    }

    if(sim->diag_offload_data.dist5D_collect) {
        print_out(VERBOSE_IO, "\nWriting 5D distribution.\n");
        int idx = sim->diag_offload_data.offload_dist5D_index;
        sprintf(path, "%sdist5d", run);
        if( hdf5_dist_write_5D(f, path, &sim->diag_offload_data.dist5D,
                               &diag_offload_array[idx]) ) {
            print_err("Warning: 5D distribution could not be written.\n");
        }
    }

    if(sim->diag_offload_data.dist6D_collect) {
        print_out(VERBOSE_IO, "\nWriting 6D distribution.\n");
        int idx = sim->diag_offload_data.offload_dist6D_index;
        sprintf(path, "%sdist6d", run);
        if( hdf5_dist_write_6D(f, path, &sim->diag_offload_data.dist6D,
                               &diag_offload_array[idx]) ) {
            print_err("Warning: 6D distribution could not be written.\n");
        }
    }
    if(sim->diag_offload_data.distrho5D_collect) {
        print_out(VERBOSE_IO, "\nWriting rho 5D distribution.\n");
        int idx = sim->diag_offload_data.offload_distrho5D_index;
        sprintf(path, "%sdistrho5d", run);
        if( hdf5_dist_write_rho5D(f, path, &sim->diag_offload_data.distrho5D,
                                  &diag_offload_array[idx]) ) {
            print_err("Warning: rho 5D distribution could not be written.\n");
        }
    }

    if(sim->diag_offload_data.distrho6D_collect) {
        print_out(VERBOSE_IO, "\nWriting rho 6D distribution.\n");
        int idx = sim->diag_offload_data.offload_distrho6D_index;
        sprintf(path, "%sdistrho6d", run);
        if( hdf5_dist_write_rho6D(f, path, &sim->diag_offload_data.distrho6D,
                                  &diag_offload_array[idx]) ) {
            print_err("Warning: rho 6D distribution could not be written.\n");
        }
    }

    if(sim->diag_offload_data.distCOM_collect) {
        print_out(VERBOSE_IO, "\nWriting COM distribution.\n");
        int idx = sim->diag_offload_data.offload_distCOM_index;
        sprintf(path, "%sdistcom", run);
        if( hdf5_dist_write_COM( f, path, &sim->diag_offload_data.distCOM,
                                 &diag_offload_array[idx]) ) {
            print_err("Warning: COM distribution could not be written.\n");
        }
    }

    if(sim->diag_offload_data.diagorb_collect) {
        print_out(VERBOSE_IO, "Writing orbit diagnostics.\n");
        int idx = sim->diag_offload_data.offload_diagorb_index;
        sprintf(path, "%sorbit", run);
        if( hdf5_orbit_write(f, path, &sim->diag_offload_data.diagorb,
                             &diag_offload_array[idx]) ) {
            print_err("Warning: Orbit diagnostics could not be written.\n");
        }
    }

    if(sim->diag_offload_data.diagtrcof_collect) {
        print_out(VERBOSE_IO, "Writing transport coefficient diagnostics.\n");
        int idx = sim->diag_offload_data.offload_diagtrcof_index;
        sprintf(path, "%stranscoef", run);
        if( hdf5_transcoef_write(f, path, &sim->diag_offload_data.diagtrcof,
                                 &diag_offload_array[idx]) ) {
            print_err("Warning: Coefficients could not be written.\n");
        }
    }

    hdf5_close(f);

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

/**
 * @brief Generate an identification number for a run
 *
 * The identification number (QID) is a 32 bit unsigned integer represented in a
 * string format, i.e., by ten characters. QID is a random integer between 0 and
 * 4 294 967 295, and it is padded with leading zeroes in string representation.
 *
 * @param qid a pointer to 11 chars wide array where generated QID is stored
 */
void hdf5_generate_qid(char* qid) {

    /* Seed random number generator with current time */
    struct timespec ts;
#ifdef __MACH__
    clock_gettime(CLOCK_MONOTONIC_RAW, &ts);
#else
    clock_gettime(CLOCK_REALTIME, &ts);
#endif
    srand48( ts.tv_nsec );

    /* Generate a 32 bit random integer by generating signed 32 bit random
     * integers with mrand48() and choosing the first one that is positive */
    long int qint = -1;
    while(qint < 0) {
        qint = mrand48();
    }

    /* Convert the random number to a string format */
    sprintf(qid, "%010lu", (long unsigned int)qint);
}
