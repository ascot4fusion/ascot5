/**
 * @file ascot5_main.c
 * @brief ASCOT5 stand-alone program
 *
 * This program reads data from input HDF5 file, simulates the given markers,
 * and writes the output data to a HDF5 file.
 *
 * The input and output files can be separate.
 *
 * Example:
 *
 *     ascot5_main --in=in --out=out
 *
 * Here "in" refers to in.h5 where input data is located and "out" to out.h5
 * where results will be stored. If no input argument is given the data is read
 * from ascot.h5. If not output argument is given the results are stored in the
 * input file.
 *
 * This program assumes that the input file contains magnetic field, electric
 * field, plasma, wall, and neutral data along with markers and options. See
 * hdf5_input.c for details. This program uses the input fields that are marked
 * as active (the HDF5 file can contain multiple instances of same input types
 * but only the active one is used here).
 *
 * The results are stored under /results/ group in output HDF5 file. The group
 * is created if one does not exists. For each run a specific "run" group is
 * created, which has the format run-XXXXXXXXXX, where "XXXXXXXXXX" is randomly
 * generated identification number (QID). The run group holds information when
 * the run was started, which input fields were used (referenced by their QIDs),
 * and at least the marker initial and end states if the simulation succeeded.
 * Also any other diagnostic data that was used is stored there.
 *
 * This program uses MPI by dividing the number of markers equally to all MPI
 * processes. The markers are not suffled so user is advised to do it beforehand
 * to ensure work is evenly distributed. A single MPI process can be simulated
 * with:
 *
 *     ascot5_main --mpi_size=size --mpi_rank=rank
 *
 * where size refers to number of MPI processes and rank is the process being
 * run (between [0, size-1]). Running the program this way does not use MPI.
 * This is intended to be used in Condor-like environments.
 *
 * You can add a description of the simulation as:
 *
 * ascot5_main --d="This is a test run"
 *
 * which is written in HDF5 file at the run group specific to this simulation.
 *
 * In addition to output data, the simulation progress may be written in
 * *.stdout files with each MPI process having dedicated file. See ascot5.h for
 * details.
 */
#include <getopt.h>
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "ascot5.h"
#include "consts.h"
#include "math.h"
#include "wall.h"
#include "diag.h"
#include "B_field.h"
#include "plasma.h"
#include "print.h"
#include "simulate.h"
#include "particle.h"
#include "endcond.h"
#include "hdf5_interface.h"
#include "offload.h"
#include "gitver.h"
#include "mpi_interface.h"

#include "ascot5_main.h"

#ifdef TRAP_FPE
#include <fenv.h>
#endif

int read_arguments(int argc, char** argv, sim_offload_data* sim);

/**
 * @brief Main function for ascot5_main
 *
 * This function calls functions that read input data from the disk, and
 * functions that initialize the offload data structs and offload arrays.
 * Actual simulation is done by calling simulate(). Once the simulation has been
 * completed, offload arrays are deallocated and the results are written to the
 * disk.
 *
 * MPI level parallelisation is done here as well as the offloading.
 *
 * @param  argc argument count of the command line arguments
 * @param  argv argument vector of the command line arguments
 *
 * @return Zero if simulation was completed with no issues
 */
int main(int argc, char** argv) {

#ifdef TRAP_FPE
    /* This will raise floating point exceptions */
    feenableexcept(FE_DIVBYZERO| FE_INVALID | FE_OVERFLOW);
    /*
     * If you are hunting a specific exception, you can disable the exceptions
     * in other parts of the code by surrounding it with: */
    /*
     * fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
     *  --- your  code here ---
     * feenableexcept(FE_DIVBYZERO  | FE_INVALID | FE_OVERFLOW);
     *
     * */
#endif

    /* Read and parse command line arguments */
    sim_offload_data sim;
    if( read_arguments(argc, argv, &sim) ) {
        abort();
        return 1;
    }

    if(sim.mpi_size > 0) {
        /* This is a pseudo-mpi run, where rank and size were set on the command
         * line. Only set root equal to rank since there are no other processes
         */
        sim.mpi_root = sim.mpi_rank;
    }
    else {
        /* Init MPI if used, or run serial */
        int mpi_rank, mpi_size, mpi_root;
        mpi_interface_init(argc, argv, &mpi_rank, &mpi_size, &mpi_root);
        sim.mpi_rank = mpi_rank;
        sim.mpi_size = mpi_size;
        sim.mpi_root = mpi_root;
    }

    print_out0(VERBOSE_MINIMAL, sim.mpi_rank, sim.mpi_root, "ASCOT5_MAIN\n");
    print_out0(VERBOSE_MINIMAL, sim.mpi_rank, sim.mpi_root,
               "Tag %s\nBranch %s\n\n", GIT_VERSION, GIT_BRANCH);
    print_out(VERBOSE_NORMAL, "Initialized MPI, rank %d, size %d.\n",
              sim.mpi_rank, sim.mpi_size);

    /* Total number of markers to be simulated */
    int n_tot;
    /* Marker input struct */
    input_particle* p;

    /* Offload data arrays that are allocated when input is read */
    real* B_offload_array;
    real* E_offload_array;
    real* plasma_offload_array;
    real* neutral_offload_array;
    real* wall_offload_array;
    int* wall_int_offload_array;
    real* boozer_offload_array;
    real* mhd_offload_array;
    real* asigma_offload_array;

    /* Read input from the HDF5 file */
    if( hdf5_interface_read_input(&sim,
                                  hdf5_input_options | hdf5_input_bfield |
                                  hdf5_input_efield  | hdf5_input_plasma |
                                  hdf5_input_neutral | hdf5_input_wall   |
                                  hdf5_input_marker  | hdf5_input_boozer |
                                  hdf5_input_mhd     | hdf5_input_asigma,
                                  &B_offload_array, &E_offload_array,
                                  &plasma_offload_array, &neutral_offload_array,
                                  &wall_offload_array,  &wall_int_offload_array,
                                  &boozer_offload_array, &mhd_offload_array,
                                  &asigma_offload_array, NULL,
                                  &p, &n_tot) ) {
        print_out0(VERBOSE_MINIMAL, sim.mpi_rank, sim.mpi_root,
                   "\nInput reading or initializing failed.\n"
                   "See stderr for details.\n");
        mpi_interface_finalize();
        abort();
        return 1;
    };

    /* Initialize marker states array ps and free marker input p */
    int n_proc; /* Number of markers allocated for this MPI process */
    particle_state* ps;
    if( prepare_markers(&sim, n_tot, p, &ps, &n_proc, B_offload_array) ) {
        goto CLEANUP_FAILURE;
    }

    /* Combine input offload arrays to one */
    offload_package offload_data;
    real* offload_array;
    int* int_offload_array;
    if( pack_offload_array(
            &sim, &offload_data, &B_offload_array, &E_offload_array,
            &plasma_offload_array, &neutral_offload_array, &wall_offload_array,
            &wall_int_offload_array, &boozer_offload_array, &mhd_offload_array,
            &asigma_offload_array, &offload_array, &int_offload_array) ) {
        goto CLEANUP_FAILURE;
    }

    /* Initialize diagnostics offload data */
    real* diag_offload_array;
    diag_init_offload(&sim.diag_offload_data, &diag_offload_array, n_tot);

    real diag_offload_array_size = sim.diag_offload_data.offload_array_length
        * sizeof(real) / (1024.0*1024.0);
    print_out0(VERBOSE_NORMAL, sim.mpi_rank, sim.mpi_root,
               "Initialized diagnostics, %.1f MB.\n", diag_offload_array_size);

    /* Write run group and inistate */
    char qid[11];
    hdf5_generate_qid(qid);
    if( write_rungroup(&sim, ps, n_tot, qid) ) {
        goto CLEANUP_FAILURE;
    }

    /* Actual simulation is done here; the results are stored in
     * diag_offload_array and pout contains end states from all processes.
     * n_gathered is the number of markers in the output array, which is n_tot
     * for most cases except when the simulation is run in condor-like manner,
     * in which case it is equal to n_proc. */
    int n_gathered;
    particle_state* pout;
    offload_and_simulate(
        &sim, n_tot, n_proc, ps, &offload_data, offload_array,
        int_offload_array, &n_gathered, &pout, diag_offload_array);

    /* Free input data */
    offload_free_offload(&offload_data, &offload_array, &int_offload_array);

    /* Write output and clean */
    if( write_output(&sim, pout, n_gathered, diag_offload_array) ) {
        goto CLEANUP_FAILURE;
    }
    diag_free_offload(&sim.diag_offload_data, &diag_offload_array);

    /* Display marker summary and free marker arrays */
    if(sim.mpi_rank == sim.mpi_root) {
        print_marker_summary(pout, n_tot);
    }
    free(pout);

    print_out0(VERBOSE_MINIMAL, sim.mpi_rank, sim.mpi_root, "\nDone.\n");
    mpi_interface_finalize();
    return 0;

/* GOTO this block to free resources in case simulation crashes */
CLEANUP_FAILURE:
    mpi_interface_finalize();
    free(p);
    free(ps);
    free(pout);
    free(offload_array);
    free(int_offload_array);
    free(diag_offload_array);
    abort();
    return 1;
}

/**
 * @brief Prepare markers for offload
 *
 * This function initializes the marker states and allocates memory for the
 * particle states. The input markers are read from the HDF5 file and
 * stored in the input_particle_states array. The initial marker states
 * are then calculated and stored in the particle_states array.
 *
 * When MPI is used, the marker states are initialized only for those markers
 * that are used in this MPI process.
 *
 * @param sim simulation offload data struct
 * @param n_tot total number of markers
 * @param pin pointer to marker input array which is deallocated here
 * @param pout pointer to marker state array created here
 * @param n_proc pointer to variable for number of markers for this process
 * @param B_offload_array pointer to magnetic field data needed for marker init
 *
 * @returns zero on success
 */
int prepare_markers(
    sim_offload_data* sim, int n_tot, input_particle* pin,
    particle_state** pout, int* n_proc, real* B_offload_array) {

    /* This sets up GC transformation order etc. so it must be called before
     * initializing markers. */
    simulate_init_offload(sim);

    /* Choose which markers are used in this MPI process. Simply put, markers
     * are divided into mpi_size sequential blocks and the mpi_rank:th block
     * is chosen for this simulation. */
    int start_index;
    mpi_my_particles(&start_index, n_proc, n_tot, sim->mpi_rank, sim->mpi_size);
    pin += start_index;

    /* Set up particlestates on host, needs magnetic field evaluation */
    print_out0(VERBOSE_NORMAL, sim->mpi_rank, sim->mpi_root,
               "\nInitializing marker states.\n");
    B_field_data Bdata;
    B_field_init(&Bdata, &sim->B_offload_data, B_offload_array);

    *pout = (particle_state*) malloc(*n_proc * sizeof(particle_state));
    for(int i = 0; i < *n_proc; i++) {
        particle_input_to_state(&(pin[i]), &((*pout)[i]), &Bdata);
    }

    print_out0(VERBOSE_NORMAL, sim->mpi_rank, sim->mpi_root,
               "Estimated memory usage %.1f MB.\n",
               (sizeof(real) * (*n_proc)) / (1024.0*1024.0));
    print_out0(VERBOSE_NORMAL, sim->mpi_rank, sim->mpi_root,
               "Marker states initialized.\n");

    /* We can now free the input markers */
    free(pin-start_index);

    return 0;
}


/**
 * @brief Prepare offload array to be offloaded
 *
 * When data is read, it is stored to input specific offload arrays which
 * are packed as a single array here (two actually as one array contains
 * integers and the other floats). The initial individual arrays are
 * deallocated.
 *
 * @param sim simulation offload data struct
 * @param offload_data empty offload package
 * @param B_offload_array magnetic field offload array
 * @param E_offload_array electric field offload array
 * @param plasma_offload_array plasma offload array
 * @param neutral_offload_array neutrals offload array
 * @param wall_offload_array wall offload array
 * @param wall_int_offload_array wall integer offload array
 * @param boozer_offload_array boozer offload array
 * @param mhd_offload_array MHD data offload array
 * @param asigma_offload_array atomic data offload array
 * @param offload_array pointer to common offload array created here
 * @param int_offload_array pointer to common offload integer array created here
 *
 * @returns zero on success
 */
int pack_offload_array(
    sim_offload_data* sim, offload_package* offload_data,
    real** B_offload_array, real** E_offload_array, real** plasma_offload_array,
    real** neutral_offload_array, real** wall_offload_array,
    int** wall_int_offload_array, real** boozer_offload_array,
    real** mhd_offload_array, real** asigma_offload_array, real** offload_array,
    int** int_offload_array) {

    /* Pack offload data into single array and free individual offload arrays */
    offload_init_offload(offload_data, offload_array, int_offload_array);

    offload_pack(offload_data, offload_array, *B_offload_array,
                 sim->B_offload_data.offload_array_length,
                 int_offload_array, NULL, 0);
    B_field_free_offload(&sim->B_offload_data, B_offload_array);

    offload_pack(offload_data, offload_array, *E_offload_array,
                 sim->E_offload_data.offload_array_length,
                 int_offload_array, NULL, 0);
    E_field_free_offload(&sim->E_offload_data, E_offload_array);

    offload_pack(offload_data, offload_array, *plasma_offload_array,
                 sim->plasma_offload_data.offload_array_length,
                 int_offload_array, NULL, 0);
    plasma_free_offload(&sim->plasma_offload_data, plasma_offload_array);

    offload_pack(offload_data, offload_array, *neutral_offload_array,
                 sim->neutral_offload_data.offload_array_length,
                 int_offload_array, NULL, 0);
    neutral_free_offload(&sim->neutral_offload_data, neutral_offload_array);

    offload_pack(offload_data, offload_array, *wall_offload_array,
                 sim->wall_offload_data.offload_array_length,
                 int_offload_array, *wall_int_offload_array,
                 sim->wall_offload_data.int_offload_array_length);
    wall_free_offload(&sim->wall_offload_data, wall_offload_array,
                      wall_int_offload_array);

    offload_pack(offload_data, offload_array, *boozer_offload_array,
                 sim->boozer_offload_data.offload_array_length,
                 int_offload_array, NULL, 0);
    boozer_free_offload(&sim->boozer_offload_data, boozer_offload_array);

    offload_pack(offload_data, offload_array, *mhd_offload_array,
                 sim->mhd_offload_data.offload_array_length,
                 int_offload_array, NULL, 0);
    mhd_free_offload(&sim->mhd_offload_data, mhd_offload_array);

    offload_pack(offload_data, offload_array, *asigma_offload_array,
                 sim->asigma_offload_data.offload_array_length,
                 int_offload_array, NULL, 0);
    asigma_free_offload(&sim->asigma_offload_data, asigma_offload_array);

    return 0;
}


/**
 * @brief Create and store run group and marker inistate.
 *
 * @param sim simulation offload data struct
 * @param ps marker state array for this process
 * @param n_tot total number of markers in this simulation
 * @param qid unique identifier for this run group
 *
 * @returns zero on success
 */
int write_rungroup(sim_offload_data* sim, particle_state* ps, int n_tot,
                   char* qid) {

    if(sim->mpi_rank == sim->mpi_root) {
        /* Initialize results group in the output file */
        print_out0(VERBOSE_IO, sim->mpi_rank, sim->mpi_root,
                   "\nPreparing output.\n");
        if( hdf5_interface_init_results(sim, qid, "run") ) {
            print_out0(VERBOSE_MINIMAL, sim->mpi_rank, sim->mpi_root,
                       "\nInitializing output failed.\n"
                       "See stderr for details.\n");
            /* Free offload data and terminate */
            return 1;
        }
        strcpy(sim->qid, qid);
    }

    /* Gather particle states so that we can write inistate */
    int n_gather;
    particle_state* ps_gather;
    mpi_gather_particlestate(ps, &ps_gather, &n_gather, n_tot,
                             sim->mpi_rank, sim->mpi_size, sim->mpi_root);

    if(sim->mpi_rank == sim->mpi_root) {
        /* Write inistate */
        if(hdf5_interface_write_state(
            sim->hdf5_out, "inistate", n_gather, ps_gather)) {
            print_out0(VERBOSE_MINIMAL, sim->mpi_rank, sim->mpi_root,
                       "\n"
                       "Writing inistate failed.\n"
                       "See stderr for details.\n"
                       "\n");
            /* Free offload data and terminate */
            return 1;
        }
        print_out0(VERBOSE_NORMAL, sim->mpi_rank, sim->mpi_root,
                   "\nInistate written.\n");
    }
    free(ps_gather);

    return 0;
}


/**
 * @brief Offload data to target, carry out the simulation, and return to host.
 *
 * @param sim simulation offload data struct
 * @param n_tot total number of markers
 * @param n_proc number of markers in this process
 * @param pin marker state array for this process (deallocated here)
 * @param offload_data packed offload data struct
 * @param offload_array packed offload array containing the input data
 * @param int_offload_array packed offload integer array containg the input data
 * @param n_gather pointer for storing the number of markers in pout (either
 *        n_tot or n_proc)
 * @param pout pointer to array containing all endstates in the simulation
 * @param diag_offload_array array to store output data
 *
 * @return zero on success
 */
int offload_and_simulate(
    sim_offload_data* sim, int n_tot, int n_proc, particle_state* pin,
    offload_package* offload_data, real* offload_array, int* int_offload_array,
    int* n_gather, particle_state** pout, real* diag_offload_array) {

    /* Divide markers among host and target */
    int n_mic = 0;
    int n_host = nprts;

    double mic_start = 0, mic_end=0, host_start=0, host_end=0;

    /* Empty message buffer before proceeding to offload */
    fflush(stdout);

    /* Allow threads to spawn threads */
    //omp_set_nested(1);

    /* Actual marker simulation happens here. Threads are spawned which
     * distribute the execution between target(s) and host. Both input and
     * diagnostic offload arrays are mapped to target. Simulation is initialized
     * at the target and completed within the simulate() function.*/
    //#pragma omp parallel sections num_threads(3)
    {
        /* Run simulation on first target */

#ifndef TARGET
        /* No target, marker simulation happens where the code execution began.
         * Offloading is only emulated. */
        //#pragma omp section
        {
            host_start = omp_get_wtime();
            simulate(0, n_host, pin+2*n_mic, sim, offload_data,
                offload_array, int_offload_array, diag_offload_array);
            host_end = omp_get_wtime();
        }
    }

    /* Code execution returns to host. */
    MPI_Barrier(MPI_COMM_WORLD);
    print_out0(VERBOSE_NORMAL, mpi_rank, "gpu %lf s, host %lf s\n",
               mic_end-mic_start, host_end-host_start);

    /* Gather output data */
    mpi_gather_particlestate(pin, pout, n_gather, n_tot, sim->mpi_rank,
                             sim->mpi_size, sim->mpi_root);
    free(pin);

    mpi_gather_diag(&sim->diag_offload_data, diag_offload_array, n_tot,
                    sim->mpi_rank, sim->mpi_size, sim->mpi_root);

    return 0;
}


/**
 * @brief Store simulation output data.
 *
 * @param sim simulation offload data
 * @param ps marker endstate array to be written
 * @param n_tot number of markers
 * @param diag_offload_array diagnostics offload data array
 *
 * @return zero on success
 */
int write_output(sim_offload_data* sim, particle_state* ps, int n_tot,
                 real* diag_offload_array){

    if(sim->mpi_rank == sim->mpi_root) {
        /* Write endstate */
        if( hdf5_interface_write_state(
                sim->hdf5_out, "endstate", n_tot, ps)) {
            print_out0(VERBOSE_MINIMAL, sim->mpi_rank, sim->mpi_root,
                   "\nWriting endstate failed.\n"
                   "See stderr for details.\n");
            return 1;
        }
        print_out0(VERBOSE_NORMAL, sim->mpi_rank, sim->mpi_root,
                   "Endstate written.\n");
    }

    /* Combine diagnostic data and write it to HDF5 file */
    print_out0(VERBOSE_MINIMAL, sim->mpi_rank, sim->mpi_root,
                   "\nCombining and writing diagnostics.\n");
    int err_writediag = 0;

    if(sim->mpi_rank == sim->mpi_root) {
        err_writediag = hdf5_interface_write_diagnostics(sim,
            diag_offload_array, sim->hdf5_out);
    }
    if(err_writediag) {
        print_out0(VERBOSE_MINIMAL, sim->mpi_rank, sim->mpi_root,
                   "\nWriting diagnostics failed.\n"
                   "See stderr for details.\n");
        return 1;
    }
    else {
        print_out0(VERBOSE_MINIMAL, sim->mpi_rank, sim->mpi_root,
                   "Diagnostics written.\n");
    }

    return 0;
}


/**
 * @brief Read command line arguments and modify sim struct accordingly
 *
 * The command line arguments are in, out, mpi_size, and mpi_rank which
 * correspond to input file, output file, number of MPI processes and
 * the rank of this process. These are stored in simulation offload data struct
 * as (default values, used if the specific argument was not given, are in
 * parenthesis):
 *
 * - sim->hdf5_in     = "in.h5" ("ascot.h5")
 * - sim->hdf5_out    = "out" (sim->hdf5_in is copied here)
 * - sim->mpi_rank    = 0
 * - sim->mpi_size    = 0
 * - sim->desc        = "No description"
 *
 * If the arguments could not be parsed, this function returns a non-zero exit
 * value.
 *
 * @param argc argument count as given to main()
 * @param argv argument vector as given to main()
 * @param sim pointer to offload data struct
 *
 * @return Zero if success
 */
int read_arguments(int argc, char** argv, sim_offload_data* sim) {
    struct option longopts[] = {
        {"in", required_argument, 0, 1},
        {"out", required_argument, 0, 2},
        {"mpi_size", required_argument, 0, 3},
        {"mpi_rank", required_argument, 0, 4},
        {"d", required_argument, 0, 5},
        {"options", required_argument, 0, 6},
        {"bfield",  required_argument, 0, 7},
        {"efield",  required_argument, 0, 8},
        {"marker",  required_argument, 0, 9},
        {"wall",    required_argument, 0, 10},
        {"plasma",  required_argument, 0, 11},
        {"neutral", required_argument, 0, 12},
        {"boozer",  required_argument, 0, 13},
        {"mhd",     required_argument, 0, 14},
        {"asigma",  required_argument, 0, 15},
        {0, 0, 0, 0}
    };

    /* Initialize default values */
    sim->hdf5_in[0]     = '\0';
    sim->hdf5_out[0]    = '\0';
    sim->mpi_rank       = 0;
    sim->mpi_size       = 0;
    strcpy(sim->description, "No description.");
    sim->qid_options[0] = '\0';
    sim->qid_bfield[0]  = '\0';
    sim->qid_efield[0]  = '\0';
    sim->qid_marker[0]  = '\0';
    sim->qid_wall[0]    = '\0';
    sim->qid_plasma[0]  = '\0';
    sim->qid_neutral[0] = '\0';
    sim->qid_boozer[0]  = '\0';
    sim->qid_mhd[0]     = '\0';
    sim->qid_asigma[0]  = '\0';
    sim->qid_nbi[0]     = '\0';

    /* Read user input */
    int c;
    int slen;  /* String length */
    while((c = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
        switch(c) {
            case 1:
                /* The .hdf5 filename can be specified with or without the
                 * trailing .h5 */
                slen = strlen(optarg);
                if ( slen > 3 && !strcmp(optarg+slen-3, ".h5") ) {
                    strcpy(sim->hdf5_in, optarg);
                    (sim->hdf5_in)[slen-3]='\0';
                }
                else {
                    strcpy(sim->hdf5_in, optarg);
                }
                break;
            case 2:
                /* The .hdf5 filename can be specified with or without the
                 * trailing .h5 */
                slen = strlen(optarg);
                if ( slen > 3 && !strcmp(optarg+slen-3, ".h5") ) {
                    strcpy(sim->hdf5_out, optarg);
                    (sim->hdf5_out)[slen-3]='\0';
                }
                else {
                    strcpy(sim->hdf5_out, optarg);
                }
                break;
            case 3:
                sim->mpi_size = atoi(optarg);
                break;
            case 4:
                sim->mpi_rank = atoi(optarg);
                break;
            case 5:
                strcpy(sim->description, optarg);
                break;
            case 6:
                strcpy(sim->qid_options, optarg);
                break;
            case 7:
                strcpy(sim->qid_bfield, optarg);
                break;
            case 8:
                strcpy(sim->qid_efield, optarg);
                break;
            case 9:
                strcpy(sim->qid_marker, optarg);
                break;
            case 10:
                strcpy(sim->qid_wall, optarg);
                break;
            case 11:
                strcpy(sim->qid_plasma, optarg);
                break;
            case 12:
                strcpy(sim->qid_neutral, optarg);
                break;
            case 13:
                strcpy(sim->qid_boozer, optarg);
                break;
            case 14:
                strcpy(sim->qid_mhd, optarg);
                break;
            case 15:
                strcpy(sim->qid_asigma, optarg);
                break;
            default:
                // Unregonizable argument(s). Tell user how to run ascot5_main
                print_out(VERBOSE_MINIMAL,
                          "\nUnrecognized argument. The valid arguments are:\n");
                print_out(VERBOSE_MINIMAL,
                          "--in input file (default: ascot.h5)\n");
                print_out(VERBOSE_MINIMAL,
                          "--out output file (default: same as input)\n");
                print_out(VERBOSE_MINIMAL,
                          "--mpi_size number of independent processes\n");
                print_out(VERBOSE_MINIMAL,
                          "--mpi_rank rank of independent process\n");
                print_out(VERBOSE_MINIMAL,
                          "--d run description maximum of 250 characters\n");
                return 1;
        }
    }

    /* Default value for input file is ascot.h5, and for output same as input
     * file. Adujust hdf5_in and hdf5_out accordingly. For output file, we don't
     * add the .h5 extension here. */
    if(sim->hdf5_in[0] == '\0' && sim->hdf5_out[0] == '\0') {
        /* No input, use default values for both */
        strcpy(sim->hdf5_in,  "ascot.h5");
        strcpy(sim->hdf5_out, "ascot.h5");
    }
    else if(sim->hdf5_in[0] == '\0' && sim->hdf5_out[0] != '\0') {
        /* Output file is given but the input file is not */
        strcpy(sim->hdf5_in,  "ascot.h5");
        strcat(sim->hdf5_out, ".h5");
    }
    else if(sim->hdf5_in[0] != '\0' && sim->hdf5_out[0] == '\0') {
        /* Input file is given but the output is not */
        strcat(sim->hdf5_in, ".h5");
        strcpy(sim->hdf5_out, sim->hdf5_in);
    }
    else {
        /* Both input and output files are given */
        strcat(sim->hdf5_in,  ".h5");
        strcat(sim->hdf5_out, ".h5");
    }
    return 0;
}

/**
 * @brief Writes a summary of what happened to the markers during simulation
 *
 * This function writes a summary of marker end conditions and possible
 * simulation-time errors. Since simulation can have billions and billions of
 * markers, we only show how many markers had specific error or end condition.
 *
 * End conditions and errors are printed in human-readable format.
 *
 * This function is called by the root MPI process only.
 *
 * @param ps array of marker states after simulation has finished
 * @param n_tot number of markers in the array
 */
void print_marker_summary(particle_state* ps, int n_tot) {

    print_out(VERBOSE_MINIMAL, "\nSummary of results:\n");

    /* Temporary arrays that are needed to store unique end conditions and
     * errors. We can have at most n_tot different values. */
    int* temp = (int*)malloc(n_tot*sizeof(int));
    int* unique = (int*)malloc(n_tot*sizeof(int));
    int* count = (int*)malloc((n_tot+1)*sizeof(int));

    /* First we find out what the end conditions are */
    for(int i=0; i<n_tot; i++) {
        temp[i] = ps[i].endcond;
    }
    math_uniquecount(temp, unique, count, n_tot);
    count[n_tot] = 0;/* This ensures the following while loops are terminated.*/

    /* Print the end conditions, if marker has multiple end conditions, separate
     * those with "and". */
    int i = 0;
    while(count[i] > 0) {
        /* Initialize */
        int endconds[32];
        for(int j=0; j<32;j++) {
            endconds[j] = 0;
        }
        char endcondstr[256];
        endcondstr[0] = '\0';

        /* Represent all end conditions with a single string and print it */
        int j = 0;
        endcond_parse(unique[i], endconds);
        while(endconds[j]) {
            if(j>0) {
                strcat(endcondstr, " and ");
            }
            char temp[256];
            endcond_parse2str(endconds[j], temp);
            strcat(endcondstr, temp);
            j++;
        }
        if(j == 0) {
            sprintf(endcondstr, "Aborted");
        }
        print_out(VERBOSE_MINIMAL, "%9d markers had end condition %s\n",
                  count[i], endcondstr);
        i++;
    }

    /* Empty line between end conditions and errors */
    print_out(VERBOSE_MINIMAL, "\n");

    /* Find all errors */
    for(int i=0; i<n_tot; i++) {
        temp[i] = (int)(ps[i].err);
    }
    math_uniquecount(temp, unique, count, n_tot);

    /* Go through each unique error and represent is a string */
    i = 0;
    while(count[i] > 0) {
        if(unique[i] == 0) {
            /* Value 0 indicates no error occurred so skip that */
            i++;
            continue;
        }
        char msg[256];
        char line[256];
        char file[256];
        error_parse2str(unique[i], msg, line, file);
        print_out(VERBOSE_MINIMAL,
                  "%9d markers were aborted with an error message:\n"
                  "          %s\n"
                  "          at line %s in %s\n",
                  count[i], msg, line, file);
        i++;
    }

    /* If count[0] equals to number of markers and their error field is zero,
     * we have no markers that were aborted. */
    if(count[0] == n_tot && unique[0] == 0) {
        print_out(VERBOSE_MINIMAL,
                  "          No markers were aborted.\n");
    }

    /* Free temporary arrays */
    free(temp);
    free(unique);
    free(count);
}
