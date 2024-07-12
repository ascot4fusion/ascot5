/**
 * @file bbnbi5.c
 * @brief Interface to BBNBI5
 *
 * BBNBI5 models neutral beam injectors and is used to evaluate shine-through
 * and beam birth-profile. Neutral markers are generated from injector geometry
 * and traced until they ionize or hit the wall. Several injectors can be
 * modelled simultaneously keeping in mind that in this case the output
 * the injector from which a particle originated is lost.
 */
#include <getopt.h>
#include <math.h>
#ifdef MPI
  #include <mpi.h>
#endif
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "ascot5.h"
#include "consts.h"
#include "gitver.h"
#include "math.h"
#include "physlib.h"
#include "print.h"
#include "simulate.h"
#include "endcond.h"
#include "random.h"
#include "particle.h"
#include "suzuki.h"
#include "B_field.h"
#include "plasma.h"
#include "neutral.h"
#include "wall.h"
#include "asigma.h"
#include "nbi.h"
#include "diag.h"
#include "bbnbi5.h"

int bbnbi_read_arguments(int argc, char** argv, sim_offload_data* sim,
                         int* nprt, real* t1, real* t2);
void bbnbi_trace_markers(particle_queue *pq, sim_data* sim);
void bbnbi_inject_markers(particle_state* p, int nprt, int ngenerated, real t0,
                          real t1, nbi_injector* inj, sim_data* sim);

/**
 * @brief Simulate NBI injection
 *
 * This function initializes neutrals and traces them until they have ionized or
 * hit the wall.
 *
 * @param sim pointer to the simulation offload data structure
 * @param nprt number of markers to be injected
 * @param t1 time instant when the injector is turned on
 * @param t2 time instant when the injector is turned off
 * @param B_offload_array pointer to the magnetic field data
 * @param plasma_offload_array pointer to the plasma data
 * @param neutral_offload_array pointer to the neutral data
 * @param wall_offload_array pointer to the wall data
 * @param wall_int_offload_array pointer to the wall int data
 * @param asigma_offload_array pointer to the atomic sigma data
 * @param nbi_offload_array pointer to the nbi data
 * @param p pointer to the marker array which is allocated here
 * @param diag_offload_array pointer to the diagnostics data
 */
void bbnbi_simulate(
    sim_offload_data* sim, int nprt, real t1, real t2, real* B_offload_array,
    real* plasma_offload_array, real* neutral_offload_array,
    real* wall_offload_array, int* wall_int_offload_array,
    real* asigma_offload_array, real* nbi_offload_array, particle_state** p,
    real* diag_offload_array) {

    /* Initialize input data */
    sim_data sim_data;
    sim_init(&sim_data, sim);
    random_init(&sim_data.random_data, time(NULL));
    B_field_init(&sim_data.B_data, &sim->B_offload_data, B_offload_array);
    plasma_init(&sim_data.plasma_data, &sim->plasma_offload_data,
                plasma_offload_array);
    neutral_init(&sim_data.neutral_data, &sim->neutral_offload_data,
                 neutral_offload_array);
    wall_init(&sim_data.wall_data, &sim->wall_offload_data, wall_offload_array,
              wall_int_offload_array);
    asigma_init(&sim_data.asigma_data, &sim->asigma_offload_data,
                asigma_offload_array);
    nbi_init(&sim_data.nbi_data, &sim->nbi_offload_data, nbi_offload_array);
    diag_init(&sim_data.diag_data, &sim->diag_offload_data, diag_offload_array);

    /* Calculate total NBI power so that we can distribute markers along
     * the injectors according to their power */
    real total_power = 0;
    for(int i=0; i < sim_data.nbi_data.ninj; i++) {
        total_power += sim_data.nbi_data.inj[i].power;
    }

    /* Initialize particle struct */
    *p = (particle_state*) malloc(nprt * sizeof(particle_state));

    /* Generate markers at the injectors */
    int nprt_generated = 0;
    for(int i = 0; i < sim_data.nbi_data.ninj; i++) {

        /* Number of markers generated is proportional to NBI power */
        int nprt_inj = ( sim_data.nbi_data.inj[i].power / total_power ) * nprt;
        if(i == sim_data.nbi_data.ninj-1) {
            /* All "remaining" markers goes to the last injector to avoid any
             * rounding issues */
            nprt_inj = nprt - nprt_generated;
        }

        /* Generates markers at the injector location and traces them until
         * they enter the region with magnetic field data */
        bbnbi_inject_markers(&((*p)[nprt_generated]), nprt_inj, nprt_generated,
                             t1, t2, &(sim_data.nbi_data.inj[i]), &sim_data);

        nprt_generated += nprt_inj;
        print_out0(VERBOSE_NORMAL, sim->mpi_rank, sim->mpi_root,
                   "Generated %d markers for injector %d.\n", nprt_inj, i+1);
    }

    /* Place markers in a queue */
    particle_queue pq;
    pq.n = 0;
    for(int i = 0; i < nprt; i++) {
        pq.n++;
    }
    pq.p = (particle_state**) malloc(pq.n * sizeof(particle_state*));
    pq.finished = 0;

    pq.next = 0;
    for(int i = 0; i < nprt; i++) {
        pq.p[pq.next++] = &((*p)[i]);

    }
    pq.next = 0;

    /* Trace neutrals until they are ionized or lost to the wall */
    #pragma omp parallel
    bbnbi_trace_markers(&pq, &sim_data);
}

/**
 * @brief Inject neutrals from an injector
 *
 * This function initializes neutral markers at the beamlet positions and
 * launches them in a (random) direction based on the injector specs.
 * The marker is traced until it enters the magnetic field, at which point
 * the particle struct is filled with only the particle data, and the struct
 * is returned.
 *
 * @param p pointer where generated markers are stored
 * @param nprt number of markers to be injected or generated
 * @param ngenerated number of markers that have already been generated
 * @param t0 time when the injector is turned on
 * @param t1 time when the injector is turned off
 * @param inj pointer to injector data
 * @param sim pointer to the sim struct with initialized data
 */
void bbnbi_inject_markers(particle_state* p, int nprt, int ngenerated, real t0,
                          real t1, nbi_injector* inj, sim_data* sim) {

    /* Set marker weights assuming a large number is created so that the energy
     * fractions of generated markers are close to the injector values */
    real f  =     1.0 * inj->efrac[0] + (1.0/2) * inj->efrac[1]
            + (1.0/3) * inj->efrac[2];
    real weight = (inj->power / inj->energy ) / ( f * nprt );

    /* Inject markers and trace their ballistic trajectories (without any
     * other physics) until they enter the plasma for the first time.     */
    #pragma omp parallel for
    for(int i = 0; i < nprt; i++) {
        real time = t0 + random_uniform(&sim->random_data) * (t1-t0);

        /* Assign initial phase-space coordinates for this marker */
        real xyz[3], vxyz[3], rpz[3], vhat[3];
        nbi_inject(xyz, vxyz, inj, &sim->random_data);
        math_xyz2rpz(xyz, rpz);
        math_unit(vxyz, vhat);

        /* Advance until the marker enters the magnetic field */
        real psi;
        real ds = 1e-3;
        a5err err = B_field_eval_psi(&psi, rpz[0], rpz[1], rpz[2], time,
                                     &sim->B_data);
        while(err) {
            xyz[0] += ds * vhat[0];
            xyz[1] += ds * vhat[1];
            xyz[2] += ds * vhat[2];
            math_xyz2rpz(xyz, rpz);
            err = B_field_eval_psi(&psi, rpz[0], rpz[1], rpz[2], time,
                                   &sim->B_data);
        }

        real vrpz[3];
        math_vec_xyz2rpz(vxyz, vrpz, rpz[1]);
        real gamma = physlib_gamma_vnorm(math_norm(vrpz));

        /* Fill the particle state with particle coordinates */
        p[i].rprt     = rpz[0];
        p[i].phiprt   = rpz[1];
        p[i].zprt     = rpz[2];
        p[i].p_r      = vrpz[0] * gamma * inj->mass;
        p[i].p_phi    = vrpz[1] * gamma * inj->mass;
        p[i].p_z      = vrpz[2] * gamma * inj->mass;
        p[i].mass     = inj->mass;
        p[i].charge   = 0.0;
        p[i].anum     = inj->anum;
        p[i].znum     = inj->znum;
        p[i].weight   = weight;
        p[i].time     = time;
        p[i].mileage  = 0.0;
        p[i].cputime  = 0.0;
        p[i].id       = ngenerated + i + 1;
        p[i].endcond  = 0;
        p[i].walltile = 0;
        p[i].err      = 0;
    }
}

/**
 * @brief Trace a neutral marker until it has ionized or hit wall
 *
 * This function is for the most part identical to simulate_fo with few
 * exceptions relevant for BBNBI.
 *
 * @param pq pointer to the marker queue containing the initial neutrals
 * @param sim pointer to the simu struct with initialized data
 */
void bbnbi_trace_markers(particle_queue *pq, sim_data* sim) {
    int cycle[NSIMD]  __memalign__;
    real hin[NSIMD]  __memalign__;
    int shinethrough[NSIMD] __memalign__;
    real remaining[NSIMD]  __memalign__;
    real threshold[NSIMD]  __memalign__;
    particle_simd_fo p, p0;

    int n_species       = plasma_get_n_species(&sim->plasma_data);
    const int* pls_anum = plasma_get_species_anum(&sim->plasma_data);
    const int* pls_znum = plasma_get_species_znum(&sim->plasma_data);

    /* Init dummy markers */
    for(int i=0; i< NSIMD; i++) {
        p.id[i] = -1;
        p.running[i] = 0;
        hin[i] = 1e-10;
        threshold[i] = random_uniform(&sim->random_data);
        remaining[i] = 1.0;
        shinethrough[i] = 0;
    }

    /* Initialize running particles */
    int n_running = particle_cycle_fo(pq, &p, &sim->B_data, cycle);
    while(n_running > 0) {

        #pragma omp simd
        for(int i=0; i< NSIMD; i++) {
            /* Store marker states */
            particle_copy_fo(&p, i, &p0, i);

            if(p.running[i]) {
                a5err err = 0;

                /* These are needed later */
                real pnorm = math_normc(p.p_r[i], p.p_phi[i], p.p_z[i]);
                real gamma = physlib_gamma_pnorm(p.mass[i], pnorm);
                real ekin  = physlib_Ekin_pnorm(p.mass[i], pnorm);

                /* Advance ballistic trajectory by converting momentum to
                 * cartesian coordinates */
                real prpz[3] = {p.p_r[i], p.p_phi[i], p.p_z[i]};
                real pxyz[3];
                math_vec_rpz2xyz(prpz, pxyz, p.phi[i]);

                real posrpz[3] = {p.r[i], p.phi[i], p.z[i]};
                real posxyz[3], fposxyz[3];
                math_rpz2xyz(posrpz, posxyz);
                fposxyz[0] = posxyz[0] + pxyz[0] * hin[i] / (gamma * p.mass[i]);
                fposxyz[1] = posxyz[1] + pxyz[1] * hin[i] / (gamma * p.mass[i]);
                fposxyz[2] = posxyz[2] + pxyz[2] * hin[i] / (gamma * p.mass[i]);

                /* Back to cylindrical coordinates (note phi is cumulative) */
                p.r[i] = sqrt(fposxyz[0]*fposxyz[0] + fposxyz[1]*fposxyz[1]);
                p.phi[i] += atan2(
                    posxyz[0] * fposxyz[1] - posxyz[1] * fposxyz[0],
                    posxyz[0] * fposxyz[0] + posxyz[1] * fposxyz[1] );
                p.z[i] = fposxyz[2];

                real cosp = cos(p.phi[i]);
                real sinp = sin(p.phi[i]);
                p.p_r[i]   =  pxyz[0] * cosp + pxyz[1] * sinp;
                p.p_phi[i] = -pxyz[0] * sinp + pxyz[1] * cosp;
                p.p_z[i]   =  pxyz[2];

                real ds = hin[i];
                p.mileage[i] += hin[i];

                /* Update background values at the new position */
                real psi, rho[2], pls_dens[MAX_SPECIES], pls_temp[MAX_SPECIES];
                err = B_field_eval_psi(
                    &psi, p.r[i], p.phi[i], p.z[i], p.time[i], &sim->B_data);
                if(!err) {
                    err = B_field_eval_rho(rho, psi, &sim->B_data);
                }
                if(!err && p.rho[i] <= 1.0 && rho[0] > 1.0) {
                    shinethrough[i] = 1;
                }
                p.rho[i] = rho[0];
                /* Update theta value */
                real axisrz[2];
                if(!err) {
                    B_field_get_axis_rz(axisrz, &sim->B_data, p.phi[i]);
                }
                p.theta[i] = atan2(p.z[i]-axisrz[1], p.r[i]-axisrz[0]);

                if(!err) {
                    err = plasma_eval_densandtemp(
                        pls_dens, pls_temp, rho[0], p.r[i], p.phi[i], p.z[i],
                        p.time[i], &sim->plasma_data);
                }

                /* Calculate ionization rate */
                real rate = 0.0;
                if(!err) {
                    real sigmav;
                    if( asigma_eval_bms(
                            &sigmav, p.znum[i], p.anum[i], ekin, p.mass[i],
                            n_species-1, pls_znum, pls_anum, pls_temp[0],
                            &(pls_dens[1]), &sim->asigma_data) ) {
                        err = 1;
                    }
                    rate = pls_dens[0] * sigmav;
                }
                remaining[i] *= exp(-rate * ds);

                /* Check for end conditions */
                if(!err) {
                    real w_coll = 0;
                    int tile = 0;
                    if(shinethrough[i]) {
                        tile = wall_hit_wall(
                            p0.r[i], p0.phi[i], p0.z[i],
                            p.r[i], p.phi[i], p.z[i], &sim->wall_data, &w_coll);
                    }
                    if(tile > 0) {
                        real w = w_coll;
                        p.time[i] = p0.time[i] + w*(p.time[i] - p0.time[i]);
                        p.r[i]    = p0.r[i]    + w*(p.r[i]    - p0.r[i]);
                        p.phi[i]  = p0.phi[i]  + w*(p.phi[i]  - p0.phi[i]);
                        p.z[i]    = p0.z[i]    + w*(p.z[i]    - p0.z[i]);

                        p.walltile[i] = tile;
                        p.endcond[i] |= endcond_wall;
                        p.running[i] = 0;
                    }
                    if(p.mileage[i] > NBI_MAX_DISTANCE) {
                        p.endcond[i] |= endcond_tlim;
                        p.running[i] = 0;
                    }
                    if(remaining[i] < threshold[i]) {
                        p.charge[i] = 1*CONST_E;
                        p.endcond[i] |= endcond_ioniz;
                        p.running[i] = 0;
                    }
                } else {
                    p.err[i] = err;
                    p.running[i] = 0;
                }
            }
        }

        /* Update markers that just finished */
        particle_simd_fo pdiag;
        #pragma omp simd
        for(int i=0; i< NSIMD; i++) {
            /* Use this as a flag for which markers to update in diagnostics */
            pdiag.running[i] = 0;
            if(!p.running[i] && p.id[i] >= 0) {
                p.time[i] += p.mileage[i];

                /* Reset these for the next marker */
                threshold[i] = random_uniform(&sim->random_data);
                remaining[i] = 1.0;
                shinethrough[i] = 0;

                /* Update the magnetic field at the marker position */
                if(!p.err[i]) {
                    real B_dB[15];
                    B_field_eval_B_dB(B_dB, p.r[i], p.phi[i], p.z[i], p.time[i],
                                      &sim->B_data);
                    p.B_r[i]        = B_dB[0];
                    p.B_r_dr[i]     = B_dB[1];
                    p.B_r_dphi[i]   = B_dB[2];
                    p.B_r_dz[i]     = B_dB[3];

                    p.B_phi[i]      = B_dB[4];
                    p.B_phi_dr[i]   = B_dB[5];
                    p.B_phi_dphi[i] = B_dB[6];
                    p.B_phi_dz[i]   = B_dB[7];

                    p.B_z[i]        = B_dB[8];
                    p.B_z_dr[i]     = B_dB[9];
                    p.B_z_dphi[i]   = B_dB[10];
                    p.B_z_dz[i]     = B_dB[11];
                }
            }
            particle_copy_fo(&p, i, &pdiag, i);
            if(!p.running[i] && p.id[i] >= 0) {
                pdiag.running[i] = 1;
            }

            /* Normalize weight with time and add hin so that we don't divide
             * with zero when updating distributions */
            pdiag.time[i]   += hin[i];
            pdiag.weight[i] /= hin[i];
        }

        /* Update distributions for markers that finished */
        diag_update_fo(&sim->diag_data, &sim->B_data, &pdiag, &p);

        /* Update running particles */
        n_running = particle_cycle_fo(pq, &p, &sim->B_data, cycle);
    }
}

/**
 * @brief Read command line arguments
 *
 * Read in command line arguments, input and output names and mpi parameters
 * are stored in sim structure as with ascot5, number of markers is passed
 * as an argument.
 *
 * If the arguments could not be parsed, this function returns a non-zero exit
 * value.
 *
 * @param argc argument count as given to main()
 * @param argv argument vector as given to main()
 * @param sim pointer to offload data struct
 * @param nprt pointer to integer where number of markers is stored
 * @param t1 pointer to store beginning of time interval
 * @param t2 pointer to store end of the time interval
 *
 * @return Zero if success
 */
int bbnbi_read_arguments(int argc, char** argv, sim_offload_data* sim,
                         int* nprt, real* t1, real* t2) {
    struct option longopts[] = {
        {"in",       required_argument, 0,  1},
        {"out",      required_argument, 0,  2},
        {"mpi_size", required_argument, 0,  3},
        {"mpi_rank", required_argument, 0,  4},
        {"d",        required_argument, 0,  5},
        {"bfield",   required_argument, 0,  6},
        {"wall",     required_argument, 0,  7},
        {"plasma",   required_argument, 0,  8},
        {"nbi",      required_argument, 0,  9},
        {"n",        required_argument, 0, 10},
        {"t1",       required_argument, 0, 11},
        {"t2",       required_argument, 0, 12},
        {0, 0, 0, 0}
    };

    /* Initialize default values */
    sim->hdf5_in[0]     = '\0';
    sim->hdf5_out[0]    = '\0';
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
    sim->mpi_rank       = 0;
    sim->mpi_size       = 0;
    *nprt               = 10000;
    *t1                 = 0.0;
    *t2                 = 0.0;
    strcpy(sim->description, "TAG");

    /* Read user input */
    int c;
    int slen;
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
                strcpy(sim->qid_bfield, optarg);
                break;
            case 7:
                strcpy(sim->qid_wall, optarg);
                break;
            case 8:
                strcpy(sim->qid_plasma, optarg);
                break;
            case 9:
                strcpy(sim->qid_nbi, optarg);
                break;
            case 10:
                *nprt = atoi(optarg);
                break;
            case 11:
                *t1 = atof(optarg);
                break;
            case 12:
                *t2 = atof(optarg);
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
                print_out(VERBOSE_MINIMAL,
                          "--n number of markers to generate (default: 10000)\n");
                print_out(VERBOSE_MINIMAL,
                          "--t1 time when injectors are turned on (default: 0.0 s)\n");
                print_out(VERBOSE_MINIMAL,
                          "--t2 time when injectors are turned off (default: 0.0 s)\n");
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
