/**
 * @file simulate_fo_fixed.c
 * @brief Simulate particles using fixed time-step
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#ifndef __NVCOMPILER
#include <immintrin.h>
#endif
#include <math.h>
#include "../ascot5.h"
#include "../physlib.h"
#include "../simulate.h"
#include "../particle.h"
#include "../wall.h"
#include "../diag.h"
#include "../B_field.h"
#include "../E_field.h"
#include "../plasma.h"
#include "../endcond.h"
#include "../math.h"
#include "../consts.h"
#include "simulate_fo_fixed.h"
#include "step/step_fo_vpa.h"
#include "mccc/mccc.h"
#include "atomic.h"

#pragma omp declare simd uniform(sim)
real simulate_fo_fixed_inidt(sim_data* sim, particle_simd_fo* p, int i);

void simulate_fo_fixed_copy_to_gpu(sim_data* sim, particle_simd_fo *p_ptr, particle_simd_fo *p0_ptr, B_field_data* Bdata, E_field_data* Edata, particle_loc*  p_loc, real* hin, real* rnd);

void simulate_fo_fixed_copy_from_gpu(sim_data* sim, particle_simd_fo *p_ptr);

/**
 * @brief Simulates particles using fixed time-step
 *
 * The simulation includes:
 * - orbit-following with Volume-Preserving Algorithm
 * - Coulomb collisions with Euler-Maruyama method
 *
 * The simulation is carried until all markers have met some
 * end condition or are aborted/rejected. The final state of the
 * markers is stored in the given marker array. Other output
 * is stored in the diagnostic array.
 *
 * The time-step is user-defined: either a directly given fixed value
 * or a given fraction of gyrotime.
 *
 * @param pq particles to be simulated
 * @param sim simulation data struct
 */
void simulate_fo_fixed(particle_queue* pq, sim_data* sim) {
    int cycle[NSIMD]  __memalign__; // Indicates if a new marker was initialized
    real hin[NSIMD]  __memalign__;  // Time step

    real cputime, cputime_last; // Global cpu time: recent and previous record

    particle_simd_fo p;  // This array holds current states
    particle_simd_fo p0; // This array stores previous states

    /* Init dummy markers */
    for(int i=0; i< NSIMD; i++) {
        p.id[i] = -1;
        p.running[i] = 0;
    }

    /* Initialize running particles */
    int n_running = particle_cycle_fo(pq, &p, &sim->B_data, cycle);

    /* Determine simulation time-step */
    #pragma omp simd
    for(int i = 0; i < NSIMD; i++) {
        if(cycle[i] > 0) {
            hin[i] = simulate_fo_fixed_inidt(sim, &p, i);

        }
    }

    cputime_last = A5_WTIME;

    /* MAIN SIMULATION LOOP
     * - Store current state
     * - Integrate motion due to background EM-field (orbit-following)
     * - Integrate scattering due to Coulomb collisions
     * - Advance time
     * - Check for end condition(s)
     * - Update diagnostics
     */
    int diag_data_field_size = sim->diag_data.diagorb.Nmrk*sim->diag_data.diagorb.Npnt;
    particle_simd_fo *p_ptr=&p;
    particle_simd_fo *p0_ptr=&p0;
    B_field_data* Bdata = &sim->B_data;
    E_field_data* Edata = &sim->E_data;
    particle_loc  p_loc;
    real rnd[3*NSIMD];

#ifdef GPU
    simulate_fo_fixed_copy_to_gpu(sim, p_ptr, p0_ptr, Bdata, Edata, &p_loc, hin, rnd);
#endif    
    while(n_running > 0) {
        /* Store marker states */
        //#pragma omp simd
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for(int i = 0; i < NSIMD; i++) {
	  particle_copy_fo(p_ptr, i, p0_ptr, i);
        }
        /*************************** Physics **********************************/

        /* Set time-step negative if tracing backwards in time */
        //#pragma omp simd
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for(int i = 0; i < NSIMD; i++) {
            if(sim->reverse_time) {
                hin[i]  = -hin[i];
            }
        }

        /* Volume preserving algorithm for orbit-following */
        if(sim->enable_orbfol) {
            if(sim->enable_mhd) {
#ifdef GPU
	      printf("step_fo_vpa_mhd NOT YET PORTED TO GPU");
	      exit(1);
#endif
                step_fo_vpa_mhd(&p, hin, &sim->B_data, &sim->E_data,
                                &sim->boozer_data, &sim->mhd_data);
            }
            else {
	      step_fo_vpa(p_ptr, hin, &sim->B_data, &sim->E_data);
            }
        }

        /* Switch sign of the time-step again if it was reverted earlier */
        //#pragma omp simd
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for(int i = 0; i < NSIMD; i++) {
            if(sim->reverse_time) {
                hin[i]  = -hin[i];
            }
        }

        /* Euler-Maruyama for Coulomb collisions */
        if(sim->enable_clmbcol) {
#if !defined(GPU) || defined(RANDOM_LCG)
	  mccc_fo_euler(p_ptr, hin, &sim->plasma_data,
#if defined(RANDOM_LCG)
			&sim->random_data,
#else
			sim->random_data,
#endif			
			&sim->mccc_data,
			rnd);
#else
	  printf("mccc_fo_euler ported on GPU only for RANDOM_LCG");
	  exit(1);	  
#endif	  
        }
        /* Atomic reactions */
        if(sim->enable_atomic) {
#ifdef GPU
	  printf("atomic_fo NOT YET PORTED TO GPU");
	  exit(1);
#endif
            atomic_fo(p_ptr, hin, &sim->plasma_data, &sim->neutral_data,
                      &sim->random_data, &sim->asigma_data,
                      &sim->enable_atomic);
        }
        /**********************************************************************/


        /* Update simulation and cpu times */
        cputime = A5_WTIME;
        //#pragma omp simd
        GPU_PARALLEL_LOOP_ALL_LEVELS
        for(int i = 0; i < NSIMD; i++) {
            if(p.running[i]){
                p.time[i]    += ( 1.0 - 2.0 * ( sim->reverse_time > 0 ) ) * hin[i];
                p.mileage[i] += hin[i];
                p.cputime[i] += cputime - cputime_last;
            }
        }
        cputime_last = cputime;

        /* Check possible end conditions */
        endcond_check_fo(p_ptr, p0_ptr, sim);

        /* Update diagnostics */
        if(!(sim->record_mode)) {
            /* Record particle coordinates */
	  diag_update_fo(&sim->diag_data, &sim->B_data, p_ptr, p0_ptr, &p_loc);
        }
        else {
#ifdef GPU
	  printf("particle_fo_to_gc NOT YET PORTED TO GPU");
	  exit(1);
#endif	  
	  /* Instead of particle coordinates we record guiding center */

            // Dummy guiding centers
            particle_simd_gc gc_f;
            particle_simd_gc gc_i;

            /* Particle to guiding center transformation */
            #pragma omp simd
            for(int i=0; i<NSIMD; i++) {
                if(p.running[i]) {
                    particle_fo_to_gc(&p, i, &gc_f, &sim->B_data);
                }
                else {
                    gc_f.id[i] = p.id[i];
                    gc_f.running[i] = 0;
                }
                if(p0.running[i]) {
                    particle_fo_to_gc(&p0, i, &gc_i, &sim->B_data);
                }
                else {
                    gc_i.id[i] = p0.id[i];
                    gc_i.running[i] = 0;
                }
            }
            diag_update_gc(&sim->diag_data, &sim->B_data, &gc_f, &gc_i);
        }

        /* Update running particles */
#ifdef GPU
	n_running = 0;
	GPU_PARALLEL_LOOP_ALL_LEVELS_REDUCTION(n_running)
	for(int i = 0; i < NSIMD; i++)
	  {
	    if(p_ptr->running[i] > 0) n_running++;
	  }
#else
	n_running = particle_cycle_fo(pq, &p, &sim->B_data, cycle);
#endif
#ifndef GPU	
        /* Determine simulation time-step for new particles */
        //#pragma omp simd
	GPU_PARALLEL_LOOP_ALL_LEVELS
        for(int i = 0; i < NSIMD; i++) {
	  if(cycle[i] > 0)
	    {
	      hin[i] = simulate_fo_fixed_inidt(sim, &p, i);
	    }
        }
#endif	
    }
    /* All markers simulated! */

#ifdef GPU
    simulate_fo_fixed_copy_from_gpu(sim, p_ptr);
    n_running = particle_cycle_fo(pq, &p, &sim->B_data, cycle);
#endif    
}

/**
 * @brief Calculates time step value
 *
 * The time step is calculated as a user-defined fraction of gyro time,
 * whose formula accounts for relativity, or an user defined value
 * is used as is depending on simulation options.
 *
 * @param sim pointer to simulation data struct
 * @param p SIMD array of markers
 * @param i index of marker for which time step is assessed
 *
 * @return Calculated time step
 */
real simulate_fo_fixed_inidt(sim_data* sim, particle_simd_fo* p, int i) {

    real h;

    /* Value defined directly by user */
    if(sim->fix_usrdef_use) {
        h = sim->fix_usrdef_val;
    }
    else {
      /* Value calculated from gyrotime */
        real Bnorm = math_normc( p->B_r[i], p->B_phi[i], p->B_z[i] );
        real pnorm = math_normc( p->p_r[i], p->p_phi[i], p->p_z[i] );
        real gyrotime = CONST_2PI/
            phys_gyrofreq_pnorm(p->mass[i], p->charge[i], pnorm, Bnorm);
        h = gyrotime/sim->fix_gyrodef_nstep;
    }

    return h;
}


void simulate_fo_fixed_copy_to_gpu(sim_data* sim, particle_simd_fo *p_ptr, particle_simd_fo *p0_ptr, B_field_data* Bdata, E_field_data* Edata, particle_loc*  p_loc, real* hin, real* rnd) {

  GPU_MAP_TO_DEVICE(
		      p_loc[0:1],\
		      p_loc->r_arr1[0:NSIMD],\
		      p_loc->r_arr2[0:NSIMD],\
		      p_loc->r_arr3[0:NSIMD],\
		      p_loc->r_arr4[0:NSIMD],\
		      p_loc->r_arr5[0:NSIMD],\
		      p_loc->i_arr1[0:NSIMD],\
		      p_loc->i_arr2[0:NSIMD],\
		      p_loc->i_arr3[0:NSIMD],\
		      p_loc->i_arr4[0:NSIMD],\
		      p_loc->i_arr5[0:NSIMD],\
		      p_loc->i_arr6[0:NSIMD],\
		      p_loc->i_arr7[0:NSIMD],\
		      p_loc->i_arr8[0:NSIMD],\
		      p_loc->i_arr9[0:NSIMD],\
  		      p_ptr[0:1],\
		      p_ptr->running        [0:NSIMD],\
		      p_ptr->r              [0:NSIMD],\
		      p_ptr->phi            [0:NSIMD],\
		      p_ptr->p_r            [0:NSIMD],\
		      p_ptr->p_phi          [0:NSIMD],\
		      p_ptr->p_z            [0:NSIMD],\
		      p_ptr->mileage        [0:NSIMD],\
		      p_ptr->z              [0:NSIMD],\
		      p_ptr->charge         [0:NSIMD],\
		      p_ptr->mass           [0:NSIMD],\
		      p_ptr->B_r            [0:NSIMD],\
		      p_ptr->B_r_dr         [0:NSIMD],\
		      p_ptr->B_r_dphi       [0:NSIMD],\
		      p_ptr->B_r_dz         [0:NSIMD],\
		      p_ptr->B_phi          [0:NSIMD],\
		      p_ptr->B_phi_dr       [0:NSIMD],\
		      p_ptr->B_phi_dphi     [0:NSIMD],\
		      p_ptr->B_phi_dz       [0:NSIMD],\
		      p_ptr->B_z            [0:NSIMD],\
		      p_ptr->B_z_dr         [0:NSIMD],\
		      p_ptr->B_z_dphi       [0:NSIMD],\
		      p_ptr->B_z_dz         [0:NSIMD],\
		      p_ptr->rho            [0:NSIMD],\
		      p_ptr->theta          [0:NSIMD],\
		      p_ptr->err            [0:NSIMD],\
		      p_ptr->time           [0:NSIMD],\
		      p_ptr->weight         [0:NSIMD],\
		      p_ptr->cputime        [0:NSIMD],\
		      p_ptr->id             [0:NSIMD],\
		      p_ptr->endcond        [0:NSIMD],\
		      p_ptr->walltile       [0:NSIMD],\
		      p_ptr->index          [0:NSIMD],\
		      p_ptr->znum           [0:NSIMD],\
		      p_ptr->anum           [0:NSIMD],\
		      p_ptr->bounces        [0:NSIMD],\
  		      p0_ptr[0:1],\
		      p0_ptr->running       [0:NSIMD],\
		      p0_ptr->r             [0:NSIMD],\
		      p0_ptr->phi           [0:NSIMD],\
		      p0_ptr->p_r           [0:NSIMD],\
		      p0_ptr->p_phi         [0:NSIMD],\
		      p0_ptr->p_z           [0:NSIMD],\
		      p0_ptr->mileage       [0:NSIMD],\
		      p0_ptr->z             [0:NSIMD],\
		      p0_ptr->charge        [0:NSIMD],\
		      p0_ptr->mass          [0:NSIMD],\
		      p0_ptr->B_r           [0:NSIMD],\
		      p0_ptr->B_r_dr        [0:NSIMD],\
		      p0_ptr->B_r_dphi      [0:NSIMD],\
		      p0_ptr->B_r_dz        [0:NSIMD],\
		      p0_ptr->B_phi         [0:NSIMD],\
		      p0_ptr->B_phi_dr      [0:NSIMD],\
		      p0_ptr->B_phi_dphi    [0:NSIMD],\
		      p0_ptr->B_phi_dz      [0:NSIMD],\
		      p0_ptr->B_z           [0:NSIMD],\
		      p0_ptr->B_z_dr        [0:NSIMD],\
		      p0_ptr->B_z_dphi      [0:NSIMD],\
		      p0_ptr->B_z_dz        [0:NSIMD],\
		      p0_ptr->rho           [0:NSIMD],\
		      p0_ptr->theta         [0:NSIMD],\
		      p0_ptr->err           [0:NSIMD],\
		      p0_ptr->time          [0:NSIMD],\
		      p0_ptr->weight        [0:NSIMD],\
		      p0_ptr->cputime       [0:NSIMD],\
		      p0_ptr->id            [0:NSIMD],\
		      p0_ptr->endcond       [0:NSIMD],\
		      p0_ptr->walltile      [0:NSIMD],\
		      p0_ptr->index         [0:NSIMD],\
		      p0_ptr->znum          [0:NSIMD],\
		      p0_ptr->anum          [0:NSIMD],\
		      p0_ptr->bounces       [0:NSIMD],\
		      hin[0:NSIMD],\
       		      sim[0:1],		\
		      sim->diag_data.dist5D.histogram[0:sim->diag_data.dist5D.n_r * sim->diag_data.dist5D.n_phi * sim->diag_data.dist5D.n_z * sim->diag_data.dist5D.n_ppara * sim->diag_data.dist5D.n_pperp * sim->diag_data.dist5D.n_time * sim->diag_data.dist5D.n_q], \
		      sim->diag_data.dist6D.histogram[0:sim->diag_data.dist6D.n_r * sim->diag_data.dist6D.n_phi * sim->diag_data.dist6D.n_z * sim->diag_data.dist6D.n_pr * sim->diag_data.dist6D.n_pphi * sim->diag_data.dist6D.n_pz * sim->diag_data.dist6D.n_time * sim->diag_data.dist6D.n_q], \
		      sim->diag_data.distrho5D.histogram[0:sim->diag_data.distrho5D.n_rho * sim->diag_data.distrho5D.n_theta * sim->diag_data.distrho5D.n_phi * sim->diag_data.distrho5D.n_ppara * sim->diag_data.distrho5D.n_pperp * sim->diag_data.distrho5D.n_time * sim->diag_data.distrho5D.n_q], \
		      sim->diag_data.distrho6D.histogram[0:sim->diag_data.distrho6D.n_rho*sim->diag_data.distrho6D.n_theta*sim->diag_data.distrho6D.n_phi*sim->diag_data.distrho6D.n_pr*sim->diag_data.distrho6D.n_pphi*sim->diag_data.distrho6D.n_pz*sim->diag_data.distrho6D.n_time*sim->diag_data.distrho6D.n_q], \
		      sim->diag_data.distCOM.histogram[0:sim->diag_data.distCOM.n_mu*sim->diag_data.distCOM.n_Ekin*sim->diag_data.distCOM.n_Ptor], \
		      sim->wall_data.w2d.wall_r[0:sim->wall_data.w2d.n],sim->wall_data.w2d.wall_z[0:sim->wall_data.w2d.n],sim->wall_data.w3d.wall_tris[0:sim->wall_data.w3d.n*9+9],sim->wall_data.w3d.tree_array[0:sim->wall_data.w3d.tree_array_size], \
		      sim->plasma_data.plasma_1D.mass       [0:MAX_SPECIES],\
		      sim->plasma_data.plasma_1D.charge     [0:MAX_SPECIES],\
		      sim->plasma_data.plasma_1D.anum       [0:MAX_SPECIES],\
		      sim->plasma_data.plasma_1D.znum       [0:MAX_SPECIES],\
		      sim->plasma_data.plasma_1D.rho        [0:sim->plasma_data.plasma_1D.n_rho],\
		      sim->plasma_data.plasma_1D.temp       [0:sim->plasma_data.plasma_1D.n_rho*sim->plasma_data.plasma_1D.n_species], \
  		      sim->plasma_data.plasma_1D.dens       [0:sim->plasma_data.plasma_1D.n_rho*sim->plasma_data.plasma_1D.n_species], \
		      sim->plasma_data.plasma_1Dt.mass      [0:MAX_SPECIES],\
		      sim->plasma_data.plasma_1Dt.charge    [0:MAX_SPECIES],\
		      sim->plasma_data.plasma_1Dt.anum      [0:MAX_SPECIES],\
		      sim->plasma_data.plasma_1Dt.znum      [0:MAX_SPECIES],\
		      sim->plasma_data.plasma_1Dt.rho       [0:sim->plasma_data.plasma_1Dt.n_rho],\
		      sim->plasma_data.plasma_1Dt.temp      [0:sim->plasma_data.plasma_1Dt.n_time*sim->plasma_data.plasma_1Dt.n_rho*sim->plasma_data.plasma_1Dt.n_species],\
		      sim->plasma_data.plasma_1Dt.dens      [0:sim->plasma_data.plasma_1Dt.n_rho*sim->plasma_data.plasma_1Dt.n_species*sim->plasma_data.plasma_1Dt.n_time],\
		      sim->plasma_data.plasma_1Dt.time      [0:sim->plasma_data.plasma_1Dt.n_time],\
		      sim->plasma_data.plasma_1DS.mass      [0:MAX_SPECIES],\
		      sim->plasma_data.plasma_1DS.charge    [0:MAX_SPECIES],\
		      sim->plasma_data.plasma_1DS.anum      [0:MAX_SPECIES],\
		      sim->plasma_data.plasma_1DS.znum      [0:MAX_SPECIES],\
		      sim->plasma_data.plasma_1DS.temp      [0:2],\
		      sim->plasma_data.plasma_1DS.dens      [0:MAX_SPECIES],\
		      sim->plasma_data.plasma_1DS.temp[0].c [0:sim->plasma_data.plasma_1DS.temp[0].n_x*NSIZE_COMP1D],\
		      sim->plasma_data.plasma_1DS.temp[1].c [0:sim->plasma_data.plasma_1DS.temp[1].n_x*NSIZE_COMP1D],\
		      Bdata[0:1],Bdata->BTC.B[0:3],Bdata->BTC.dB[0:9],\
		      Bdata->BSTS.axis_r, Bdata->BSTS.axis_r.c [0:Bdata->BSTS.axis_r.n_x                                                           ],\
		      Bdata->BSTS.axis_z, Bdata->BSTS.axis_z.c [0:Bdata->BSTS.axis_z.n_x                                                           ],\
		      Bdata->BSTS.psi,    Bdata->BSTS.psi.c    [0:Bdata->BSTS.psi.n_x   *Bdata->BSTS.psi.n_y   *Bdata->BSTS.psi.n_z   *NSIZE_COMP3D],\
		      Bdata->BSTS.B_r,    Bdata->BSTS.B_r.c    [0:Bdata->BSTS.B_r.n_x   *Bdata->BSTS.B_r.n_y   *Bdata->BSTS.B_r.n_z   *NSIZE_COMP3D],\
		      Bdata->BSTS.B_z,    Bdata->BSTS.B_z.c    [0:Bdata->BSTS.B_z.n_x   *Bdata->BSTS.B_z.n_y   *Bdata->BSTS.B_z.n_z   *NSIZE_COMP3D],\
		      Bdata->BSTS.B_phi,  Bdata->BSTS.B_phi.c  [0:Bdata->BSTS.B_phi.n_x *Bdata->BSTS.B_phi.n_y *Bdata->BSTS.B_phi.n_z *NSIZE_COMP3D],\
		      Bdata->B3DS.psi,    Bdata->B3DS.psi.c    [0:Bdata->B3DS.psi.n_x   *Bdata->B3DS.psi.n_y                          *NSIZE_COMP2D],\
		      Bdata->B3DS.B_r,    Bdata->B3DS.B_r.c    [0:Bdata->B3DS.B_r.n_x   *Bdata->B3DS.B_r.n_y   *Bdata->B3DS.B_r.n_z   *NSIZE_COMP3D],\
		      Bdata->B3DS.B_phi,  Bdata->B3DS.B_phi.c  [0:Bdata->B3DS.B_phi.n_x *Bdata->B3DS.B_phi.n_y *Bdata->B3DS.B_phi.n_z *NSIZE_COMP3D],\
		      Bdata->B3DS.B_z,    Bdata->B3DS.B_z.c    [0:Bdata->B3DS.B_z.n_x   *Bdata->B3DS.B_z.n_y   *Bdata->B3DS.B_z.n_z   *NSIZE_COMP3D],\
		      Bdata->B2DS.psi,    Bdata->B2DS.psi.c    [0:Bdata->B2DS.psi.n_x   *Bdata->B2DS.psi.n_y                          *NSIZE_COMP2D],\
  		      Bdata->B2DS.B_r,    Bdata->B2DS.B_r.c    [0:Bdata->B2DS.B_r.n_x   *Bdata->B2DS.B_r.n_y                          *NSIZE_COMP2D],\
		      Bdata->B2DS.B_phi,  Bdata->B2DS.B_phi.c  [0:Bdata->B2DS.B_phi.n_x *Bdata->B2DS.B_phi.n_y                        *NSIZE_COMP2D],\
		      Bdata->B2DS.B_z,    Bdata->B2DS.B_z.c    [0:Bdata->B2DS.B_z.n_x   *Bdata->B2DS.B_z.n_y                          *NSIZE_COMP2D],\
		      Bdata->BGS.psi_coeff[0:13],				\
		      Edata[0:1],Edata->type,Edata->ETC,Edata->E1DS,Edata->ETC.Exyz[0:1],Edata->E1DS.dV,Edata->E1DS.dV.c[0:Edata->E1DS.dV.n_x*NSIZE_COMP1D], \
		      rnd[0:3*NSIMD] \
			)
    for (int i=0;i<MAX_SPECIES;i++) {
GPU_MAP_TO_DEVICE(
				  sim->plasma_data.plasma_1DS.dens[i].c[0:sim->plasma_data.plasma_1DS.dens[i].n_x*NSIZE_COMP1D] )
    }
}

void simulate_fo_fixed_copy_from_gpu(sim_data* sim, particle_simd_fo *p_ptr){

  GPU_UPDATE_FROM_DEVICE(
      p_ptr[0:1],p_ptr->running[0:NSIMD],p_ptr->r[0:NSIMD],p_ptr->phi[0:NSIMD],p_ptr->p_r[0:NSIMD],p_ptr->p_phi[0:NSIMD],p_ptr->p_z[0:NSIMD],p_ptr->mileage[0:NSIMD], \
  p_ptr->z[0:NSIMD],p_ptr->charge[0:NSIMD],p_ptr->mass[0:NSIMD],p_ptr->B_r[0:NSIMD],p_ptr->B_r_dr[0:NSIMD],p_ptr->B_r_dphi[0:NSIMD],p_ptr->B_r_dz[0:NSIMD], \
  p_ptr->B_phi[0:NSIMD],p_ptr->B_phi_dr[0:NSIMD],p_ptr->B_phi_dphi[0:NSIMD],p_ptr->B_phi_dz[0:NSIMD],p_ptr->B_z[0:NSIMD],p_ptr->B_z_dr[0:NSIMD],p_ptr->B_z_dphi[0:NSIMD], \
  p_ptr->B_z_dz[0:NSIMD],p_ptr->rho[0:NSIMD],p_ptr->theta[0:NSIMD],p_ptr->err[0:NSIMD],p_ptr->time[0:NSIMD],p_ptr->weight[0:NSIMD],p_ptr->cputime[0:NSIMD], \
      p_ptr->id[0:NSIMD],p_ptr->endcond[0:NSIMD],p_ptr->walltile[0:NSIMD],p_ptr->index[0:NSIMD],p_ptr->znum[0:NSIMD],p_ptr->anum[0:NSIMD],p_ptr->bounces[0:NSIMD] )

    GPU_MAP_FROM_DEVICE(
			      sim[0:1]  )
}
