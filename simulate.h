/**
 * @file simulate.h
 * @brief Simulation module interface
 *
 * Interface for calling different simulation modules, @see B_field.h
 */
#ifndef SIMULATE_H
#define SIMULATE_H

#include "ascot5.h"
#include "B_field.h"
#include "E_field.h"
#include "plasma_1d.h"
#include "distributions.h"
#include "wall.h"
#include "diag.h"

typedef struct {
    char hdf5fn[256];
    B_field_offload_data B_offload_data;
    E_field_offload_data E_offload_data;
    plasma_1d_offload_data plasma_offload_data;
    wall_offload_data wall_offload_data;
    diag_offload_data diag_offload_data;

    int sim_mode;
    int enable_ada;
    int record_GOasGC;

    int fix_usrdef_use;
    real fix_usrdef_val;
    int fix_stepsPerGO;

    real ada_tol_orbfol;
    real ada_tol_clmbcol;
    real ada_max_drho;
    real ada_max_dphi;
    real ada_max_acc;

    int enable_orbfol;
    int enable_clmbcol;
    
    int endcond_active;
    real endcond_maxSimTime;
    real endcond_maxCpuTime;
    real endcond_minRho;
    real endcond_maxRho;
    real endcond_minEkin;
    real endcond_minEkinPerTe;
    real endcond_maxTorOrb;
    real endcond_maxPolOrb;
    
} sim_offload_data;

typedef struct {
    
    B_field_data B_data;
    E_field_data E_data;
    plasma_1d_data plasma_data;
    wall_data wall_data;
    diag_data diag_data;

    int sim_mode;
    int enable_ada;
    int record_GOasGC;

    int fix_usrdef_use;
    real fix_usrdef_val;
    int fix_stepsPerGO;

    real ada_tol_orbfol;
    real ada_tol_clmbcol;
    real ada_max_drho;
    real ada_max_dphi;
    real ada_max_acc;

    int enable_orbfol;
    int enable_clmbcol;
    
    int endcond_active;
    real endcond_maxSimTime;
    real endcond_maxCpuTime;
    real endcond_minRho;
    real endcond_maxRho;
    real endcond_minEkin;
    real endcond_minEkinPerTe;
    real endcond_maxTorOrb;
    real endcond_maxPolOrb;
    
} sim_data;


#pragma omp declare target
void sim_init(sim_data* sim, sim_offload_data* offload_data);
void simulate(int id, int n_particles, input_particle* p,
              sim_offload_data* offload_data,
              real* B_offload_array,
              real* E_offload_array,
              real* plasma_offload_array,
              real* wall_offload_array,
              real* diag_offload_array);

void simulate_begin(int id, int n_particles, input_particle* p,
		    sim_offload_data* offload_data, sim_data* sim,
		    real* B_offload_array,
		    real* E_offload_array,
		    real* plasma_offload_array,
		    real* wall_offload_array,
		    real* diag_offload_array);

void simulate_continue(int id, int n_particles, input_particle* p,
		       sim_data* sim);

void simulate_hybrid(int id, int n_particles, input_particle* p,
		     sim_data* sim);

void simulate_end(sim_data* sim);
#pragma omp end declare target

#endif
