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
#include "plasma.h"
#include "neutral.h"
#include "wall.h"
#include "diag.h"
#include "offload.h"

enum {
    simulate_mode_fo = 1,
    simulate_mode_gc = 2,
    simulate_mode_hybrid = 3,
    simulate_mode_ml = 4
};

typedef struct {
    char hdf5_in[256];
    char hdf5_out[256];
    char outfn[256];
    char qid[256];

    int mpi_rank;
    int mpi_size;

    B_field_offload_data B_offload_data;
    E_field_offload_data E_offload_data;
    plasma_offload_data plasma_offload_data;
    neutral_offload_data neutral_offload_data;
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
    plasma_data plasma_data;
    neutral_data neutral_data;
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
void simulate(int id, int n_particles, particle_state* p,
              sim_offload_data* sim_offload,
              offload_package* offload_data,
              real* offload_array,
              real* diag_offload_array);
#pragma omp end declare target

#endif
