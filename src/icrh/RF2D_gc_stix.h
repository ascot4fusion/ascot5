/**
 * @author Pablo Oyola oyola@pppl.gov
 * @file RF2D_gc_stix.h
 * @brief Header file for the kick-operator for the Stix diffusion operator
 *
 * Contains the declaration of the input structure for the evaluation
 * of the RF fields for the Stix diffusion operator.
 */

#ifndef RF2D_GC_STIX_H
#define RF2D_GC_STIX_H

#include "../offload.h"
#include "../ascot5.h"
#include "../error.h"
#include "../spline/interp.h" 
#include "../B_field.h"
#include "../boozer.h"
#include "../particle.h"
#include "RF_stix_particle_history.h"
#include <hdf5.h>

typedef unsigned char uint8;

// Defining the object structure of the RF fields.
typedef struct RF2D_gc_stix{
    // Physical parameters.
    real rmin;  /** Minimum radial coordinate */
    real rmax;  /** Maximum radial coordinate */
    real zmin;  /** Minimum vertical coordinate */
    real zmax;  /** Maximum vertical coordinate */
    int nr;     /** Number of radial grid points */
    int nz;     /** Number of vertical grid points */

    // 2D wave field.
    interp2D_data* Eplus_2;  /** Interpolation object for |E_+^2|*/
    interp2D_data* Eminus_2; /** Interpolation object for |E_-^2|*/
    interp2D_data* E2cross;  /** Interpolation object for Re(E_+^2 cnjt(E_-)) */
    interp2D_data* kperp;    /** Local value of the perpendicular vector */

    // Outer resources.
    B_field_data* bdata;       /** Pointer to the B_data structure containing the magnetic field data */
    
    // Wave properties.
    real* omega; /** Frequency/ies of the wave(s) */
    int* ntor;  /** Local value of the parallel wavenumber(s) */
    int nwaves;  /** Number of waves in the Stix diffusion operator */

    // Physical flags.
    int include_Eminus;       /** Flag to include E_- in the Stix diffusion operator */
    int include_stochastic;   /** Flag to include stochastic diffusion in the Stix operator */
    int include_vpara_kick;   /** Flag to include vpara kick in the Stix operator */
    int include_phase_factor; /** Flag to include phase factor in the Stix operator */

    // Internal flags.
    int enabled; /** Flag to check if the object has been initialized */
    
    // Resonance locations.
    int n_max_res;      /** Maximum number of resonance locations */
    real** R_resonances; /** Radial location of the cold resonance condition */
    int* nres;          /** Number of resonance locations */
    int** res_nums;      /** Resonance levels. */
} RF2D_gc_stix;


a5err RF2D_gc_stix_init_from_file(RF2D_gc_stix* stix_data, hid_t f, char* qid,
                        int lhigh, B_field_data* bdata); // Done
a5err RF2D_gc_stix_init(RF2D_gc_stix* stix_data,
                        real* Eplus_re, real* Eplus_im,
                        real* Eminus_re, real* Eminus_im,
                        real* kperp, real* costerm, real* sinterm,
                        real rmin, real rmax, real zmin, real zmax,
                        int nr, int nz, int nwaves,
                        int include_Eminus, int include_stochastic,
                        int include_vpara_kick, int include_phase_factor,
                        B_field_data* bdata);

void RF2D_gc_stix_free(RF2D_gc_stix* stix_data);
void RF2D_gc_stix_offload(RF2D_gc_stix* stix_data);
real guess_qm(particle_queue* pq);

// Internal useful functions.
a5err RF2D_gc_stix_compute_cold_resonances(RF2D_gc_stix* stix_data,
                                           B_field_data* bdata, int n_max_res,
                                           real qm); 

void RF2D_gc_stix_scatter(RF2D_gc_stix* stix, RF_particle_history* hist, 
                          particle_simd_gc* p, real* h, real* rnd, uint8* used); 


GPU_DECLARE_TARGET_SIMD_UNIFORM(stix_data)
real RF2D_gc_stix_get_interaction_time(RF2D_gc_stix* stix_data, 
                                       RF_particle_history* hist, 
                                       int iwave, int l);
GPU_DECLARE_TARGET_SIMD_UNIFORM_END


GPU_DECLARE_TARGET_SIMD_UNIFORM(stix_data)
a5err RF2D_gc_stix_eval_fields(real r, real phi, real z, real t,
                               RF2D_gc_stix* stix_data,
                               real* Eplus_2, real* Eminus_2,
                               real* E2cross, real* kperp);
GPU_DECLARE_TARGET_SIMD_UNIFORM_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(stix_data, used, value)
void set_rndusage(RF2D_gc_stix* stix_data, uint8* used, int imrk, uint8 value);
GPU_DECLARE_TARGET_SIMD_UNIFORM_END

#endif