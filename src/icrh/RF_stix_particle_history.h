#ifndef RF_STIX_PARTICLE_HISTORY_H
#define RF_STIX_PARTICLE_HISTORY_H
#include "../offload.h"
#include "../ascot5.h"
#include "../error.h"
#include "../spline/interp.h" 
#include "../B_field.h"
#include "../boozer.h"
#include "../particle.h"


#define RF_N_HISTORY 3 // Number of history points for the RF particle history

typedef struct RF_particle_history{
    real dt[RF_N_HISTORY];    /** Particle history: time step*/
    real bnorm[RF_N_HISTORY]; /** Particle history: magnetic field strength */
    real rhopara[RF_N_HISTORY];  /** Particle history: parallel momentum */
    real R[RF_N_HISTORY];  /** Particle history: radial coordinate */
    real *resn;  /** Resonance locations for the Stix diffusion operator ( - kpara * vpara + omega / ntor)*/
    real *resp;  /** Resonance locations for the Stix diffusion operator ( - kpara * vpara + omega / ntor) */
    int nwaves; /** Number of waves in the Stix diffusion operator */
    int lhigh;  /** Maximum number of resonances to store. */
    real* omega; /** Frequencies of the waves */
    int* ntor;  /** Toroidal mode numbers of the waves */
    real qm;  /** Particle history: magnetic moment */
} RF_particle_history;

GPU_DECLARE_TARGET_SIMD_UNIFORM(hist, nwaves, omega, ntor, lhigh)
void RF_particle_history_init(RF_particle_history* hist, particle_simd_gc* p, int nwaves, real h,
                              int imrk, real* omega, int* ntor, int lhigh);
GPU_DECLARE_TARGET_SIMD_UNIFORM_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(hist)
void RF_particle_eval_nkicks(RF_particle_history* hist, particle_simd_gc* p, 
                             int imrk, int iwave, int *nkicks,
                             int *lres);
GPU_DECLARE_TARGET_SIMD_UNIFORM_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(hist)
void RF_particle_history_free(RF_particle_history* hist);
GPU_DECLARE_TARGET_SIMD_UNIFORM_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(hist)
void RF_particle_history_update(RF_particle_history* hist, particle_simd_gc* p, int imrk, real h);
GPU_DECLARE_TARGET_SIMD_UNIFORM_END

#endif // RF_STIX_PARTICLE_HISTORY_H