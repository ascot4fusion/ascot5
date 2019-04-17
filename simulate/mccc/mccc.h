/**
 * @file mccc.h
 * @brief Header file for mccc package
 */
#ifndef MCCC_H
#define MCCC_H

#include "../../ascot5.h"
#include "../../B_field.h"
#include "../../plasma.h"
#include "../../particle.h"
#include "../../random.h"
#include "mccc_wiener.h"

/**
 * @brief Defines minimum energy boundary condition
 *
 * This times local electron temperature is minimum energy boundary. If guiding
 * center energy goes below this, it is mirrored to prevent collision
 * coefficients from diverging.
 */
#define MCCC_CUTOFF 0.1

/**
 * @brief Parameters and data required to evaluate Coulomb collisions
 */
typedef struct {
    int usetabulated;   /**< Use tabulated values for special functions    */
    int include_energy; /**< Let collisions change energy                  */
    int include_pitch;  /**< Let collisions change pitch                   */
    int include_gcdiff; /**< Let collisions change guiding center position */
} mccc_data;

#pragma omp declare target

void mccc_init(mccc_data* mdata, int include_energy, int include_pitch,
               int include_gcdiff);
int mccc_eval_coefs(real ma, real qa, real r, real phi, real z, real t,
                    real* va, int nv, plasma_data* pdata, B_field_data* Bdata,
                    real* F, real* Dpara, real* Dperp, real* K, real* nu,
                    real* Q, real* dQ, real* dDpara, real* clog, real* mu0,
                    real* mu1, real* dmu0);
void mccc_fo_euler(particle_simd_fo* p, real* h, B_field_data* Bdata,
                   plasma_data* pdata, random_data* rdata, mccc_data* mdata);
void mccc_gc_euler(particle_simd_gc* p, real* h, B_field_data* Bdata,
                   plasma_data* pdata, random_data* rdata, mccc_data* mdata);
void mccc_gc_milstein(particle_simd_gc* p, real* hin, real* hout, real tol,
                      mccc_wienarr* wienarr, B_field_data* Bdata,
                      plasma_data* pdata, random_data* rdata, mccc_data* mdata);
#pragma omp end declare target

#endif
