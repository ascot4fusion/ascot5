/**
 * @file asigma.h
 * @brief Header file for asigma.c
 *
 * Contains a list declaring all atomic reaction data types, and declaration
 * of asigma_offload_data and asigma_data structs.
 */
#ifndef ASIGMA_H
#define ASIGMA_H

#include "ascot5.h"
#include "error.h"
#include "asigma/asigma_loc.h"

/**
 * @brief Atomic reaction data types
 *
 * Atomic reaction data types are used in the atomic reaction data interface
 * (asigma.c) to direct function calls to correct atomic reaction data
 * instances. Each atomic reaction data instance must have a
 * corresponding type.
 */
typedef enum asigma_type {
    asigma_type_loc, /**< Atomic reaction data from local files */
} asigma_type;

/**
 * @brief Reaction types for atomic reaction data
 *
 * The reaction type of atomic reactions is one of the reaction indentifier
 * parameters. It specifies the nature of the reaction and the form of the
 * reaction probability data.
 */
typedef enum asigma_reac_type {
    sigma_ioniz      = 1,  /* sigma(E), ionization (charge-increasing)    */
    sigma_recomb     = 2,  /* sigma(E), recombination (charge-decreasing) */
    sigma_CX         = 3,  /* sigma(E), charge exchange                   */
    sigmav_ioniz     = 4,  /* sigmav(E,T), ionization (charge-increasing) */
    sigmav_recomb    = 5,  /* sigmav(E,T), recombination (charge-decr.)   */
    sigmav_CX        = 6,  /* sigmav(E,T), charge exchange                */
    sigmav_BMS       = 7,  /* sigmav(E,Te,ne), beam-stopping coefficient  */
    sigmaveff_ioniz  = 8,  /* sigmav(n,T), eff. ioniz. (charge-incr.)     */
    sigmaveff_recomb = 9,  /* sigmav(n,T), eff. recomb. (charge-decr.)    */
    sigmaveff_CX     = 10  /* sigmav(n,T), effective charge exchange      */
} asigma_reac_type;

/**
 * @brief Atomic reaction simulation data
 *
 * This struct holds data necessary for simulation. The struct is initialized
 * from input data in asigma_init().
 */
typedef struct {
    asigma_loc_data* asigma_loc; /**< Local-files data or NULL if not active */
    asigma_type type;            /**< Atomic reaction data                   */
} asigma_data;

void asigma_free(asigma_data* data);
void asigma_offload(asigma_data* data);
void asigma_extrapolate(int extrapolate);
DECLARE_TARGET_SIMD_UNIFORM(asigmadata)
a5err asigma_eval_sigma(
    real* sigma, int z_1, int a_1, int z_2, int a_2, real E_coll_per_amu,
    asigma_reac_type reac_type, asigma_data* asigmadata);
DECLARE_TARGET_SIMD_UNIFORM(asigmadata)
a5err asigma_eval_sigmav(
    real* sigmav, int z_1, int a_1, real m_1, int z_2, int a_2,
    real E, real T_e, real T_0, real n_i, asigma_reac_type reac_type,
    asigma_data* asigmadata);
DECLARE_TARGET_SIMD_UNIFORM(asigmadata)
a5err asigma_eval_cx(
    real* ratecoeff, int z_1, int a_1, real E, real mass, int nspec,
    const int* znum, const int* anum, real T_0, real* n_0,
    asigma_data* asigmadata);
DECLARE_TARGET_SIMD_UNIFORM(asigmadata)
a5err asigma_eval_bms(
    real* ratecoeff, int z_1, int a_1, real E, real mass, int nion,
    const int* znum, const int* anum, real T_e, real* n_i,
    asigma_data* asigmadata);

#endif
