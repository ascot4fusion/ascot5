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
#include "offload.h"
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
    sigmav_EII = 5, /* sigmav(E,T), electron-impact ionization                */
    sigmav_CX  = 6, /* sigmav(E,T), charge exchange                           */
    sigmav_BMS = 7, /* sigmav(E,Te,ne), beam-stopping coefficient             */
} asigma_reac_type;

/**
 * @brief Atomic reaction simulation data
 *
 * This struct holds data necessary for simulation. The struct is initialized
 * from input data in asigma_init().
 */
typedef struct {
    asigma_type type; /**< Atomic reaction data type wrapped by this struct */
    asigma_loc_data asigma_loc; /**< Local-files data or NULL if not active */
} asigma_data;

void asigma_free(asigma_data* data);
void asigma_offload(asigma_data* data);
void asigma_extrapolate(int extrapolate);

GPU_DECLARE_TARGET_SIMD_UNIFORM(asigmadata, z_1, a_1, mass, nspec, znum, anum)
a5err asigma_eval_cx(
    real* ratecoeff, int z_1, int a_1, real E, real mass, int nspec,
    const int* znum, const int* anum, real T_0, real* n_0,
    asigma_data* asigmadata);
GPU_DECLARE_TARGET_SIMD_UNIFORM(asigmadata, z_1, a_1, mass, nion, znum, anum)
a5err asigma_eval_bms(
    real* ratecoeff, int z_1, int a_1, real E, real mass, int nion,
    const int* znum, const int* anum, real T_e, real* n_i,
    asigma_data* asigmadata);
GPU_DECLARE_TARGET_SIMD_UNIFORM(asigmadata)
a5err asigma_eval_eii(
    real* ratecoeff, int z_1, int a_1, real E, real mass, real Te, real ne,
    asigma_data* asigmadata);

#endif
