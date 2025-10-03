/**
 * @file asigma.h
 * Atomic reaction data interface
 *
 * This is an interface through which atomic reaction data is initialized
 * and accessed. Reading, for example from local files, is done elsewhere.
 *
 * The name asigma is short for atomicsigma. The word sigma refers to
 * cross-section, a fundamental type of reaction probability data. Note
 * that the data is not necessarily in the form of pure cross-sections.
 * It might be in some derivative form, such as rate coefficients.
 *
 * To add a new atomic reaction data instance, make sure the functions
 * are implemented and called from this interface, and that asigma.h
 * contains enum type for the new instance.
 *
 * The interface checks which instance given data corresponds to from the
 * "type"-field in asigma_offload_data or asigma_data that is given
 * as an argument, and calls the relevant function for that instance.
 */
#ifndef ASIGMA_H
#define ASIGMA_H

#include "ascot5.h"
#include "error.h"
#include "interp.h"

/**
 * Atomic reaction data types.
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
 * Local-files atomic reaction simulation data.
 */
typedef struct {
    int N_reac;               /**< Number of reactions.                       */
    int* z_1;                 /**< Atomic number of test particle.            */
    int* a_1;                 /**< Mass number of test particle.              */
    int* z_2;                 /**< Atomic number of bulk particle.            */
    int* a_2;                 /**< Mass number of bulk particle.              */
    int* reac_type;           /**< Reaction type.                             */
    interp1D_data* sigma;     /**< Spline of cross-sections.                  */
    interp2D_data* sigmav;    /**< Spline of rate coefficients.               */
    interp3D_data* BMSsigmav; /**< Spline of BMS rate coefficients.           */
} asigma_loc_data;

/**
 * @brief Reaction types for atomic reaction data
 *
 * The reaction type of atomic reactions is one of the reaction indentifier
 * parameters. It specifies the nature of the reaction and the form of the
 * reaction probability data.
 */
typedef enum asigma_reac_type {
    sigma_ioniz      = 1,  /**< sigma(E), ionization (charge-increasing).     */
    sigma_recomb     = 2,  /**< sigma(E), recombination (charge-decreasing).  */
    sigma_CX         = 3,  /**< sigma(E), charge exchange.                    */
    sigmav_ioniz     = 4,  /**< sigmav(E,T), ionization (charge-increasing).  */
    sigmav_recomb    = 5,  /**< sigmav(E,T), recombination (charge-decr.).    */
    sigmav_CX        = 6,  /**< sigmav(E,T), charge exchange.                 */
    sigmav_BMS       = 7,  /**< sigmav(E,Te,ne), beam-stopping coefficient.   */
    sigmaveff_ioniz  = 8,  /**< sigmav(n,T), eff. ioniz. (charge-incr.).      */
    sigmaveff_recomb = 9,  /**< sigmav(n,T), eff. recomb. (charge-decr.).     */
    sigmaveff_CX     = 10  /**< sigmav(n,T), effective charge exchange.       */
} asigma_reac_type;

/**
 * @brief Atomic reaction simulation data
 *
 * This struct holds data necessary for simulation. The struct is initialized
 * from input data in asigma_init().
 */
typedef struct {
    asigma_loc_data* asigma_loc; /**< Local-files data or NULL if not active. */
    asigma_type type;            /**< Atomic reaction data.                   */
} asigma_data;

/**
 * Free allocated resources.
 *
 * @param data Pointer to the data struct.
 */
void asigma_free(asigma_data* atomic);

/**
 * Offload data to the accelerator.
 *
 * @param atomic pointer to the data struct.
 */
void asigma_offload(asigma_data* atomic);

/**
 * Toggle extrapolation when evaluating cross sections.
 *
 * In this context the extrapolation means values outside the abscissae
 * are set to zero instead of raising an error.
 *
 * @param extrapolate Flag whether to extrapolate.
 */
void asigma_extrapolate(int extrapolate);

/**
 * Evaluate atomic reaction cross-section.
 *
 * This function evaluates the cross-section (sigma) for the atomic reaction
 * corresponding to the reaction identifiers given as parameters at the
 * given mass-normalized collision energy.
 *
 * This is a SIMD function.
 *
 * @param sigma Pointer to evaluated cross-section.
 * @param z_1 Atomic number of fast particle.
 * @param a_1 Atomic mass number of fast particle.
 * @param z_2 Atomic number of bulk particle.
 * @param a_2 Atomic mass number of bulk particle.
 * @param E_coll_per_amu Energy per amu corresponding to collision speed.
 * @param reac_type Reaction type.
 * @param asigma_data Pointer to atomic data struct.
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise.
 */
DECLARE_TARGET_SIMD_UNIFORM(atomic)
a5err asigma_eval_sigma(
    real* sigma, int z_1, int a_1, int z_2, int a_2, real E_coll_per_amu,
    asigma_reac_type reac_type, asigma_data* atomic);

/**
 * Evaluate atomic reaction rate coefficient.
 *
 * This function evaluates the rate coefficient (<sigma*v>) for the atomic
 * reaction corresponding to the reaction identifiers given as parameters
 * at the given fast particle energy and bulk plasma conditions.
 *
 * This is a SIMD function.
 *
 * @param sigmav Pointer to evaluated rate coefficient.
 * @param z_1 Atomic number of fast particle.
 * @param a_1 Atomic mass number of fast particle.
 * @param m_1 Mass of fast particle.
 * @param z_2 Atomic number of bulk particle.
 * @param a_2 Atomic mass number of bulk particle.
 * @param E Energy of fast particle.
 * @param T_e Electron temperature of bulk plasma.
 * @param T_0 Temperature of bulk neutrals.
 * @param n_i Ion density of bulk plasma.
 * @param reac_type Reaction type.
 * @param asigma_data Pointer to atomic data struct.
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise.
 */
DECLARE_TARGET_SIMD_UNIFORM(atomic)
a5err asigma_eval_sigmav(
    real* sigmav, int z_1, int a_1, real m_1, int z_2, int a_2,
    real E, real T_e, real T_0, real n_i, asigma_reac_type reac_type,
    asigma_data* atomic);

/**
 * Evaluate charge exchange rate coefficient.
 *
 * This is a SIMD function.
 *
 * @param ratecoeff Pointer to evaluated rate coefficient.
 * @param z_1 Atomic number of fast particle.
 * @param a_1 Atomic mass number of fast particle.
 * @param E Energy of fast particle.
 * @param mass Mass of fast particle.
 * @param nspec Number of bulk neutral species.
 * @param znum Atomic numbers of bulk particles.
 * @param anum Atomic mass numbers of bulk particles.
 * @param T_0 Temperature of bulk neutrals.
 * @param n_0 Densities of bulk neutrals.
 * @param asigma_data Pointer to atomic data struct.
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise.
 */
DECLARE_TARGET_SIMD_UNIFORM(atomic)
a5err asigma_eval_cx(
    real* ratecoeff, int z_1, int a_1, real E, real mass, int nspec,
    const int* znum, const int* anum, real T_0, real* n_0,
    asigma_data* atomic);

/**
 * Evaluate beam stopping rate coefficient.
 *
 * This is a SIMD function.
 *
 * @param ratecoeff Pointer to evaluated rate coefficient.
 * @param z_1 Atomic number of fast particle.
 * @param a_1 Atomic mass number of fast particle.
 * @param E Energy of fast particle.
 * @param mass Mass of fast particle.
 * @param nion Number of bulk ion species.
 * @param znum Atomic numbers of bulk particles.
 * @param anum Atomic mass numbers of bulk particles.
 * @param T_e Electron temperature of bulk plasma.
 * @param n_i Densities of bulk ions.
 * @param asigma_data Pointer to atomic data struct.
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise.
 */
DECLARE_TARGET_SIMD_UNIFORM(atomic)
a5err asigma_eval_bms(
    real* ratecoeff, int z_1, int a_1, real E, real mass, int nion,
    const int* znum, const int* anum, real T_e, real* n_i,
    asigma_data* atomic);

#endif
