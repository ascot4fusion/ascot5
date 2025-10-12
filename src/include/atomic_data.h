/**
 * @file atomic_data.h
 * Data structures for atomic data.
 */
#ifndef ATOMIC_DATA_H
#define ATOMIC_DATA_H

#include "interp_data.h"
#include <stddef.h>

/**
 * Reaction types for atomic reaction data.
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
 * Local-files atomic reaction simulation data.
 */
typedef struct {
    size_t N_reac;            /**< Number of reactions.                       */
    int* z_1;                 /**< Atomic number of test particle.            */
    int* a_1;                 /**< Mass number of test particle.              */
    int* z_2;                 /**< Atomic number of bulk particle.            */
    int* a_2;                 /**< Mass number of bulk particle.              */
    asigma_reac_type* reac_type;/**< Reaction type.                           */
    Spline1D* sigma;     /**< Spline of cross-sections.                  */
    Spline2D* sigmav;    /**< Spline of rate coefficients.               */
    Spline3D* BMSsigmav; /**< Spline of BMS rate coefficients.           */
} Atomic;

#endif
