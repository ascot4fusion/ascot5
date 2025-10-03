/**
 * @file plasma.h
 * Plasma interface.
 *
 * This is an interface through which plasma data is initialized and accessed.
 * Reading e.g. from disk is done elsewhere.
 *
 * To add a new plasma instance, make sure these functions are implemented and
 * called from this interface, and that plasma.h contains enum type for the new
 * instance.
 *
 * The interface checks which instance given data corresponds to from the
 * "type"-field in plasma_offload_data or plasma_data that is given as an
 * argument, and calls the relevant function for that instance.
 *
 * Plasma species are referred by their index number: if a function returns e.g.
 * an array containing species densities these are ordered as electrons first
 * by different ion species in no particular order. However, the order is
 * invariant and consistent across all functions.
 */
#ifndef PLASMA_H
#define PLASMA_H

#include "ascot5.h"
#include "error.h"
#include "offload.h"

/**
 * @brief Plasma data types
 */
typedef enum plasma_type
{
    plasma_type_1D,  /**< Linear-interpolated 1D plasma data                */
    plasma_type_1Dt, /**< Linear-interpolated time-dependent 1D plasma data */
} plasma_type;

/**
 * 1D plasma parameters on the target.
 */
typedef struct
{
    int nrho;     /**< number of rho values in the data             */
    int nspecies; /**< number of plasma species including electrons */
    int *anum;    /**< ion species atomic number                    */
    int *znum;    /**< ion species charge number                    */
    real *mass;   /**< plasma species masses [kg]                   */
    real *charge; /**< plasma species charges [C]                   */
    real *rho;    /**< pointer to start of rho values               */
    real *temp;   /**< pointer to start of temperatures             */
    real *dens;   /**< pointer to start of densities                */
    real *vtor;   /**< pointer to start of toroidal rotation        */
} PlasmaLinear1D;

/**
 * @brief 1D plasma parameters on the target
 */
typedef struct
{
    int nrho;     /**< number of rho values in the data                */
    int ntime;    /**< number of time points                           */
    int nspecies; /**< number of plasma species including electrons    */
    int *anum;    /**< ion species atomic number                       */
    int *znum;    /**< ion species charge number                       */
    real *mass;   /**< plasma species masses [kg]                      */
    real *charge; /**< plasma species charges [C]                      */
    real *rho;    /**< pointer to start of rho values                  */
    real *time;   /**< pointer to start of time values                 */
    real *temp;   /**< pointer to start of temperatures                */
    real *dens;   /**< pointer to start of densities                   */
    real *vtor;   /**< pointer to start of toroidal rotation           */
} PlasmaDynamic1D;

/**
 * @brief Plasma simulation data
 *
 * The intended usage is that only single type of data is used at a time. This
 * is declared using the `type` field.
 */
typedef struct
{
    PlasmaLinear1D *linear1d;   /**< 1D data or NULL if not active           */
    PlasmaDynamic1D *dynamic1d; /**< 1D data or NULL if not active           */
    plasma_type type;           /**< Plasma data type wrapped by this struct */
} plasma_data;

/**
 * @brief Free allocated resources
 *
 * @param data pointer to the data struct
 */
void plasma_free(plasma_data *plasma);

/**
 * @brief Offload data to the accelerator.
 *
 * @param data pointer to the data struct
 */
void plasma_offload(plasma_data *plasma);

/**
 * @brief Evaluate plasma temperature
 *
 * This function evaluates the temperature of a given plasma species at the
 * given position.
 *
 * This is a SIMD function.
 *
 * @param temp array where evaluated temperature [J] is stored
 * @param rho normalized poloidal flux coordinate
 * @param r R-coordinate [m]
 * @param phi phi-coordinate [rad]
 * @param z z-coordinate [m]
 * @param t time coordinate [s]
 * @param species index of plasma species, 1 refers to electrons
 * @param plasma pointer to plasma data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(plasma)
a5err plasma_eval_temp(
    real *temp, real rho, real r, real phi, real z, real t, int species,
    plasma_data *plasma);
DECLARE_TARGET_END

/**
 * @brief Evaluate plasma density
 *
 * This function evaluates the density of a plasma species at the given
 * radial coordinate.
 *
 * This is a SIMD function.
 *
 * @param dens array where evaluated density will be stored
 * @param rho normalized poloidal flux coordinate
 * @param r R-coordinate [m]
 * @param phi phi-coordinate [rad]
 * @param z z-coordinate [m]
 * @param t time coordinate [s]
 * @param species index of plasma species, 1 refers to electrons
 * @param pls_data pointer to plasma data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(plasma)
a5err plasma_eval_dens(
    real *dens, real rho, real r, real phi, real z, real t, int species,
    plasma_data *plasma);
DECLARE_TARGET_END

/**
 * @brief Evaluate plasma density and temperature for all species
 *
 * This function evaluates the density and temperature of all plasma species at
 * the given radial coordinate.
 *
 * This is a SIMD function.
 *
 * @param dens pointer where density [m^-3] will be stored
 * @param temp pointer where temperature [J] will be stored
 * @param rho normalized poloidal flux coordinate
 * @param r R-coordinate [m]
 * @param phi phi-coordinate [rad]
 * @param z z-coordinate [m]
 * @param t time coordinate [s]
 * @param pls_data pointer to plasma data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(plasma)
a5err plasma_eval_densandtemp(
    real *dens, real *temp, real rho, real r, real phi, real z, real t,
    plasma_data *plasma);
DECLARE_TARGET_END

/**
 * @brief Evalate plasma flow along the field lines
 *
 * @param vflow pointer where the flow value is stored [m/s]
 * @param rho particle rho coordinate [1]
 * @param r particle R coordinate [m]
 * @param phi particle toroidal coordinate [rad]
 * @param z particle z coordinate [m]
 * @param t particle time coordinate [s]
 * @param pls_data pointer to plasma data
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(plasma)
a5err plasma_eval_flow(
    real *vflow, real rho, real r, real phi, real z, real t,
    plasma_data *plasma);
DECLARE_TARGET_END

/**
 * @brief Get the number of plasma species
 *
 * Retrieve the number of how many plasma species the data contains. The number
 * includes electrons.
 *
 * This is a SIMD function.
 *
 * @param pls_data pointer to plasma data struct
 *
 * @return The number of plasma species
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(plasma)
int plasma_get_n_species(plasma_data *plasma);
DECLARE_TARGET_END

/**
 * @brief Get mass of all plasma species
 *
 * Retrieve species' mass.
 *
 * This is a SIMD function.
 *
 * @param pls_data pointer to plasma data struct
 *
 * @return Pointer to array containing the requested mass values [kg]
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(plasma)
const real *plasma_get_species_mass(plasma_data *plasma);
DECLARE_TARGET_END

/**
 * @brief Get charge of all plasma species
 *
 * Retrieve species' charge.
 *
 * This is a SIMD function.
 *
 * @param pls_data pointer to plasma data struct
 *
 * @return Pointer to array containing the requested charge values [C]
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(plasma)
const real *plasma_get_species_charge(plasma_data *plasma);
DECLARE_TARGET_END

/**
 * @brief Get charge number of ion species
 *
 * This is a SIMD function.
 *
 * @param pls_data pointer to plasma data struct
 *
 * @return Pointer to array containing the charge numbers
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(plasma)
const int *plasma_get_species_znum(plasma_data *plasma);
DECLARE_TARGET_END

/**
 * @brief Get atomic mass number of ion species
 *
 * This is a SIMD function.
 *
 * @param pls_data pointer to plasma data struct
 *
 * @return Pointer to array containing the atomic mass numbers
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(plasma)
const int *plasma_get_species_anum(plasma_data *plasma);
DECLARE_TARGET_END

#endif
