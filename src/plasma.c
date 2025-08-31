/**
 * @file plasma.c
 * @brief Plasma interface
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
 *
 * Maximum number of plasma species is defined by MAX_SPECIES in ascot5.h.
 */
#include <stdio.h>
#include <stdlib.h>
#include "ascot5.h"
#include "error.h"
#include "print.h"
#include "plasma.h"
#include "plasma/plasma_1D.h"
#include "plasma/plasma_1DS.h"
#include "consts.h"

/**
 * @brief Free allocated resources
 *
 * @param data pointer to the data struct
 */
void plasma_free(plasma_data* data) {
    switch(data->type) {
        case plasma_type_1D:
            plasma_1D_free(data->plasma_1D);
            break;
        case plasma_type_1Dt:
            plasma_1Dt_free(data->plasma_1Dt);
            break;
        case plasma_type_1DS:
            plasma_1DS_free(data->plasma_1DS);
            break;
    }
}

/**
 * @brief Offload data to the accelerator.
 *
 * @param data pointer to the data struct
 */
void plasma_offload(plasma_data* data) {
    switch(data->type) {
        case plasma_type_1D:
            plasma_1D_offload(data->plasma_1D);
            break;
        case plasma_type_1Dt:
            plasma_1Dt_offload(data->plasma_1Dt);
            break;
        case plasma_type_1DS:
            plasma_1DS_offload(data->plasma_1DS);
            break;
    }
}

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
 * @param pls_data pointer to plasma data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err plasma_eval_temp(real* temp, real rho, real r, real phi, real z, real t,
                       int species, plasma_data* pls_data) {
    a5err err = 0;

    switch(pls_data->type) {
        case plasma_type_1D:
            err = plasma_1D_eval_temp(temp, rho, species, pls_data->plasma_1D);
            break;

        case plasma_type_1Dt:
            err = plasma_1Dt_eval_temp(temp, rho, t, species,
                                       pls_data->plasma_1Dt);
            break;

        case plasma_type_1DS:
            err = plasma_1DS_eval_temp(temp, rho, species,
                                       pls_data->plasma_1DS);
            break;
    }
    if(err) {
        /* In case of error, return some reasonable value to avoid further
           complications */
        temp[0] = 1e20;
    }

    return err;
}

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
a5err plasma_eval_dens(real* dens, real rho, real r, real phi, real z, real t,
                       int species, plasma_data* pls_data) {
    a5err err = 0;

    switch(pls_data->type) {
        case plasma_type_1D:
            err = plasma_1D_eval_dens(dens, rho, species,
                                      pls_data->plasma_1D);
            break;

        case plasma_type_1Dt:
            err = plasma_1Dt_eval_dens(dens, rho, t, species,
                                       pls_data->plasma_1Dt);
            break;

        case plasma_type_1DS:
            err = plasma_1DS_eval_dens(dens, rho, species,
                                       pls_data->plasma_1DS);
            break;
    }

    if(err) {
        /* In case of error, return some reasonable value to avoid further
           complications */
        dens[0] = 1e20;
    }
    return err;
}

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
a5err plasma_eval_densandtemp(real* dens, real* temp, real rho,
                              real r, real phi, real z, real t,
                              plasma_data* pls_data) {
    a5err err = 0;

    switch(pls_data->type) {
        case plasma_type_1D:
            err = plasma_1D_eval_densandtemp(dens, temp,
                                             rho, pls_data->plasma_1D);
            break;

        case plasma_type_1Dt:
            err = plasma_1Dt_eval_densandtemp(dens, temp,
                                              rho, t, pls_data->plasma_1Dt);
            break;

        case plasma_type_1DS:
            err = plasma_1DS_eval_densandtemp(dens, temp,
                                              rho, pls_data->plasma_1DS);
            break;

        default:
            /* Unregonized input. Produce error. */
            err = error_raise( ERR_UNKNOWN_INPUT, __LINE__, EF_PLASMA );
            break;
    }

    if(err) {
        /* In case of error, return some reasonable values to avoid further
           complications */
        for(int i=0; i<MAX_SPECIES; i++) {
            dens[i] = 1e20;
            temp[i] = 1e3;
        }
    }

    return err;
}

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
a5err plasma_eval_flow(
    real* vflow, real rho, real r, real phi, real z, real t,
    plasma_data* pls_data) {
    a5err err = 0;

    switch(pls_data->type) {
        case plasma_type_1D:
            err = plasma_1D_eval_flow(
                vflow, rho, r, pls_data->plasma_1D);
            break;

        case plasma_type_1DS:
            err = plasma_1DS_eval_flow(
                vflow, rho, r, pls_data->plasma_1DS);
            break;

        case plasma_type_1Dt:
            err = plasma_1Dt_eval_flow(
                vflow, rho, t, r, pls_data->plasma_1Dt);
            break;

        default:
            /* Unregonized input. Produce error. */
            err = error_raise( ERR_UNKNOWN_INPUT, __LINE__, EF_PLASMA );
            break;
    }

    if(err) {
        /* In case of error, return some reasonable values to avoid further
           complications */
        *vflow = 0;
    }

    return err;
}

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
int plasma_get_n_species(plasma_data* pls_data) {
    int n = 0;
    switch(pls_data->type) {
        case plasma_type_1D:
            n = pls_data->plasma_1D->n_species;
            break;

        case plasma_type_1Dt:
            n = pls_data->plasma_1Dt->n_species;
            break;

        case plasma_type_1DS:
            n = pls_data->plasma_1DS->n_species;
            break;
    }

    return n;
}

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
const real* plasma_get_species_mass(plasma_data* pls_data) {
    const real* mass = NULL;
    switch(pls_data->type) {
        case plasma_type_1D:
            mass = pls_data->plasma_1D->mass;
            break;

        case plasma_type_1Dt:
            mass = pls_data->plasma_1Dt->mass;
            break;

        case plasma_type_1DS:
            mass = pls_data->plasma_1DS->mass;
            break;
    }

    return mass;
}

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
const real* plasma_get_species_charge(plasma_data* pls_data) {
    const real* charge = NULL;
    switch(pls_data->type) {
        case plasma_type_1D:
            charge = pls_data->plasma_1D->charge;
            break;

        case plasma_type_1Dt:
            charge = pls_data->plasma_1Dt->charge;
            break;

        case plasma_type_1DS:
            charge = pls_data->plasma_1DS->charge;
            break;
    }

    return charge;
}

/**
 * @brief Get charge number of ion species
 *
 * This is a SIMD function.
 *
 * @param pls_data pointer to plasma data struct
 *
 * @return Pointer to array containing the charge numbers
 */
const int* plasma_get_species_znum(plasma_data* pls_data) {
    const int* znum = NULL;
    switch(pls_data->type) {
        case plasma_type_1D:
            znum = pls_data->plasma_1D->znum;
            break;

        case plasma_type_1Dt:
            znum = pls_data->plasma_1Dt->znum;
            break;

        case plasma_type_1DS:
            znum = pls_data->plasma_1DS->znum;
            break;
    }

    return znum;
}

/**
 * @brief Get atomic mass number of ion species
 *
 * This is a SIMD function.
 *
 * @param pls_data pointer to plasma data struct
 *
 * @return Pointer to array containing the atomic mass numbers
 */
const int* plasma_get_species_anum(plasma_data* pls_data) {
    const int* anum = NULL;
    switch(pls_data->type) {
        case plasma_type_1D:
            anum = pls_data->plasma_1D->anum;
            break;

        case plasma_type_1Dt:
            anum = pls_data->plasma_1Dt->anum;
            break;

        case plasma_type_1DS:
            anum = pls_data->plasma_1DS->anum;
            break;
    }

    return anum;
}
