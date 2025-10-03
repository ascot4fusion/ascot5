/**
 * @file rfof.h
 * Contains the function to be called during the simulation when using
 * ICRH.
 *
 * Requires librfof which contains the Fortran routines.
 */
#include "B_field.h"
#include "ascot5.h"
#include "particle.h"
#include <stdlib.h>

#ifndef RFOF_H
#define RFOF_H

/**
 * Reusable struct for storing marker specific data during the simulation loop.
 *
 * The data in this struct is altered during the simulation. Only pointers are
 * stored as the actual data is stored on the Fortran side.
 */
typedef struct rfof_marker
{
    void *p[NSIMD]; /**< The marker struct in a format required by librfof */
    void *history_array[NSIMD]; /**< Stores values of the resonance function for
                                     estimating the next time-step */
    void *diag_array[NSIMD]; /**< C equivalents of Fortran diagnostics pointers
                                  which are required but unused at the moment */
    int nrow[NSIMD]; /**< Number of rows in an resonance history matrix    */
    int ncol[NSIMD]; /**< Number of columns in an resonance history matrix */

} rfof_marker;

/** RFOF simulation input data.
 *
 * Immutable input data shared between all markers. The actual data is stored
 * in the Fortran side and this struct only stores the pointers.
 */
typedef struct
{
    void *rfof_input_params; /**< Pointer to rfof_input_param struct on
                                  the fortran side.                 */
    void *rfglobal;          /**< Wave field; same for all markers. */
} rfof_data;

/**
 * RFOF simulation output data.
 */
typedef struct rfof_output
{
    double dmu;          /**< Change in magnetic moment due to ICRH kick.     */
    double dvpar;        /**< Change in parallel velocity component due to ICRH
                              kick.                                           */
    double de;           /**< Change in energy due to a single ICRH kick [J]. */
    double deCumulative; /**< Change in energy due to possibly several ICRH
                              kicks during an orbit time step [J].            */
    double dpitch;       /**< Change in pitch due to ICRH kick.               */
    double maxAcc;       /**< Maximum acceleration allowed by RFOF.           */
    double RFdt;         /**< time step suggested by RFOF.                    */
} rfof_output;

/**
 * Initialise input data.
 *
 * Reads the ICRH (RFOF) inputs (xml, xsd, ASCII) and initialises the wave
 * field.
 *
 * @param rfof pointer to the RFOF data structure.
 */
void rfof_init(rfof_data *rfof);

/**
 * Deallocate the rfof_input_param struct on the fortran side.
 *
 * There exists only one copy of this struct and therefore it is to be
 * deallocated in the simulate.c after the loop is completed.
 */
void rfof_free(rfof_data *rfof);

/**
 * Initialises resonance history, diagnostics, and the marker struct.
 *
 * This function is to be called before the simulation loop.
 *
 * @param rfof_mrk Pointer to the local RFOF marker data.
 * @param rfof Pointer to the shared RFOF data.
 */
void rfof_set_up(rfof_marker *rfof_mrk, rfof_data *rfof_data);

/**
 * Deallocate the data structs used by the RFOF marker simulation data.
 *
 * This function is to be called after the simulation loop.
 *
 * @param rfof_mrk pointer to the RFOF marker simulation data.
 */
void rfof_tear_down(rfof_marker *rfof_mrk);

/**
 * Clear resonance history of an RFOF marker.
 *
 * History should be cleared whenever a marker finishes simulation. The new
 * marker cannot receive ICRH kicks during the first two time steps as its
 * resonance history must have at least two data points stored to estimate
 * the resonance location.
 *
 * @param rfof_mrk pointer to the RFOF marker simulation data.
 * @param imrk index in the NSIMD array for the marker whose history is cleared.
 */
void rfof_clear_history(rfof_marker *rfof_mrk, int imrk);

/**
 * Check if the marker is in resonance and apply kick.
 *
 * 1. Updates the fields of the rfof_marker based on the given input
 * ascot_marker.
 *
 * 2. Calls the "kick" function, which
 *      a) Checks resonance condition and
 *      b) if in resonance, kicks marker and updates velocity of the rfof marker
 *          and consequently also ascot marker (as only pointers are passed when
 *          creating the rfof marker).
 *
 * The time step can fail if the marker overshoots the resonance.
 *
 * @param p Pointer to marker simulation struct.
 * @param hin Current time step.
 * @param hout Suggestion for the next time step with negative sign
 * indicating a failed step.
 * @param rfof_mrk Pointer to the rfof marker simulation struct.
 * @param rfof_data Pointer to the shared rfof data.
 * @param bfield Pointer to the magnetic field data needed to evaluate psi.
 */
void rfof_resonance_check_and_kick_gc(
    particle_simd_gc *p, real *hin, real *hout_rfof, rfof_marker *rfof_mrk,
    rfof_data *rfof_data, B_field_data *bfield);

/**
 * @brief Explicitly set the coordinates in a marker struct.
 *
 * This routine is only used for testing and in libascot. Only the first marker
 * in the NSIMD array of rfof_marker is manipulated.
 *
 * @param rfof_mrk Pointer to the RFOF marker data.
 * @param id Marker identifier.
 * @param weight Marker weight [prt/s].
 * @param R Major radius [m].
 * @param phi Toroidal angle [rad].
 * @param z z coordinate [m].
 * @param psi Poloidal flux [Wb/rad].
 * @param charge Particle charge [C].
 * @param mass Particle mass [kg].
 * @param Ekin Kinetic energy [J].
 * @param vnorm Velocity [m/s].
 * @param mu Magnetic moment [J/T].
 * @param Pphi Canonical toroidal angular momentum [kg*m/s].
 * @param vpar Perpendicular velocity [m/s].
 * @param vperp Parallel velocity [m/s].
 * @param gyrof Gyrofrequency [rad/s].
 * @param vdriftRho Drift velocity in radial direction.
 * @param acc Accleration factor, should be 1.0 always.
 * @param is_accelerated Flag indicating whether acceleration is used, false
 * always.
 * @param is_already_allocated Flag whether to allocate a new marker struct,
 * should be false always.
 */
void rfof_set_marker_manually(
    rfof_marker *rfof_mrk, int *id, real *weight, real *R, real *phi, real *z,
    real *psi, real *charge, real *mass, real *Ekin, real *vnorm, real *mu,
    real *Pphi, real *vpar, real *vperp, real *gyrof, real *vdriftRho,
    real *acceleration, int *is_accelerated, int *is_already_allocated);

/**
 * Calculate the local E+ and E- values of the ICRH field.
 *
 * The definitions of E+ and E- sometimes differ. In this context,
 * E+ = E_LH * (cos(phi) + i sin(phi)),
 * where E_LH is the magnitude of the left-hand polarised (rotating) component
 * and phi is its phase. That is, E_LH and phi are real. Often, however, E_LH is
 * called E+ which creates confusion. Afterall, in this function, E+ is a
 * complex number.
 *
 * @param e_plus_real Re("E+"") component of the local wave field.
 * @param e_minus_real Re("E-"") component of the local wave field.
 * @param e_plus_imag Im("E+"") component of the local wave field.
 * @param e_minus_imag Im("E-"") component of the local wave field.
 * @param R Major radius coordinate [m].
 * @param z z-coordinate [m].
 * @param rfof Pointer to the RFOF data structure.
 */
void rfof_eval_rf_wave(
    real *e_plus_real, real *e_minus_real, real *e_plus_imag,
    real *e_minus_imag, real R, real z, rfof_data *rfof);

/**
 * Evaluate the value of resonance function (zero at the resonance).
 *
 * This function finds the closest resonance and in addition to evaluating the
 * resonance function it also returns the corresponding harmonic value.
 *
 * The resonance is evaluated for the first marker in the NSIMD array in
 * rfof_marker fields.
 *
 * @param cptr_marker Void pointer to the rfof_marker.
 * @param rfof_data Pointer to the RFOF data structure.
 * @param omega_res Evaluated value of the resonance function.
 * @param nharm The number of the closest harmonic found.
 */
void rfof_eval_resonance_function(
    real *omega_res, int *nharm, rfof_marker *rfof_mrk, rfof_data *rfof);

#endif
