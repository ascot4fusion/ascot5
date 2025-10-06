/**
 * @file libascot.c
 * @brief Library of Ascot5 functions for external use.
 *
 * Functions in this file allows to evaluate input data and quantities using
 * the same methods as is used in actual simulation.
 */
#ifndef ASCOT_H
#define ASCOT_H

#include "defines.h"
#include "atomic.h"
#include "bfield.h"
#include "boozer.h"
#include "efield.h"
#include "mhd.h"
#include "neutral.h"
#include "plasma.h"
#include "particle.h"
#include "simulate.h"
#include "wall.h"
#include "fusion_source.h"

#define STORE(index, val, ptr)                                                 \
    if ((ptr))                                                                 \
    (ptr)[(index)] = (val)
/** Store value in output array if the output array is allocated. */


void ascot_solve_distribution(int nmarkers, particle_state* p, sim_data* sim);

/**
 * Calculate fusion source from two arbitrary ion distributions.
 *
 * Inputs and outputs are expected to have same physical (R, phi, z)
 * dimensions.
 *
 * @param sim Pointer to simulation data.
 * @param reaction Fusion reaction type, see the description.
 * @param n Number of Monte Carlo samples to be used.
 * @param react1 Reactant 1 distribution data.
 * @param react2 Reactant 2 distribution data.
 * @param mult Factor by which the output is scaled.
 * @param prod1_data Distribution data for product 1 output.
 * @param prod2_data Distribution data for product 2 output.
 */
void ascot_solve_fusion(
    sim_data *sim, afsi_data *data, size_t n, histogram *prod1, histogram *prod2);

/**
 * Simulate NBI injection.
 *
 * This function initializes neutrals and traces them until they have ionized or
 * hit the wall.
 *
 * @param sim Pointer to the simulation data structure.
 * @param nprt Number of markers to be injected.
 * @param t1 Time instant when the injector is turned on.
 * @param t2 Time instant when the injector is turned off.
 * @param p Pointer to the marker array which is allocated here.
 */
void ascot_solve_nbi(
    sim_data *sim, int nprt, real t1, real t2, particle_state **p);

/**
 * @brief Interpolate input quantities at the given coordinates.
 *
 * @param bfield Magnetic field data.
 * @param efield Electric field data.
 * @param plasma plasma data.
 * @param neutral Neutral data.
 * @param boozer Boozer data.
 * @param mhd MHD data.
 * @param npnt Number of evaluation points.
 * @param modenumber Evaluate mhd perturbation only for this mode.
 *        If modenumber < 0, evaluate all modes.
 * @param R R coordinates where the inputs are interpolated [m].
 * @param phi phi coordinates where the inputs are interpolated [rad].
 * @param z z coordinates where the inputs are interpolated [m].
 * @param t time coordinates where the inputs are interpolated [s].
 * @param B Magnetic field vector (Br, Bphi, Bz) [T].
 * @param Bjac Magnetic field Jacobian (dBr/dr, dBr/dphi, dBr/dz, dBphi/dr,
 *        dBphi/dphi, dBphi/dz, dBz/dr, dBz/dphi, dBz/dz) [T/m].
 * @param psi Poloidal flux and its derivatives (psi, dpsi/dr, dpsi/dphi,
 *        dpsi/dz) [Wb/rad, Wb/rad m].
 * @param rho Square root of the normalized poloidal flux and its
 *        derivative (rho, drho/dpsi) [1, rad/Wb].
 * @param E Electric field vector (Er, Ephi, Ez) [V/m].
 * @param n Density of plasma species (ne, ni1, ni2, ...) [m^-3].
 * @param T Temperature of plasma species (Te, Ti) [eV].
 * @param n0 Density of neutral species (n1, n2, ...) [m^-3].
 * @param T0 Temperature of neutral species (T1, T2, ...) [eV].
 * @param theta Poloidal Boozer coordinate and its derivatives (theta,
 *        dtheta/dr, dtheta/dphi, dtheta/dz) [rad, rad/m].
 * @param zeta output array [rad].
 * @param alpha Magnetic perturbation eigenfunction and its derivatives
 *        (Phi, dPhi/dr, dPhi/dphi, dPhi/dz) [1]
 * @param Phi Electric perturbation eigenfunction and its derivatives
 *        (Phi, dPhi/dr, dPhi/dphi, dPhi/dz) [1].
 * @param mhd_br Magnetic field perturbation components due to MHD
 *        (Br, Bphi, Bz) [T].
 * @param mhd_er Electric field perturbation components due to MHD
 *        (Er, Ephi, Ez) [V/m].
 * @param mhd_phi Electric field perturbation potential [V/m].
 */
void ascot_interpolate(
    Bfield *bfield, Efield *efield, Plasma *plasma,
    Neutral *neutral, Boozer *boozer, Mhd *mhd,
    // Atomic* atomic,
    int npnt, int modenumber, real R[npnt], real phi[npnt], real z[npnt],
    real t[npnt], real B[3][npnt], real Bjac[9][npnt], real psi[4][npnt],
    real rho[2][npnt], real E[3][npnt], real n[][npnt], real T[2][npnt],
    real n0[][npnt], real T0[][npnt], real theta[4][npnt], real zeta[4][npnt],
    real alpha[5][npnt], real Phi[5][npnt], real mhd_b[3][npnt],
    real mhd_e[3][npnt], real mhd_phi[npnt]);

/**
 * @brief Get magnetic axis at given coordinates.
 *
 * @param bfield Magnetic field data
 * @param nphi Number of evaluation points.
 * @param phi phi coordinates of the evaluation points [rad].
 * @param Raxis output array for axis R coordinates.
 */
void ascot_eval_axis(
    Bfield *bfield, int nphi, real phi[nphi], real axisRz[2][nphi]);

/**
 * Evaluate collision coefficients
 *
 * @param sim_data initialized simulation data struct
 * @param Neval number of evaluation points
 * @param R R coordinates of the evaluation points [m]
 * @param phi phi coordinates of the evaluation points [rad]
 * @param z z coordinates of the evaluation points [m]
 * @param t time coordinates of the evaluation points [s]
 * @param Nv number of velocity grid points
 * @param va test particle velocities at the evaluation point [m/s]
 * @param ma test particle mass [kg]
 * @param qa test particle charge [C]
 * @param coefficients The collision coefficients (F, Dpara, Dperp, K, nu, Q,
 *        dQ, dDpara, clog) [].
 * @param clog The Coulomb logarithm.
 * @param mu The special functions (mu0, mu1, dmu0)
 */
void ascot_eval_collcoefs(
    Bfield *bfield, Plasma *plasma, int npnt, real ma, real qa,
    real R[npnt], real phi[npnt], real z[npnt], real t[npnt], real va[npnt],
    real coefficients[6][npnt], real clog[][npnt], real mu[3][npnt]);


/**
 * Find one psi minimum using the gradient descent method for 3D field
 * inside a sector (phimin < phi < phimax). Not guaranteed to find the global
 * minimum inside the given sector.
 *
 * @param bfield magnetic field data
 * @param psi value of psi on axis if this function did not fail
 * @param rzphi initial (R,z,phi) position where also the result is stored
 * @param step the step size
 * @param tol the current position is accepted if the distance (in meters)
 * between this and the previous point is below this value
 * @param maxiter maximum number of iterations before failure
 * @param ascent if true the algorithm instead ascends to find psi0 (> psi1)
 */
void ascot_find_psi_on_axis_3d(
    Bfield *bfield, int maxiter, int ascent, real phimin, real phimax,
    real step, real tol, real psi[1], real rzphi[3]);



#endif
