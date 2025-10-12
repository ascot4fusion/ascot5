/**
 * @file ascot.h
 * Ascot functions and structures for external use.
 *
 * Ascot is a high-performance orbit-following code for fusion plasma physics
 * and engineering. In addition to marker tracing, it provides tools to
 * interpolate and process input data used in the simulation. The code can also
 * calculate fusion or NBI source, and magnetic field based on the coil
 * geometry.
 *
 * All performance sensitive components are implemented within. The code uses
 * OpenMP for parallelization on CPUs and OpenMP/OpenACC for GPU
 * parallelization. For CPUs, markers are parallelized both using threads and
 * vector operations.
 *
 * MPI level parallelization is left to the user of this library.
 */
#ifndef ASCOT_H
#define ASCOT_H

#include "datatypes.h"
#include "defines.h"
#include <stddef.h>
#include <stdint.h>

/**
 * Trace markers and solve the test particle Fokker-Planck equation.
 *
 * @param sim Simulation data with necessary inputs and diagnostics present.
 * @param nmrk Number of markers to be simulated.
 * @param mrk Markers to be simulated.
 */
void ascot_solve_distribution(Simulation *sim, size_t nmrk, State mrk[nmrk]);

/**
 * Calculate fusion source from two arbitrary ion distributions.
 *
 * Inputs and outputs are expected to have same (R, phi, z) or (rho, theta, phi)
 * dimensions.
 *
 * @param sim Simulation data with bfield and plasma inputs.
 * @param source Fusion source data that specifies the reactions and the
 *        reactants.
 * @param nsample Number of Monte Carlo samples to be used per one spatial cell.
 * @param product1 Output distribution for the first output species.
 * @param product2 Output distribution for the second output species.
 */
void ascot_solve_fusion(
    Simulation *sim, FusionSource *source, size_t nsample, histogram *product1,
    histogram *product2);

/**
 * Simulate NBI injection and calculate shinethrough and birth profile of NBI
 * ions.
 *
 * This function initializes neutrals at beamlet positions and traces them until
 * they have ionized or hit the wall.
 *
 * It is not necessary to have bfield to extend all the way to the beamlet as
 * the markers are assumed to remain neutrals on ballistic trajectories before
 * they enter the magnetic field data regime for the first time.
 *
 * Several injectors can be modelled simultaneously while keeping in mind that
 * the output does not contain any information on which injector a marker
 * originated.
 *
 * @param sim The simulation data with bfield, plasma, (wall), and any used
 *        diagnostics.
 * @param ninj Number of injectors.
 * @param injectors Neutral beam injector(s).
 * @param tlim Time interval over which the markers are injected.
 *        The marker initial time (at the point of injection) is sampled from
 *        this interval.
 * @param nmrk Number of markers to inject in total.
 * @param mrk Array where generated markers are stored.
 *        Contains all injected markers, not just those that end up ionized.
 */
void ascot_solve_nbi(
    Simulation *sim, size_t ninj, Nbi injectors[ninj], real tlim[2],
    size_t nmrk, State mrk[nmrk]);

/**
 * Evaluate magnetic field due to a coil at given points.
 *
 * The magnetic field is evaluated using Biot-Savart law. The coil geometry is
 * not assumed to be closed: the first and last points must coincide if the coil
 * is to be closed.
 *
 * @param npnt Number of query points.
 * @param ncoil Number of points in coil geometry.
 * @param coilxyz Coil geometry (x, y, z) coordinates [m].
 * @param x Query point x coordinate [m].
 * @param y Query point y coordinate [m].
 * @param z Query point z coordinate [m].
 * @param Bx Evaluated magnetic field x-component [T].
 * @param By Evaluated magnetic field y-component [T].
 * @param Bz Evaluated magnetic field z-component [T].
 */
void ascot_solve_field(
    size_t npnt, size_t ncoil, real coilxyz[3][ncoil], real xyz[3][npnt],
    real bxyz[3][npnt]);

/**
 * Interpolate input quantities at the given coordinates.
 *
 * @param bfield Magnetic field data.
 * @param efield Electric field data.
 * @param plasma Plasma data.
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
    Bfield *bfield, Efield *efield, Plasma *plasma, Neutral *neutral,
    Boozer *boozer, Mhd *mhd,
    // Atomic* atomic,
    size_t npnt, int modenumber, real r[npnt], real phi[npnt], real z[npnt],
    real t[npnt], real b[3][npnt], real bjac[9][npnt], real psi[4][npnt],
    real rho[2][npnt], real E[3][npnt], real n[][npnt], real T[2][npnt],
    real n0[][npnt], real T0[][npnt], real theta[4][npnt], real zeta[4][npnt],
    real alpha[5][npnt], real Phi[5][npnt], real mhd_b[3][npnt],
    real mhd_e[3][npnt], real mhd_phi[npnt]);

/**
 * Evaluate collision coefficients
 *
 * @param Simulation initialized simulation data struct
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
    Bfield *bfield, Plasma *plasma, size_t npnt, real ma, real qa, real R[npnt],
    real phi[npnt], real z[npnt], real t[npnt], real va[npnt],
    real coefficients[6][npnt], real clog[][npnt], real mu[3][npnt]);

/**
 * Return size of the ``real`` type.
 *
 * This can be used to determine if the code was compiled in single or double
 * precision.
 */
size_t ascot_sizeof_real(void);

/**
 * Map (rho, theta, phi) to (R, z) coordinates using Newton method.
 *
 * Due to padding, rho might not be exactly zero on the axis so we
 * return the axis position for small values of queried rho
 *
 * @param bfield Magnetic field data.
 * @param npnt Number of query points.
 * @param maxiter Maximum number of iterations in Newton algorithm before abort.
 * @param tol Maximum difference between two consecutive iterations before the
 *        value is accepted.
 * @param t Time instant when the inputs are interpolated [s].
 * @param rho Queried normalized poloidal flux [1].
 * @param theta Queried poloidal angle [rad].
 * @param phi Queried toroidal angle [rad].
 * @param rz Output array for (R, z) coordinates [m].
 */
void ascot_map_rhotheta_to_rz(
    Bfield *bfield, size_t npnt, size_t maxiter, real tol, real t,
    real rho[npnt], real theta[npnt], real phi[npnt], real rz[2][npnt]);

/**
 * Find psi on axis using the gradient descent method.
 *
 * Note that the psi value is not returned in case this algorithm fails.
 *
 * @param bfield magnetic field data
 * @param psi value of psi on axis if this function did not fail
 * @param rz initial (R,z) position where also the result is stored
 * @param step the step size
 * @param tol the current position is accepted if the distance (in meters)
 * between this and the previous point is below this value
 * @param maxiter maximum number of iterations before failure
 * @param ascent if true the algorithm instead ascends to find psi0 (> psi1)
 */
void ascot_find_psi_on_axis_2d(
    Bfield *bfield, size_t maxiter, size_t ascent, real step, real tol,
    real psi[1], real rz[2]);

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
    Bfield *bfield, size_t maxiter, int ascent, real phimin, real phimax,
    real step, real tol, real psi[1], real rzphi[3]);

#endif
