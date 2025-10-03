/**
 * @file B_STS.h
 * Spline-interpolated stellarator magnetic field implementation.
 *
 * Both poloidal flux and magnetic field is evaluated by interpolating
 * uniform 3D tables in which the data is provided. In contrast to tokamak
 * fields, the poloidal component of the evaluated magnetic field comes
 * exclusively from the tabulated data, with no contribution from the poloidal
 * flux.
 *
 * The magnetic axis location may change along the torus, and for this the
 * location of the axis is tabulated and linearly interpolated as a function of
 * the toroidal angle.
 */
#ifndef B_STS_H
#define B_STS_H
#include "B_field.h"
#include "ascot5.h"
#include "error.h"
#include "offload.h"

/**
 * Initialize the 3D spline-interpolated tokamak magnetic field data.
 *
 * Assigns the fields in the struct with the provided values, and initializes
 * and allocates the spline-interpolants.
 *
 * @param bfield The struct to initialize.
 * @param pnr Number of R grid points in ``psi``.
 * @param pnz Number of z grid points in ``psi``.
 * @param bnr Number of R grid points in ``br``, ``bz``, and ``bphi``.
 * @param bnz Number of z grid points in ``br``, ``bz``, and ``bphi``.
 * @param bnphi Number of phi grid points in ``br``, ``bz``, and ``bphi``.´
 * @param naxis Number of grid points in ``axisr`` and ``axisz``.
 * @param prlim Limits of the uniform R abscissa in ``psi`` [m].
 * @param pzlim Limits of the uniform z abscissa in ``psi`` [m].
 * @param brlim Limits of the uniform R abscissa in ``br``, ``bz``,
 *        and ``bphi`` [m].
 * @param bzlim Limits of the uniform z abscissa in ``br``, ``bz``,
 *        and ``bphi`` [m].
 * @param bphilim Limits of the uniform phi abscissa in ``br``, ``bz``,
 *        and ``bphi`` [rad].
 * @param axislim Limits of the uniform phi abscissa in ``axisr`` and
 *        ``axisz`` [m].
 * @param axisr Tabulated magnetic axis R coordinates [m].
 * @param axisz Tabulated magnetic axis z coordinates [m].
 * @param psilimits Poloidal flux at axis and separatrix [Wb/rad].
 * @param psi Tabulated values of poloidal flux [Wb/rad].
 *        Layout: (Ri, phij, zk) = [k*nr*nphi + j*nr + i] (C order).
 * @param br Tabulated values of R component of B [T].
 *        Layout: (Ri, phij, zk) = [k*nr*nphi + j*nr + i] (C order).
 * @param bz Tabulated values of z component of B [T].
 *        Layout: (Ri, phij, zk) = [k*nr*nphi + j*nr + i] (C order).
 * @param bphi Tabulated values of phi component of B [T].
 *        Layout: (Ri, phij, zk) = [k*nr*nphi + j*nr + i] (C order).
 *
 * @return  Zero if the initialization succeeded.
 */
int BfieldStellarator_init(
    BfieldStellarator *bfield, int pnr, int pnz, int pnphi, int bnr, int bnz,
    int bnphi, int naxis, real prlim[2], real pzlim[2], real pphilim[2],
    real brlim[2], real bzlim[2], real bphilim[2], real axislim[2],
    real axisr[naxis], real axisz[naxis], real psilimits[2],
    real psi[pnr * pnz * pnphi], real br[bnr * bnz * bnphi],
    real bz[bnr * bnz * bnphi], real bphi[bnr * bnz * bnphi]);

/**
 * @brief Free allocated resources
 *
 * Spline and linear interpolants are freed.
 *
 * @param bfield The struct whose fields are deallocated.
 */
void BfieldStellarator_free(BfieldStellarator *bfield);

/**
 * @brief Offload data to the accelerator.
 *
 * @param bfield The struct to offload.
 */
void BfieldStellarator_offload(BfieldStellarator *bfield);

/**
 * Evaluate poloidal flux.
 *
 * @param psi Evaluated poloidal flux [Wb/rad].
 * @param r R coordinate of the query point [m].
 * @param phi phi coordinate of the query point [rad].
 * @param z z coordinate of the query point [m].
 * @param bfield The magnetic field data.
 *
 * @return Zero if the evaluation succeeded.
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(bfield)
a5err BfieldStellarator_eval_psi(
    real psi[1], real r, real phi, real z, BfieldStellarator *bfield);
DECLARE_TARGET_END

/**
 * Evaluate poloidal flux and its derivatives.
 *
 * @param psi_dpsi Evaluated poloidal flux and it's derivatives [Wb/rad].
 *        Layout: [psi, dpsi/dr, dpsi/dphi, dpsi/dz].
 * @param r R coordinate of the query point [m].
 * @param phi phi coordinate of the query point [rad].
 * @param z z coordinate of the query point [m].
 * @param bfield The magnetic field data.
 *
 * @return Zero if the evaluation succeeded.
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(bfield)
a5err BfieldStellarator_eval_psi_dpsi(
    real psi_dpsi[4], real r, real phi, real z, BfieldStellarator *bfield);
DECLARE_TARGET_END

/**
 * @brief Evaluate normalized poloidal flux and its derivatives.
 *
 * @param rho_drho Evaluated normalized poloidal flux and it's derivatives [1].
 *        Layout: [rho, drho/dr, drho/dphi, drho/dz].
 * @param r R coordinate of the query point [m].
 * @param phi phi coordinate of the query point [rad].
 * @param z z coordinate of the query point [m].
 * @param bfield The magnetic field data.
 *
 * @return Zero if the evaluation succeeded.
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(bfield)
a5err BfieldStellarator_eval_rho_drho(
    real rho_drho[4], real r, real phi, real z, BfieldStellarator *bfield);
DECLARE_TARGET_END

/**
 * Evaluate magnetic field vector.
 *
 * @param b Evaluated magnetic field vector [T].
 *        Layout: [br, bphi, bz].
 * @param r R coordinate of the query point [m].
 * @param phi phi coordinate of the query point [rad].
 * @param z z coordinate of the query point [m].
 * @param bfield The magnetic field data.
 *
 * @return Zero if the evaluation succeeded.
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(bfield)
a5err BfieldStellarator_eval_b(
    real b[3], real r, real phi, real z, BfieldStellarator *bfield);
DECLARE_TARGET_END

/**
 * Evaluate magnetic field vector and its derivatives.
 *
 * @param b_db Evaluated magnetic field vector and its derivatives [T].
 *        Layout: [br, dbr/dr, dbr/dphi, bz, dbz/dz, bphi, dbphi/dr, dbphi/dphi,
 *        dbphi/dz, bz, dbz/dr, dbz/dphi, dbz/dz].
 * @param r R coordinate of the query point [m].
 * @param phi phi coordinate of the query point [rad].
 * @param z z coordinate of the query point [m].
 * @param bfield The magnetic field data.
 *
 * @return Zero if the evaluation succeeded.
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(bfield)
a5err BfieldStellarator_eval_b_db(
    real b_db[12], real r, real phi, real z, BfieldStellarator *bfield);
DECLARE_TARGET_END

/**
 * Evaluate the magnetic axis (R, z) coordinates.
 *
 * Uses linear interpolation to find the axis location from tabulated values at
 * the given toroidal position.
 *
 * @param axisrz Evaluated axis coordinates [m].
 * @param phi phi coordinate of the query point [rad].
 * @param bfield The magnetic field data.
 *
 * @return Zero if the evaluation succeeded.
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(bfield)
a5err BfieldStellarator_eval_axisrz(
    real axisrz[2], real phi, BfieldStellarator *bfield);
DECLARE_TARGET_END
#endif
