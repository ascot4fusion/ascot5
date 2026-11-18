/**
 * @file bfield_linear3d.h
 * Linear-interpolated magnetic field implementation.
 *
 * Poloidal flux is evaluated by interpolating an uniform 2D table while
 * magnetic field components are evaluated by interpolating uniform 3D
 * tables. The 2D and the 3D table have separate abscissae. The poloidal
 * component of the evaluated magnetic field includes the contribution from the
 * gradient of the poloidal flux.
 */
#ifndef BFIELD_LINEAR3D_H
#define BFIELD_LINEAR3D_H
#include "bfield.h"
#include "defines.h"
#include "parallel.h"

/**
 * Initialize the linearly interpolated 3D tokamak magnetic field data.
 *
 * Assigns the fields in the struct with the provided values. The tabulated
 * data (psi, br, bz, bphi) is not copied: the interpolants point directly to
 * the given arrays, which the caller must keep alive (and eventually free)
 * for as long as the struct is in use.
 *
 * @param bfield The struct to initialize.
 * @param pnr Number of R grid points in ``psi``.
 * @param pnz Number of z grid points in ``psi``.
 * @param bnr Number of R grid points in ``br``, ``bz``, and ``bphi``.
 * @param bnz Number of z grid points in ``br``, ``bz``, and ``bphi``.
 * @param bnphi Number of phi grid points in ``br``, ``bz``, and ``bphi``.
 * @param prlim Limits of the uniform R abscissa in ``psi`` [m].
 * @param pzlim Limits of the uniform z abscissa in ``psi`` [m].
 * @param brlim Limits of the uniform R abscissa in ``br``, ``bz``,
 *        and ``bphi`` [m].
 * @param bzlim Limits of the uniform z abscissa in ``br``, ``bz``,
 *        and ``bphi`` [m].
 * @param bphilim Limits of the uniform phi abscissa in ``br``, ``bz``,
 *        and ``bphi`` [rad].
 * @param axisrz Magnetic axis (R, z) coordinates [m].
 * @param psilimits Poloidal flux at axis and separatrix [Wb/rad].
 * @param psi Tabulated values of poloidal flux [Wb/rad].
 *        Layout: (Ri, zj) = [i*nz + j] (C order).
 * @param br Tabulated values of R component of B [T].
 *        Layout: (Ri, phij, zk) = [i*nz*nphi + j*nz + k] (C order).
 * @param bz Tabulated values of z component of B [T].
 *        Layout: (Ri, phij, zk) = [i*nz*nphi + j*nz + k] (C order).
 * @param bphi Tabulated values of phi component of B [T].
 *        Layout: (Ri, phij, zk) = [i*nz*nphi + j*nz + k] (C order).
 *
 * @return Zero on success.
 */
int BfieldLinear3D_init(
    BfieldLinear3D *bfield, size_t pnr, size_t pnz, size_t bnr, size_t bnz,
    size_t bnphi, real prlim[2], real pzlim[2], real brlim[2], real bzlim[2],
    real bphilim[2], real axisrz[2], real psilimits[2], real psi[pnr * pnz],
    real br[bnr * bnz * bnphi], real bz[bnr * bnz * bnphi],
    real bphi[bnr * bnz * bnphi]);

/**
 * Free allocated resources.
 *
 * This is a no-op since the tabulated data is owned by the caller of
 * BfieldLinear3D_init.
 *
 * @param bfield The struct whose fields are deallocated.
 */
void BfieldLinear3D_free(BfieldLinear3D *bfield);

/**
 * Offload data to the accelerator.
 *
 * @param bfield The struct to offload.
 */
void BfieldLinear3D_offload(BfieldLinear3D *bfield);

GPU_DECLARE_TARGET_SIMD_UNIFORM(bfield)
/**
 * Evaluate poloidal flux.
 *
 * @param psi Evaluated poloidal flux [Wb/rad].
 * @param r R coordinate of the query point [m].
 * @param z z coordinate of the query point [m].
 * @param bfield The magnetic field data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t BfieldLinear3D_eval_psi(
    real psi[1], real r, real z, BfieldLinear3D *bfield);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(bfield)
/**
 * Evaluate poloidal flux and its derivatives.
 *
 * @param psi_dpsi Evaluated poloidal flux and its derivatives [Wb/rad].
 *        Layout: [psi, dpsi/dr, dpsi/dphi, dpsi/dz].
 * @param r R coordinate of the query point [m].
 * @param z z coordinate of the query point [m].
 * @param bfield The magnetic field data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t BfieldLinear3D_eval_psi_dpsi(
    real psi_dpsi[4], real r, real z, BfieldLinear3D *bfield);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(bfield)
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
err_t BfieldLinear3D_eval_b(
    real b[3], real r, real phi, real z, BfieldLinear3D *bfield);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(bfield)
/**
 * Evaluate magnetic field vector and its derivatives.
 *
 * @param b_db Evaluated magnetic field vector and its derivatives [T].
 *        Layout: [br, bphi, bz, dbr/dr, dbr/dphi, dbrdz, dbphi/dr, dbphi/dphi,
 *        dbphi/dz, dbz/dr, dbz/dphi, dbz/dz].
 * @param r R coordinate of the query point [m].
 * @param phi phi coordinate of the query point [rad].
 * @param z z coordinate of the query point [m].
 * @param bfield The magnetic field data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t BfieldLinear3D_eval_b_db(
    real b_db[12], real r, real phi, real z, BfieldLinear3D *bfield);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(bfield)
/**
 * Evaluate the magnetic axis (R, z) coordinates.
 *
 * Returns the position stored in the struct.
 *
 * @param rz Evaluated axis coordinates [m].
 * @param bfield The magnetic field data.
 *
 * @return Zero if the evaluation succeeded.
 */
err_t BfieldLinear3D_eval_axisrz(real rz[2], BfieldLinear3D *bfield);
DECLARE_TARGET_END

#endif
