/**
 * @file B_field.h
 * Magnetic field interface
 *
 * This is an interface through which magnetic field data is initialized and
 * accessed.
 *
 * To add a new magnetic field instance, make sure these functions are
 * implemented and called from this interface, and that B_field.h contains enum
 * type for the new instance.
 *
 * The interface checks which instance given data corresponds to from
 * B_field_offload_data.type and B_field_data.type from the structure that is
 * given as an argument, and calls the relevant function for that instance.
 */
#ifndef B_FIELD_H
#define B_FIELD_H

#include "ascot5.h"
#include "error.h"
#include "interp.h"
#include "linint.h"
#include "offload.h"

/**
 * @brief Magnetic field types
 *
 * Magnetic field types are used in the magnetic field interface to direct
 * function calls to correct magnetic field instances. Each magnetic field
 * instance must have a corresponding type.
 */
typedef enum B_field_type
{
    B_field_type_cartesian,   /**< Cartesian field                            */
    B_field_type_analytical,  /**< Analytical field                           */
    B_field_type_spline2d,    /**< Spline 2D field                            */
    B_field_type_spline3d,    /**< Spline 3D field                            */
    B_field_type_stellarator, /**< Stellarator field                          */
} B_field_type;

/**
 * Cartesian magnetic field data.
 */
typedef struct
{
    real axisrz[2]; /**< Pseudo magnetic axis (R, z) coordinates [m].         */
    real psival;    /**< Pseudo poloidal flux value [Wb/rad].                 */
    real rhoval;    /**< Pseudo normalized poloidal flux value [1].           */
    real bxyz[3];   /**< Magnetic field at origo: [bx, by, bz] [T].           */

    /** Magnetic field Jacobian [T/m].
     *
     * Layout: [dB_x/dx, dB_x/dy, dB_x/dz, dB_y/dx, dB_y/dy, dB_y/dz, dB_z/dx,
     * dB_z/dy, dB_z/dz].
     */
    real jacobian[9];
} BfieldCartesian;

/**
 * Analytic magnetic field data.
 */
typedef struct
{
    int nripple;        /**< Number of toroidal field coils.                  */
    real bphi;          /**< Toroidal field strength at ``rmajor`` [T].       */
    real rmajor;        /**< Plasma major radius [m].                         */
    real rminor;        /**< Plasma minor radius [m].                         */
    real psiscaling;    /**< Scaling factor to scale poloidal flux [Wb/rad].  */
    real ripplescaling; /**< Ripple scaling factor.                           */
    real rippledamping; /**< Ripple minor radius damping factor [m].          */
    real axisrz[2];     /**< Magnetic axis (R, z) coordinates [m].            */
    real psilimits[2];  /**< Poloidal flux at axis and separatrix [Wb/rad].   */

    /**
     * Coefficients for evaluating poloidal flux.
     * Layout: [c1, c2, ..., c12, A].
     */
    real coefficients[13];
} BfieldAnalytical;

/**
 * Spline-interpolated 2D magnetic field data.
 */
typedef struct
{
    real axisrz[2];     /**< Magnetic axis (R, z) coordinates [m].            */
    real psilimits[2];  /**< Poloidal flux at axis and separatrix [Wb/rad].   */
    interp2D_data psi;  /**< Poloidal flux spline-interpolant [Wb/rad].       */
    interp2D_data br;   /**< Spline-interpolant for R component of B [T].     */
    interp2D_data bz;   /**< Spline-interpolant for z component of B [T].     */
    interp2D_data bphi; /**< Spline-interpolant for phi component of B [T].   */
} BfieldSpline2D;

/**
 * Spline-interpolated 3D magnetic field data.
 */
typedef struct
{
    real axisrz[2];     /**< Magnetic axis (R, z) coordinates [m].            */
    real psilimits[2];  /**< Poloidal flux at axis and separatrix [Wb/rad].   */
    interp2D_data psi;  /**< Poloidal flux spline-interpolant [Wb/rad].       */
    interp3D_data br;   /**< Spline-interpolant for R component of B [T].     */
    interp3D_data bz;   /**< Spline-interpolant for z component of B [T].     */
    interp3D_data bphi; /**< Spline-interpolant for phi component of B [T].   */
} BfieldSpline3D;

/**
 * Spline-interpolated stellarator magnetic field data.
 */
typedef struct
{
    real psilimits[2];   /**< Poloidal flux at axis and separatrix [Wb/rad].  */
    linint1D_data axisr; /**< Linear interpolant for axis R coordinate [m].   */
    linint1D_data axisz; /**< Linear interpolant for axis z coordinate [m].   */
    interp3D_data psi;   /**< Poloidal flux spline-interpolant [Wb/rad].      */
    interp3D_data br;    /**< Spline-interpolant for R component of B [T].    */
    interp3D_data bz;    /**< Spline-interpolant for z component of B [T].    */
    interp3D_data bphi;  /**< Spline-interpolant for phi component of B [T].  */
} BfieldStellarator;

/**
 * Magnetic field simulation data.
 *
 * The intended usage is that only single type of data is used at a time. This
 * is declared using the `type` field.
 */
typedef struct
{
    BfieldCartesian *cartesian;     /**< Cartesian data.                      */
    BfieldAnalytical *analytical;   /**< Analytical data.                     */
    BfieldSpline2D *spline2d;       /**< Spline 2D data.                      */
    BfieldSpline3D *spline3d;       /**< Spline 3D data.                      */
    BfieldStellarator *stellarator; /**< Stellarator data.                    */
    B_field_type type;              /**< Current magnetic field type.         */
} B_field_data;

/**
 * Free allocated resources
 *
 * Spline-interpolants are freed.
 *
 * @param bfield The struct whose fields are deallocated.
 */
void B_field_free(B_field_data *bfield);

/**
 * Offload data to the accelerator.
 *
 * @param bfield The struct to offload.
 */
void B_field_offload(B_field_data *bfield);

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
a5err B_field_eval_psi(
    real psi[1], real r, real phi, real z, real t, B_field_data *bfield);
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
a5err B_field_eval_psi_dpsi(
    real psi_dpsi[4], real r, real phi, real z, real t, B_field_data *bfield);
DECLARE_TARGET_END

/**
 * Evaluate normalized poloidal flux and its derivative with respect to poloidal
 * flux.
 *
 * @param rho Evaluated normalized poloidal flux and it's derivatives [1].
 *        Layout: [rho, drho/dpsi].
 * @param r R coordinate of the query point [m].
 * @param phi phi coordinate of the query point [rad].
 * @param z z coordinate of the query point [m].
 * @param bfield The magnetic field data.
 *
 * @return Zero if the evaluation succeeded.
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(bfield)
a5err B_field_eval_rho(real rho[2], real psi, B_field_data *bfield);
DECLARE_TARGET_END

/**
 * Evaluate normalized poloidal flux and its derivatives.
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
a5err B_field_eval_rho_drho(
    real rho_drho[4], real r, real phi, real z, B_field_data *bfield);
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
a5err B_field_eval_B(
    real b[3], real r, real phi, real z, real t, B_field_data *bfield);
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
a5err B_field_eval_B_dB(
    real b_db[15], real r, real phi, real z, real t, B_field_data *bfield);
DECLARE_TARGET_END

/**
 * Evaluate the magnetic axis (R, z) coordinates.
 *
 * Returns the position stored in the struct.
 *
 * @param axisrz Evaluated axis coordinates [m].
 * @param bfield The magnetic field data.
 *
 * @return Zero if the evaluation succeeded.
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(bfield)
a5err B_field_get_axis_rz(real rz[2], B_field_data *bfield, real phi);
DECLARE_TARGET_END

#endif
