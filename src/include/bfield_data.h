/**
 * @file bfield_data.h
 * Magnetic field data structures.
 */
#ifndef BFIELD_DATA_H
#define BFIELD_DATA_H

#include "interp_data.h"

/**
 * Magnetic field types.
 *
 * These are used to direct function calls to correct implementations.
 */
typedef enum Bfield_type
{
    BFIELD_CARTESIAN = 1, /**< Corresponds to BfieldCartesian.                */
    BFIELD_ANALYTICAL,    /**< Corresponds to BfieldAnalytical.               */
    BFIELD_SPLINE2D,      /**< Corresponds to BfieldSpline2D.                 */
    BFIELD_SPLINE3D,      /**< Corresponds to BfieldSpline3D.                 */
    BFIELD_STELLARATOR,   /**< Corresponds to BfieldStellarator.              */
} Bfield_type;

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
    size_t nripple;     /**< Number of toroidal field coils.                  */
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
    real axisrz[2];    /**< Magnetic axis (R, z) coordinates [m].             */
    real psilimits[2]; /**< Poloidal flux at axis and separatrix [Wb/rad].    */
    Spline2D psi;      /**< Poloidal flux spline-interpolant [Wb/rad].        */
    Spline2D br;       /**< Spline-interpolant for R component of B [T].      */
    Spline2D bz;       /**< Spline-interpolant for z component of B [T].      */
    Spline2D bphi;     /**< Spline-interpolant for phi component of B [T].    */
} BfieldSpline2D;

/**
 * Spline-interpolated 3D magnetic field data.
 */
typedef struct
{
    real axisrz[2];    /**< Magnetic axis (R, z) coordinates [m].             */
    real psilimits[2]; /**< Poloidal flux at axis and separatrix [Wb/rad].    */
    Spline2D psi;      /**< Poloidal flux spline-interpolant [Wb/rad].        */
    Spline3D br;       /**< Spline-interpolant for R component of B [T].      */
    Spline3D bz;       /**< Spline-interpolant for z component of B [T].      */
    Spline3D bphi;     /**< Spline-interpolant for phi component of B [T].    */
} BfieldSpline3D;

/**
 * Spline-interpolated stellarator magnetic field data.
 */
typedef struct
{
    real psilimits[2]; /**< Poloidal flux at axis and separatrix [Wb/rad].    */
    Linear1D axisr;    /**< Linear interpolant for axis R coordinate [m].     */
    Linear1D axisz;    /**< Linear interpolant for axis z coordinate [m].     */
    Spline3D psi;      /**< Poloidal flux spline-interpolant [Wb/rad].        */
    Spline3D br;       /**< Spline-interpolant for R component of B [T].      */
    Spline3D bz;       /**< Spline-interpolant for z component of B [T].      */
    Spline3D bphi;     /**< Spline-interpolant for phi component of B [T].    */
} BfieldStellarator;

/**
 * Magnetic field simulation data.
 *
 * The intended usage is that only single type of data is used at a time. This
 * is declared using the ``type`` field.
 */
typedef struct
{
    BfieldCartesian *cartesian;     /**< Cartesian data.                      */
    BfieldAnalytical *analytical;   /**< Analytical data.                     */
    BfieldSpline2D *spline2d;       /**< Spline 2D data.                      */
    BfieldSpline3D *spline3d;       /**< Spline 3D data.                      */
    BfieldStellarator *stellarator; /**< Stellarator data.                    */
    Bfield_type type;               /**< Current magnetic field type.         */
} Bfield;

#endif
