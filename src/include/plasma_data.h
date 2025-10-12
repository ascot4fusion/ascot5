/**
 * @file plasma_data.h
 * Plasma data types.
 */
#ifndef PLASMA_DATA_H
#define PLASMA_DATA_H

/**
 * Plasma data types.
 */
typedef enum plasma_type
{
    PLASMA_LINEAR1D = 1, /**< Corresponds to PlasmaLinear1D.                  */
    PLASMA_DYNAMIC1D,    /**< Corresponds to PlasmaDynamic1D.                 */
} Plasma_type;

/**
 * Parameters for linearly-interpolated 1D plasma.
 */
typedef struct
{
    size_t nrho;     /**< Number of rho grid points.                          */
    size_t nspecies; /**< Number of plasma species (including electrons).     */
    int *anum;       /**< Atomic mass number of the ion species.              */
    int *znum;       /**< Charge number of the ion species.                   */
    real *mass;      /**< Plasma species (electrons first) mass [kg].         */
    real *charge;    /**< Plasma species (electrons first) charge [C].        */
    real *rho;       /**< Rho grid values [1].                                */
    real *vtor;      /**< Plasma toroidal rotation [rad/s].                   */

    /**
     * Electron and ion temperatures [J].
     *
     * First ``nrho`` elements are electron temperatures, followed by ``nrho``
     * ion temperatures (same for all ion species).
     */
    real *temp;

    /**
     * Electron and ion densities [m^-3].
     *
     * First ``nrho`` elements are electron densities, followed by ion densities
     * (in same order as they are listed in ``anum``, ``znum``, etc.).
     */
    real *dens;
} PlasmaLinear1D;

/**
 * Parameters for linearly-interpolated time-dependent 1D plasma.
 */
typedef struct
{
    size_t nrho;     /**< Number of rho grid points.                          */
    size_t ntime;    /**< Number of time grid points.                         */
    size_t nspecies; /**< Number of plasma species (including electrons).     */
    int *anum;       /**< Atomic mass number of the ion species.              */
    int *znum;       /**< Charge number of the ion species.                   */
    real *mass;      /**< Plasma species (electrons first) mass [kg].         */
    real *charge;    /**< Plasma species (electrons first) charge [C].        */
    real *rho;       /**< Rho grid values [1].                                */
    real *time;      /**< Time grid values [s].                               */
    real *vtor;      /**< Plasma toroidal rotation [rad/s].                   */

    /**
     * Electron and ion temperatures [J].
     *
     * First ``nrho`` x ``ntime`` elements are electron temperatures, followed
     * by ``nrho`` x ``ntime`` ion temperatures (same for all ion species).
     * Otherwise the layout is (rhoi, tj) = [j*nrho + i] (C order).
     */
    real *temp;

    /**
     * Electron and ion densities [m^-3].
     *
     * First ``nrho`` x ``ntime`` elements are electron densities, followed by
     * ion densities (in same order as they are listed in ``anum``, ``znum``,
     * etc.). Otherwise the layout is (rhoi, tj) = [j*nrho + i] (C order).
     */
    real *dens;
} PlasmaDynamic1D;

/**
 * Plasma simulation data.
 *
 * The intended usage is that only single type of data is used at a time. This
 * is declared using the ``type`` field.
 */
typedef struct
{
    PlasmaLinear1D *linear1d;   /**< Linearly-interpolated 1D plasma.         */
    PlasmaDynamic1D *dynamic1d; /**< Linear time-dependent 1D plasma.         */
    Plasma_type type;           /**< Current plasma type.                     */
} Plasma;

#endif
