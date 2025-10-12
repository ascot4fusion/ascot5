/**
 * @file efield_data.h
 * Electric field data types.
 */
#ifndef EFIELD_DATA_H
#define EFIELD_DATA_H

#include "bfield_data.h"
#include "interp_data.h"

/**
 * Electric field types.
 *
 * These are used to direct function calls to correct implementations.
 */
typedef enum Efield_type
{
    EFIELD_CARTESIAN = 1, /**< Corresponds to EfieldCartesian.                */
    EFIELD_POTENTIAL1D,   /**< Corresponds to EfieldPotential1D.              */
} Efield_type;

/**
 * Cartesian electric field data.
 */
typedef struct
{
    real exyz[3]; /**< Electric field vector in Cartesian basis [V/m].        */
} EfieldCartesian;

/**
 * Radial electric field potential data.
 */
typedef struct
{
    Linear1D dvdrho; /**< Gradient of the potential with respect to rho [V].  */
} EfieldPotential1D;

/**
 * Electric field simulation data.
 *
 * The intended usage is that only single type of data is used at a time. This
 * is declared using the ``type`` field.
 */
typedef struct
{
    EfieldCartesian *cartesian;     /**< Cartesian data.                      */
    EfieldPotential1D *potential1d; /**< Radial potential data.               */
    Efield_type type;               /**< Current electric field type.         */
} Efield;

#endif
