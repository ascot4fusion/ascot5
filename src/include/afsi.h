/**
 * @file afsi.h
 * ASCOT Fusion Source Integrator AFSI.
 */
#ifndef AFSI_H
#define AFSI_H

#include "ascot5.h"
#include "boschhale.h"
#include "consts.h"
#include "diag.h"
#include "random.h"
#include "simulate.h"
#include <math.h>

/**
 * Valid momentum space basis.
 */
typedef enum
{
    PPARPPERP,
    EKINXI
} mom_space_basis;

/**
 * Wrapper around input data structures.
 */
typedef struct
{
    int type1;          /**< Distribution type (1:beam, 2:thermal).           */
    int type2;          /**< Distribution type (1:beam, 2:thermal).           */
    int thermal1;       /**< Thermal species index for reactant 1.            */
    int thermal2;       /**< Thermal species index for reactant 2.            */
    histogram *beam1;   /**< Distribution data for reactant 1.                */
    histogram *beam2;   /**< Distribution data for reactant 2.                */
    real *r;            /**< Radial coordinate at the grid center [m].        */
    real *phi;          /**< Toroidal coordinate at the grid center [rad].    */
    real *z;            /**< Axial coordinate at the grid center [m].         */
    real *vol;          /**< Grid cell volume [m^3].                          */
    size_t volshape[3]; /**< Dimensions of r, phi, z, and volume.             */
    Reaction reaction;  /**< The fusion reaction that is modelled.            */
    real mult;          /**< Multiplication factor which is 0.5 if species is
                             interacting with itself, 1.0 otherwise.          */
} afsi_data;

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
void afsi_run(
    sim_data *sim, afsi_data *data, int n, histogram *prod1, histogram *prod2);
#endif
