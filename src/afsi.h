/**
 * @file afsi.c
 * @brief ASCOT Fusion Source Integrator AFSI
 */
#ifndef AFSI_H
#define AFSI_H

#include <math.h>
#include "ascot5.h"
#include "consts.h"
#include "simulate.h"
#include "random.h"
#include "boschhale.h"
#include "diag/hist.h"

/**
 * @brief Valid momentum space basis.
 */
typedef enum {
    PPARPPERP,
    EKINXI
} mom_space_basis;

/**
 * @brief Wrapper around input data structures
 */
typedef struct {
    int type1;          /**< Distribution type (1:beam, 2:thermal)            */
    int type2;          /**< Distribution type (1:beam, 2:thermal)            */
    int thermal1;       /**< Thermal species index for reactant 1             */
    int thermal2;       /**< Thermal species index for reactant 2             */
    histogram* beam1;   /**< Distribution data for reactant 1                 */
    histogram* beam2;   /**< Distribution data for reactant 2                 */
    real* r;            /**< Radial coordinate at the grid center [m]         */
    real* phi;          /**< Toroidal coordinate at the grid center [rad]     */
    real* z;            /**< Axial coordinate at the grid center [m]          */
    real* vol;          /**< Grid cell volume [m^3]                           */
    size_t volshape[3]; /**< Dimensions of r, phi, z, and volume              */
    Reaction reaction;  /**< The fusion reaction that is modelled             */
    real mult;          /**< Multiplication factor which is 0.5 if species is
                             interacting with itself, 1.0 otherwise           */
} afsi_data;

void afsi_run(sim_data* sim, afsi_data* data, int n,
              histogram* prod1, histogram* prod2);

void afsi_run_new(sim_data* sim, afsi_data* afsi, int n,
            real* rvec, real* phivec, real* zvec, real* prod2);
#endif


// key files: afsi.c (contains afsi_run), afsi.h (makes in communication afsi.c with libascot.so and therefore the python side of the program)
// afsi.py thermal (initializes the output histograms prod1 and prod2 and runs afsi_run)

// we could either write a new afsi_run function in afsi.c -> add new afsi_run in afsi.h -> new function in afsi.py that calls new afsi_run -> recompile libascot.so, recompile ascot2py.py or push and pull the modified ascot2py
// OR modify existing afsi_run function in afsi.c -> no needed to modify afsi.h -> proper modify afsi.py to create the histograms prod1 and prod2 but necessary to get distribution first -> only needed to recompile libascot.so (not sure about this)