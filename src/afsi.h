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
    int type;        /**< Distribution type (1:beam, 2:thermal)               */
    histogram* beam; /**< Distribution data                                   */
    int thermal;     /**< Thermal species index                               */
} afsi_data;

void afsi_run(sim_data* sim, Reaction reaction, int n,
              afsi_data* react1, afsi_data* react2, real mult,
              histogram* prod1, histogram* prod2);
#endif
