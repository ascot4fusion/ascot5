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
#include "diag/dist_5D.h"

/**
 * @brief Structure for passing in 2D thermal temperature and density
 */
typedef struct {
    int n_r;          /**< number of r bins       */
    real min_r;       /**< value of lowest r bin  */
    real max_r;       /**< value of highest r bin */

    int n_phi;        /**< number of r bins       */
    real min_phi;     /**< value of lowest r bin  */
    real max_phi;     /**< value of highest r bin */

    int n_z;          /**< number of z bins       */
    real min_z;       /**< value of lowest z bin  */
    real max_z;       /**< value of highest z bin */

    real* temperature; /**< pointer to start of histogram array */
    real* density;     /**< pointer to start of histogram array */
} afsi_thermal_data;

/**
 * @brief Wrapper around input data structures
 */
typedef struct {
    int type;                        /**< Distribution type (1:5D, 2:thermal) */
    dist_5D_data* dist_5D;           /**< Distribution data                   */
    afsi_thermal_data* dist_thermal; /**< Thermal data                        */
} afsi_data;

void afsi_run(sim_data* sim, Reaction reaction, int n,
              afsi_data* react1, afsi_data* react2, real mult,
              dist_5D_data* prod1, dist_5D_data* prod2);
void afsi_test_dist(dist_5D_data* dist1);
void afsi_test_thermal();

#endif
