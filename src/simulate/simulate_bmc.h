/**
 * @file simulate_bmc.h
 * @brief Header file for simulate_bmc.c
 */
#ifndef SIMULATE_BMC_H
#define SIMULATE_BMC_H

#include <stdlib.h>
#include "../ascot5.h"
#include "../simulate.h"
#include "../bmc_mesh.h"

#if HERMITE_KNOTS == 5
    /** @brief Roots of the (phycisist's) Hermite polynomial multiplied with
     *  sqrt(2) */
    #define HERMITE_K {-2.8569700138728056541623, -1.3556261799742658658305, \
                   0.0, 1.3556261799742658658305,  2.8569700138728056541623}
    /** @brief Weights of the Gauss-Hermite quadrature divided by sqrt(pi) */
    #define HERMITE_W {0.0112574113277206889333702,0.2220759220056126443999631,\
             8.0/15.0, 0.2220759220056126443999631,0.0112574113277206889333702}
#elif HERMITE_KNOTS == 10
    #define HERMITE_K { -4.859462828332313, -3.581823483551928, \
    -2.484325841638955, -1.465989094391159, -0.484935707515498, \
     0.484935707515498,  1.465989094391159,  2.484325841638955, \
     3.581823483551928,  4.859462828332313}
    #define HERMITE_W {   4.310652630718501e-6, 7.580709343120864e-4, \
    1.911158050076969e-2, 1.354837029802594e-1, 3.446423349320045e-1, \
    3.446423349320045e-1, 1.354837029802594e-1, 1.911158050076969e-2, \
    7.580709343120864e-4, 4.310652630718501e-6}
#endif

void simulate_bmc_gc(
    sim_data* sim, bmc_mesh* mesh, real h, real time, size_t start, size_t stop,
    real* r, real* phi, real* z, real* ppara, real* pperp, int* fate);

#endif