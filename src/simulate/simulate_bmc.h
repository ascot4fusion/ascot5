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
    /** @brief Roots of the (phycisist's) Hermite polunomial multiplied with
     *  sqrt(2) */
    #define HERMITE_K {-2.8569700138728056541623, -1.3556261799742658658305, \
        0.0, 1.3556261799742658658305, 2.8569700138728056541623}
    /** @brief Weights of the Gauss-Hermite quadrature divided by sqrt(pi) */
    #define HERMITE_W {0.0112574113277206889333702,0.2220759220056126443999631,\
        8.0/15.0, 0.2220759220056126443999631, 0.0112574113277206889333702}
#elif HERMITE_KNOTS == 10
    #define HERMITE_K {-4.859463, -3.581823, -2.484326, -1.465989, -0.484936, \
                       0.484936, 1.465989, 2.484326, 3.581823, 4.859463}
    #define HERMITE_W {0.000011, 0.001900, 0.047906, 0.339607, 0.863890, \
                       0.863890, 0.339607, 0.047906, 0.001900, 0.000011}
#endif

void simulate_bmc_gc(
    sim_data* sim, bmc_mesh* mesh, real h, real mass, real charge,
    int anum, int znum, real time, size_t start, size_t stop,
    real* r, real* phi, real* z, real* ppara, real* pperp, int* fate);

#endif