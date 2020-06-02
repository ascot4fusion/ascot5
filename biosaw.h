/**
 * @file biosaw.h
 * @brief Header file for biosaw.c
 */
#ifndef BIOSAW_H
#define BIOSAW_H

#include "ascot5.h"

void biosaw_calc_B(int n, real* x, real* y, real* z,
                   int coil_n, real* coil_x, real* coil_y, real* coil_z,
                   real* Bx, real* By, real* Bz);

#endif
