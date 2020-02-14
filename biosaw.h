/**
 * @file biosaw.h
 * @brief Header file for biosaw.c
 */
#ifndef BIOSAW_H
#define BIOSAW_H

#include "ascot5.h"

void biosaw_calc_B(real* x, int coil_n, real* coil_x, real* coil_y,
                   real* coil_z, real* B);

#endif
