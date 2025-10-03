/**
 * @file biosaw.h
 * Tool to evaluate magnetic field from coil geometry with Biot-Savart law.
 */
#ifndef BIOSAW_H
#define BIOSAW_H

#include "ascot5.h"

/**
 * Evaluate magnetic field due to a coil at given points.
 *
 * The magnetic field is evaluated using Biot-Savart law.
 *
 * @param n Number of query points.
 * @param coil_n number of points in coil geometry
 * @param x x-coordinate of a query point [m]
 * @param y y-coordinate of a query point [m]
 * @param z z-coordinate of a query point [m]
 * @param coil_x coil geometry x-coordinate [m]
 * @param coil_y coil geometry y-coordinate [m]
 * @param coil_z coil geometry z-coordinate [m]
 * @param Bx evaluated magnetic field x-component [T]
 * @param By evaluated magnetic field y-component [T]
 * @param Bz evaluated magnetic field z-component [T]
 */
void biosaw_calc_B(
    int n, int coil_n, real *x, real *y, real *z, real *coil_x, real *coil_y,
    real *coil_z, real *bx, real *by, real *bz);

#endif
