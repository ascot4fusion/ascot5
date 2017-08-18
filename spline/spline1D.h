/**
 * @file spline1D.h
 * @brief Header file for spline1D.c
 */
#ifndef SPLINE1D_H
#define SPLINE1D_H
#include "../ascot5.h"

#pragma omp declare target
void spline1D(real* f, int n, int bc, real* c);
#pragma omp end declare target
#endif
