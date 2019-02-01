/**
 * @file spline.h
 * @brief Header file for spline_expl.c and spline_comp.c
 */
#ifndef SPLINE_H
#define SPLINE_H
#include "../ascot5.h"

#pragma omp declare target
void spline(real* f, int n, int bc, real* c);
void splinecomp(real* f, int n, int bc, real* c);
#pragma omp end declare target
#endif
