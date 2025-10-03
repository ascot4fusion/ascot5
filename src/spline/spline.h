/**
 * @file spline.h
 * @brief Header file for splineexpl.c and splinecomp.c
 */
#ifndef SPLINE_H
#define SPLINE_H
#include "ascot5.h"

void splineexpl(real* f, int n, int bc, real* c);
void splinecomp(real* f, int n, int bc, real* c);
#endif
