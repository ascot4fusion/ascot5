/**
 * @file spline1Dcomp.h
 * @brief Header file for spline1Dcomp.c
 */
#ifndef SPLINE1DCOMP_H
#define SPLINE1DCOMP_H
#include "../ascot5.h"

#pragma omp declare target
void spline1Dcomp(real* f, int n, int bc, real* c);
#pragma omp end declare target
#endif
