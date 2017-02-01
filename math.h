/**
 * @file math.h
 * @brief Header file for math.c
 */
#ifndef MATH_H
#define MATH_H
#include "ascot5.h"

/** @brief Pi */
#define math_pi 3.1415926535897932384626

/** @brief Maximum recursion depth for the simpson integral rule */
#define math_maxSimpsonDepth 20

/** @brief Calculate dot product a[3] dot b[3] */
#define math_dot(a,b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])

/** @brief Calculate cross product c[3] = a[3] cross b[3] */
#define math_cross(a,b,c) (c[0]=a[1]*b[2]-a[2]*b[1],c[1]=a[2]*b[0]-a[0]*b[2],c[2] = a[0]*b[1] - a[1]*b[0])

/** @brief Calculate norm |a[3]| */
#define math_norm(a) (sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]))

/** @brief Calculate norm from components a1, a2, a3 */
#define math_normc(a1, a2, a3) (sqrt(a1*a1+a2*a2+a3*a3))

#pragma omp declare target
#pragma omp declare simd
void math_xyz2rpz(real* xyz, real* rpz);
#pragma omp declare simd
void math_rpz2xyz(real* rpz, real* xyz);
#pragma omp declare simd
void math_vec_rpz2xyz(real* rpz, real* xyz, real phi);
#pragma omp declare simd
void math_vec_xyz2rpz(real* rpz, real* xyz, real phi);
#pragma omp declare simd
real math_normal_rand();
int math_ipow(int a, int p);
double math_simpson_helper(double (*f)(double), double a, double b, double eps, double S, double fa, double fb, double fc, int bottom);
double math_simpson(double (*f)(double), double a, double b, double epsilon);
void math_linspace(real* vec, real a, real b, int n);
#pragma omp end declare target


#endif
