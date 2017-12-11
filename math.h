/**
 * @file math.h
 * @brief Header file for math.c
 */
#ifndef MATH_H
#define MATH_H
#include "ascot5.h"

/** @brief Pi */
#define math_pi 3.1415926535897932384626
/** @brief One degree in radians */
#define math_degrad 0.0174532925199432957692
/** @brief One radian in degrees */
#define math_raddeg 57.295779513082320876798

/** @brief Maximum recursion depth for the simpson integral rule */
#define math_maxSimpsonDepth 20

/** @brief Copies elements of vector b to vector a */
#define math_copy(a,b) (a[0]=b[0],a[1]=b[1],a[2]=b[2])

/** @brief Copies elements of matrix b to matrix a */
#define math_matcopy(a,b) (a[0]=b[0],a[1]=b[1],a[2]=b[2],a[3]=b[3],a[4]=b[4],a[5]=b[5],a[6]=b[6],a[7]=b[7],a[8]=b[8])

/** @brief Calculate dot product a[3] dot b[3] */
#define math_dot(a,b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])

/** @brief Calculate cross product c[3] = a[3] cross b[3] */
#define math_cross(a,b,c) (c[0]=a[1]*b[2]-a[2]*b[1],c[1]=a[2]*b[0]-a[0]*b[2],c[2] = a[0]*b[1] - a[1]*b[0])

/** @brief Elementwise vector sum a = a + b */
#define math_sumew(a,b) (a[0]=a[0]+b[0],a[1]=a[1]+b[1],a[2]=a[2]+b[2])

/** @brief Multiply vector elements with a scalar a = b*a */
#define math_prod(a,b) (a[0]=a[0]*b,a[1]=a[1]*b,a[2]=a[2]*b)

/** @brief Elementwise matrix sum a = a + b */
#define math_matsumew(a,b) (a[0]=a[0]+b[0],a[1]=a[1]+b[1],a[2]=a[2]+b[2],a[3]=a[3]+b[3],a[4]=a[4]+b[4],a[5]=a[5]+b[5],a[6]=a[6]+b[6],a[7]=a[7]+b[7],a[8]=a[8]+b[8])

/** @brief Multiply matrix elements with a scalar a = b*a */
#define math_matprod(a,b) (a[0]=a[0]*b,a[1]=a[1]*b,a[2]=a[2]*b,a[3]=a[3]*b,a[4]=a[4]*b,a[5]=a[5]*b,a[6]=a[6]*b,a[7]=a[7]*b,a[8]=a[8]*b)

/** @brief Calculate norm |a[3]| */
#define math_norm(a) (sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]))

/** @brief Calculate norm from components a1, a2, a3 */
#define math_normc(a1, a2, a3) (sqrt(a1*a1+a2*a2+a3*a3))

/** @brief Convert degrees to radians */
#define math_deg2rad(a) (a * math_degrad)

/** @brief Convert radians to degrees */
#define math_rad2deg(a) (a * math_raddeg)

#pragma omp declare target
#pragma omp declare simd simdlen(8)
void math_unit(real* vec, real* vec_unit);
#pragma omp declare simd simdlen(8)
void math_xyz2rpz(real* xyz, real* rpz);
#pragma omp declare simd simdlen(8)
void math_rpz2xyz(real* rpz, real* xyz);
#pragma omp declare simd simdlen(8)
void math_vec_rpz2xyz(real* rpz, real* xyz, real phi);
#pragma omp declare simd simdlen(8)
void math_vec_xyz2rpz(real* xyz, real* rpz, real phi);
#pragma omp declare simd simdlen(8)
void math_grad_xyz2rpz(real* xyz, real* rpz, real r, real phi);
#pragma omp declare simd simdlen(8)
void math_grad_rpz2xyz(real* xyz, real* rpz, real r, real phi);
#pragma omp declare simd simdlen(8)
void math_jac_rpz2xyz(real* rpz, real* xyz, real r, real phi);
#pragma omp declare simd simdlen(8)
void math_jac_xyz2rpz(real* xyz, real* rpz, real r, real phi);
#pragma omp declare simd simdlen(8)
void math_matmul(real* matA, real* matB, int d1, int d2, int d3, real* matC);
#pragma omp declare simd simdlen(8)
real math_normal_rand();
int math_ipow(int a, int p);
double math_simpson_helper(double (*f)(double), double a, double b, double eps, double S, double fa, double fb, double fc, int bottom);
double math_simpson(double (*f)(double), double a, double b, double epsilon);
void math_linspace(real* vec, real a, real b, int n);
#pragma omp end declare target


#endif
