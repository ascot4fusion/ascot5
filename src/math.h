/**
 * @file math.h
 * @brief Header file for math.c
 */
#ifndef MATH_H
#define MATH_H

#include <math.h>
#include "ascot5.h"

/** @brief One degree in radians */
#define math_degrad 0.0174532925199432957692
/** @brief One radian in degrees */
#define math_raddeg 57.295779513082320876798

/** @brief Maximum recursion depth for the simpson integral rule */
#define math_maxSimpsonDepth 20

/** @brief Copies elements of vector b to vector a */
#define math_copy(a,b) do { a[0]=b[0];a[1]=b[1];a[2]=b[2]; } while(0)

/** @brief Copies elements of matrix b to matrix a */
#define math_matcopy(a,b) do { a[0]=b[0];a[1]=b[1];a[2]=b[2];a[3]=b[3]; \
        a[4]=b[4];a[5]=b[5];a[6]=b[6];a[7]=b[7];a[8]=b[8]; } while(0)

/** @brief Calculate dot product a[3] dot b[3] */
#define math_dot(a,b) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])

/** @brief Calculate cross product for 3D vectors c = a x b */
#define math_cross(a,b,c) do { c[0]=a[1]*b[2]-a[2]*b[1];                \
        c[1]=a[2]*b[0]-a[0]*b[2];c[2]=a[0]*b[1]-a[1]*b[0]; } while(0)

 /** @brief the mixed or triple product of three vectors (a cross b) dot c
  * mathematica:
  * a={a0,a1,a2};b={b0,b1,b2};c={c0,c1,c2};FullSimplify[Dot[Cross[a,b],c]]
  * -(a2 b1 c0) + a1 b2 c0 + a2 b0 c1 - a0 b2 c1 - a1 b0 c2 + a0 b1 c2
  * */
#define math_scalar_triple_product(a,b,c) \
       ( - (a[2] * b[1] * c[0] ) \
         + (a[1] * b[2] * c[0] ) \
         + (a[2] * b[0] * c[1] ) \
         - (a[0] * b[2] * c[1] ) \
         - (a[1] * b[0] * c[2] ) \
         + (a[0] * b[1] * c[2] )     )

/** @brief Elementwise vector sum a = a + b */
#define math_sumew(a,b) do { a[0]=a[0]+b[0];a[1]=a[1]+b[1];a[2]=a[2]+b[2] } while(0)

/** @brief Multiply vector elements with a scalar a = b*a */
#define math_prod(a,b) do { a[0]=a[0]*b;a[1]=a[1]*b;a[2]=a[2]*b; } while(0)

/** @brief Elementwise matrix sum a = a + b */
#define math_matsumew(a,b) do { a[0]=a[0]+b[0];a[1]=a[1]+b[1];a[2]=a[2]+b[2]; \
        a[3]=a[3]+b[3];a[4]=a[4]+b[4];a[5]=a[5]+b[5];a[6]=a[6]+b[6];    \
        a[7]=a[7]+b[7];a[8]=a[8]+b[8]; } while(0)

/** @brief Multiply matrix a with scalar b: a = b*a */
#define math_matprod(a,b) do { a[0]=a[0]*b;a[1]=a[1]*b;a[2]=a[2]*b;  \
        a[3]=a[3]*b;a[4]=a[4]*b;a[5]=a[5]*b;a[6]=a[6]*b;a[7]=a[7]*b; \
        a[8]=a[8]*b } while(0)

/** @brief Calculate norm of 3D vector a */
#define math_norm(a) (sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]))

/** @brief Calculate norm of 3D vector from its components a1, a2, a3 */
#define math_normc(a1, a2, a3) (sqrt(a1*a1+a2*a2+a3*a3))

/** @brief Calculate unit vector b from a 3D vector a */
#define math_unit(a, b) do {real _n=sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]); \
        b[0]=a[0]/_n;b[1]=a[1]/_n;b[2]=a[2]/_n; } while(0)

/** @brief Convert cartesian coordinates xyz to cylindrical coordinates rpz */
#define math_xyz2rpz(xyz, rpz) do { rpz[0]=sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]); \
        rpz[1]=atan2(xyz[1],xyz[0]);rpz[2]=xyz[2]; } while(0)

/** @brief Convert cylindrical coordinates rpz to cartesian coordinates xyz */
#define math_rpz2xyz(rpz, xyz) do { xyz[0]=rpz[0]*cos(rpz[1]); \
        xyz[1]=rpz[0]*sin(rpz[1]);xyz[2]=rpz[2]; } while(0)

/** @brief Transform vector from cylindrical to cartesian basis: vrpz -> vxyz,
    phi is the toroidal angle in radians */
#define math_vec_rpz2xyz(vrpz, vxyz, phi) do {     \
        vxyz[0]=vrpz[0]*cos(phi)-vrpz[1]*sin(phi); \
        vxyz[1]=vrpz[0]*sin(phi)+vrpz[1]*cos(phi); \
        vxyz[2]=vrpz[2]; } while(0)

/** @brief Transform vector from cartesian to cylindrical basis: vxyz -> vrpz,
    phi is the toroidal angle in radians */
#define math_vec_xyz2rpz(vxyz, vrpz, phi) do {      \
        vrpz[0]=vxyz[0]*cos(phi)+vxyz[1]*sin(phi);  \
        vrpz[1]=-vxyz[0]*sin(phi)+vxyz[1]*cos(phi); \
        vrpz[2]=vxyz[2]; } while(0)

 /** @brief Direct expansion of 3x3 matrix determinant
  */
#define math_determinant3x3( \
    x1, x2, x3,              \
    y1, y2, y3,              \
    z1, z2, z3)              \
    (x1) * ( (y2)*(z3) - (y3)*(z2) ) + \
    (x2) * ( (y3)*(z1) - (y1)*(z3) ) + \
    (x3) * ( (y1)*(z2) - (y2)*(z1) )


/** @brief Convert degrees to radians */
#define math_deg2rad(a) (a * math_degrad)

/** @brief Convert radians to degrees */
#define math_rad2deg(a) (a * math_raddeg)

#pragma omp declare target
#pragma omp declare simd
void math_jac_rpz2xyz(real* rpz, real* xyz, real r, real phi);
#pragma omp declare simd
void math_jac_xyz2rpz(real* xyz, real* rpz, real r, real phi);
#pragma omp declare simd
void math_matmul(real* matA, real* matB, int d1, int d2, int d3, real* matC);
#pragma omp declare simd
real math_normal_rand();
#pragma omp declare simd
int math_ipow(int a, int p);
double math_simpson(double (*f)(double), double a, double b, double epsilon);
void math_linspace(real* vec, real a, real b, int n);
#pragma omp declare simd
int math_point_on_plane(real q[3], real t1[3], real t2[3], real t3[3]);
#pragma omp declare simd
void math_barycentric_coords_triangle(
    real AP[3], real AB[3], real AC[3], real n[3], real *s, real *t);
void math_uniquecount(int* in, int* unique, int* count, int n);
#pragma omp declare simd
real* math_rsearch(const real key, const real* base, int num);
#pragma omp declare simd uniform(rv,zv,n)
int math_point_in_polygon(real r, real z, real* rv, real* zv, int n);
#pragma omp end declare target


#endif
