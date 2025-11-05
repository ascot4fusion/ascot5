/**
 * @file math.h
 * Mathematical utility functions.
 */
#ifndef MATH_H
#define MATH_H

#include "defines.h"
#include "parallel.h"
#include <math.h>

/**
 * Find the bin index on a uniform grid.
 */
#define math_bin_index(x, nx, xmin, xmax) floor(nx *(x - xmin) / (xmax - xmin))

/**
 * Copies elements of vector b to vector a.
 */
#define math_copy(a, b)                                                        \
    do                                                                         \
    {                                                                          \
        a[0] = b[0];                                                           \
        a[1] = b[1];                                                           \
        a[2] = b[2];                                                           \
    } while (0)

/**
 * Copies elements of matrix b to matrix a.
 */
#define math_matcopy(a, b)                                                     \
    do                                                                         \
    {                                                                          \
        a[0] = b[0];                                                           \
        a[1] = b[1];                                                           \
        a[2] = b[2];                                                           \
        a[3] = b[3];                                                           \
        a[4] = b[4];                                                           \
        a[5] = b[5];                                                           \
        a[6] = b[6];                                                           \
        a[7] = b[7];                                                           \
        a[8] = b[8];                                                           \
    } while (0)

/**
 * Calculate dot product a[3] dot b[3].
 */
#define math_dot(a, b) (a[0] * b[0] + a[1] * b[1] + a[2] * b[2])

/**
 * Calculate cross product for 3D vectors c = a x b.
 */
#define math_cross(a, b, c)                                                    \
    do                                                                         \
    {                                                                          \
        c[0] = a[1] * b[2] - a[2] * b[1];                                      \
        c[1] = a[2] * b[0] - a[0] * b[2];                                      \
        c[2] = a[0] * b[1] - a[1] * b[0];                                      \
    } while (0)

/**
 * The mixed or triple product of three vectors (a cross b) dot c.
 *
 * mathematica:
 * a={a0,a1,a2};b={b0,b1,b2};c={c0,c1,c2};FullSimplify[Dot[Cross[a,b],c]]
 * -(a2 b1 c0) + a1 b2 c0 + a2 b0 c1 - a0 b2 c1 - a1 b0 c2 + a0 b1 c2
 */
#define math_scalar_triple_product(a, b, c)                                    \
    (-(a[2] * b[1] * c[0]) + (a[1] * b[2] * c[0]) + (a[2] * b[0] * c[1]) -     \
     (a[0] * b[2] * c[1]) - (a[1] * b[0] * c[2]) + (a[0] * b[1] * c[2]))

/**
 * Elementwise vector sum a = a + b.
 */
#define math_sumew(a, b)                                                       \
    do                                                                         \
    {                                                                          \
        a[0] = a[0] + b[0];                                                    \
        a[1] = a[1] + b[1];                                                    \
        a[2] = a[2] + b[2]                                                     \
    } while (0)

/**
 * Multiply vector elements with a scalar a = b*a.
 */
#define math_prod(a, b)                                                        \
    do                                                                         \
    {                                                                          \
        a[0] = a[0] * b;                                                       \
        a[1] = a[1] * b;                                                       \
        a[2] = a[2] * b;                                                       \
    } while (0)

/**
 * Elementwise matrix sum a = a + b
 */
#define math_matsumew(a, b)                                                    \
    do                                                                         \
    {                                                                          \
        a[0] = a[0] + b[0];                                                    \
        a[1] = a[1] + b[1];                                                    \
        a[2] = a[2] + b[2];                                                    \
        a[3] = a[3] + b[3];                                                    \
        a[4] = a[4] + b[4];                                                    \
        a[5] = a[5] + b[5];                                                    \
        a[6] = a[6] + b[6];                                                    \
        a[7] = a[7] + b[7];                                                    \
        a[8] = a[8] + b[8];                                                    \
    } while (0)

/**
 * Multiply matrix a with scalar b: a = b*a.
 */
#define math_matprod(a, b)                                                     \
    do                                                                         \
    {                                                                          \
        a[0] = a[0] * b;                                                       \
        a[1] = a[1] * b;                                                       \
        a[2] = a[2] * b;                                                       \
        a[3] = a[3] * b;                                                       \
        a[4] = a[4] * b;                                                       \
        a[5] = a[5] * b;                                                       \
        a[6] = a[6] * b;                                                       \
        a[7] = a[7] * b;                                                       \
        a[8] = a[8] * b                                                        \
    } while (0)

/**
 * Calculate norm of 3D vector a.
 */
#define math_norm(a) (sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]))

/**
 * Calculate norm of 3D vector from its components a1, a2, a3.
 */
#define math_normc(a1, a2, a3) (sqrt(a1 * a1 + a2 * a2 + a3 * a3))

/**
 * Calculate unit vector b from a 3D vector a.
 */
#define math_unit(a, b)                                                        \
    do                                                                         \
    {                                                                          \
        real _n = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);               \
        b[0] = a[0] / _n;                                                      \
        b[1] = a[1] / _n;                                                      \
        b[2] = a[2] / _n;                                                      \
    } while (0)

/**
 * Convert cartesian coordinates xyz to cylindrical coordinates rpz.
 */
#define math_xyz2rpz(xyz, rpz)                                                 \
    do                                                                         \
    {                                                                          \
        rpz[0] = sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1]);                      \
        rpz[1] = atan2(xyz[1], xyz[0]);                                        \
        rpz[2] = xyz[2];                                                       \
    } while (0)

/**
 * Convert cylindrical coordinates rpz to cartesian coordinates xyz.
 */
#define math_rpz2xyz(rpz, xyz)                                                 \
    do                                                                         \
    {                                                                          \
        xyz[0] = rpz[0] * cos(rpz[1]);                                         \
        xyz[1] = rpz[0] * sin(rpz[1]);                                         \
        xyz[2] = rpz[2];                                                       \
    } while (0)

/**
 * Transform vector from cylindrical to cartesian basis: vrpz -> vxyz.
 * phi is the toroidal angle in radians.
 */
#define math_vec_rpz2xyz(vrpz, vxyz, phi)                                      \
    do                                                                         \
    {                                                                          \
        vxyz[0] = vrpz[0] * cos(phi) - vrpz[1] * sin(phi);                     \
        vxyz[1] = vrpz[0] * sin(phi) + vrpz[1] * cos(phi);                     \
        vxyz[2] = vrpz[2];                                                     \
    } while (0)

/**
 * Transform vector from cartesian to cylindrical basis: vxyz -> vrpz.
 * phi is the toroidal angle in radians.
 */
#define math_vec_xyz2rpz(vxyz, vrpz, phi)                                      \
    do                                                                         \
    {                                                                          \
        vrpz[0] = vxyz[0] * cos(phi) + vxyz[1] * sin(phi);                     \
        vrpz[1] = -vxyz[0] * sin(phi) + vxyz[1] * cos(phi);                    \
        vrpz[2] = vxyz[2];                                                     \
    } while (0)

/**
 * Direct expansion of 3x3 matrix determinant.
 */
#define math_determinant3x3(x1, x2, x3, y1, y2, y3, z1, z2, z3)                \
    (x1) * ((y2) * (z3) - (y3) * (z2)) + (x2) * ((y3) * (z1) - (y1) * (z3)) +  \
        (x3) * ((y1) * (z2) - (y2) * (z1))

DECLARE_TARGET
/**
 * Compute the modulus of two real numbers.
 *
 * @param x The dividend
 * @param y The divisor
 * @return The modulus (remainder) of x and y
 */
real fmod(real x, real y);
DECLARE_TARGET_END

DECLARE_TARGET_SIMD
/**
 * Convert a Jacobian from cylindrical to cartesian coordinates.
 *
 * This function converts a Jacobian located at angle phi and radius r from
 * cylindrical to cartesian coordinates. The input Jacobian is an array with
 *
 * [Br dBr/dr dBr/dphi dBr/dz  Bphi dBphi/dr dBphi/dphi dBphi/dz  Bz dBz/dr
 * dBz/dphi dBz/dz]
 *
 * an the output is
 *
 * [Bx dBx/dx dBx/dy dBx/dz  By dBy/dx dBy/dy dBy/dz  Bz dBz/dx dBz/dy dBz/dz].
 *
 * @param rpz input rpz coordinates in a 3-length array
 * @param xyz output xyz coordinates in a 3-length array
 * @param r   r coordinate of the gradient
 * @param phi phi coordinate of the gradient
 */
void math_jac_rpz2xyz(real *rpz, real *xyz, real r, real phi);

GPU_DECLARE_TARGET_SIMD
/**
 * Convert a Jacobian from cartesian to cylindrical coordinates.
 *
 * This function converts a Jacobian located at angle phi and radius r from
 * cartesian to cylindrical coordinates. The input Jacobian is an array with
 *
 * [Bx dBx/dx dBx/dy dBx/dz  By dBy/dx dBy/dy dBy/dz  Bz dBz/dx dBz/dy dBz/dz]
 *
 * an the output is
 *
 * [Br dBr/dr dBr/dphi dBr/dz  Bphi dBphi/dr dBphi/dphi dBphi/dz  Bz dBz/dr
 * dBz/dphi dBz/dz].
 *
 * @param xyz Output xyz coordinates in a 3-length array.
 * @param rpz Input rpz coordinates in a 3-length array.
 * @param r   r coordinate of the gradient.
 * @param phi phi coordinate of the gradient.
 */
void math_cart2cyl_gradient(
    real *gcyl, const real *gcart, const real *bcart, real r, real phi);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD
/**
 * Matrix multiplication.
 *
 * This function performs matrix multiplication. The matrices are given as
 * arrays, the indexing being row-wise.
 *
 * @param matA Input array representing a d1 x d2 matrix.
 * @param matB Input array representing a d2 x d3 matrix.
 * @param d1   Input scalar dimension (columns in matA and matC).
 * @param d2   Input scalar dimension (rows in matA and columns in matB).
 * @param d3   Input scalar dimension (rows in matB and matC).
 * @param matC Output array representing a d1 x d3 matrix.
 */
inline void
math_matrix_multiplication(real *C, real *A, real *B, const size_t dim[3])
{
    for (size_t i = 0; i < dim[0]; ++i)
    {
        for (size_t j = 0; j < dim[2]; ++j)
        {
            real sum = 0.0;
            for (size_t k = 0; k < dim[1]; ++k)
                sum += A[k * dim[0] + i] * B[j * dim[1] + k];
            C[i * dim[2] + j] = sum;
        }
    }
}
DECLARE_TARGET_END

void test_matrix_multiplication(real *C, real *A, real *B, const size_t dim[3]);

/**
 * Generate normally distributed random numbers.
 *
 * This function generates normally distributed random numbers with expected
 * value of 0 and variance of 1 using the polar method [1]. The rand48 random
 * number generator should be initialized with srand48() before calling this
 * function.
 *
 * [1] G. Marsaglia, T. A. Bray. A convenient method for generating normal
 *     variables. Siam Review 6.3 (1964): 260-264.
 *     http://epubs.siam.org/doi/abs/10.1137/1006063:w
 *
 * @todo Currently only one of the generated numbers is returned; a second
 *       independent variate is Y = v2 * sqrt(-2*log(s) / s).
 * @todo Try other random number generators such as those in Intel MKL.
 */
DECLARE_TARGET_SIMD
real math_normal_rand();
DECLARE_TARGET_SIMD

/**
 * Adaptive Simpsons rule for integral.
 *
 * This function uses recursive splitting of the interval until either maximum
 * number of intervals is reached or given relative error tolerance is obtained.
 *
 * @param f pointer to the integrand function (of one variable).
 * @param a lower limit for the interval.
 * @param b upper limit for the interval.
 * @param eps absolute tolerance.
 *
 * @bug Probably does not work with SIMD
 */
double math_simpson(double (*f)(double), double a, double b, double epsilon);

/**
 * Generate linearly spaced vector.
 *
 * Generates linearly space vector with n elements whose end points are a and b.
 * If n = 1, return b.
 *
 * @param vec Length n array where result is stored.
 * @param a Start point.
 * @param b End point.
 * @param n Number of elements.
 */
void math_linspace(real *vec, real a, real b, int n);

/**
 * Find unique numbers and their frequency in given array.
 *
 * All argument arrays are assumed to have the same length as the array holding
 * the input numbers. If not all numbers are unique, the remaining entries in
 * unique and count arguments are zero.
 *
 * @param in Input array of length n.
 * @param unique Array where unique numbers are stored, must be of length n.
 * @param count Count[i] is how many times unique[i] appears in array in,
 *        must be of length n.
 * @param n Length of the argument arrays.
 */
DECLARE_TARGET_SIMD
int math_point_on_plane(real q[3], real t1[3], real t2[3], real t3[3]);

/**
 * Find if a point is on a given plane.
 *
 * Let t1=(x1, y1, z1), t2=(x2, y2, z2), and t3=(x3, y3, z3) be non-collinear
 * points on the plane. Let q =(x,  y,  z) be an arbitrary point. It lies in
 * the plane if:
 *
 * det(  |x -x1  y -y1 z -z1|
 *       |x2-x1  y2-y1 z2-z1|
 *       |x3-x1  y3-y1 z3-z1| ) == 0
 *
 * det(  |x -x1  y -y1 z -z1|
 *       |x -x2  y -y2 z -z2|
 *       |x -x3  y -y3 z -z3| ) == 0
 *
 * @param q xyz coordinates of a query point.
 * @param t1 xyz coordinates of the first point defining the plane.
 * @param t2 xyz coordinates of the second point defining the plane.
 * @param t3 xyz coordinates of the third point defining the plane.
 */
DECLARE_TARGET_SIMD
void math_barycentric_coords_triangle(
    real AP[3], real AB[3], real AC[3], real n[3], real *s, real *t);

/**
 * Find barycentric coordinates for a given point.
 *
 * Let a,b,c be the vertices of a triangle and p an arbitrary point.
 * The s,t are the barycentric coordinates of p if:
 * p = (1-s-t)*a + s*b + t*c
 *
 * @param AP point position vector [x,y,z].
 * @param AB triangle edge AB vector.
 * @param AC triangle edge AC vector.
 * @param n normal vector [nx, ny, nz] of the triangle.
 * @param s evaluated s coordinate.
 * @param t evaluated t coordinate.
 */
void math_uniquecount(int *in, int *unique, int *count, int n);

/**
 * Search for array element preceding a key value.
 *
 * This function takes an array of real values and a key value, and finds the
 * array element which precedes the key value. The array to be searched has to
 * be in ascending sorted order (according to comparison function rcomp).
 *
 * @param key Real value that serves as key for the search.
 * @param base Pointer to the first object of the array to be searched.
 * @param num Number of elements in array.
 *
 * @return Pointer to array element preceding the key, or NULL if search failed.
 */
DECLARE_TARGET_SIMD
real *math_rsearch(const real key, const real *base, int num);

/**
 * Check if coordinates are within polygon.
 *
 * This function checks if the given coordinates are within a 2D polygon using
 * a modified axis crossing method [1]. Origin is moved to the coordinates and
 * the number of wall segments crossing the positive r-axis are calculated. If
 * this is odd, the point is inside the polygon.
 *
 * [1] D.G. Alciatore, R. Miranda. A Winding Number and Point-in-Polygon
 *     Algorithm. Technical report, Colorado State University, 1995.
 *     http://www.engr.colostate.edu/~dga/dga/papers/point_in_polygon.pdf
 *
 * @param r r coordinate.
 * @param z z coordinate.
 * @param rv Array of r points for the polygon.
 * @param zv Array of z points for the polygon.
 * @param n Number of points in the polygon.
 */
DECLARE_TARGET_SIMD_UNIFORM(rv, zv, n)
int math_point_in_polygon(real r, real z, real *rv, real *zv, int n);

#endif
