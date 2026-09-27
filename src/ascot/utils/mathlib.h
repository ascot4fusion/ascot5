/**
 * @file mathlib.h
 * Mathematical utility functions.
 */
#ifndef MATH_H
#define MATH_H

#include "defines.h"
#include "parallel.h"
#include <math.h>

/**
 * Find the bin index on a uniform grid.
 *
 * Undefined behavior if x < xmin or x > xmax.
 *
 * @param x Value to find the bin index for.
 * @param nx Number of bins.
 * @param xmin Minimum value of the grid.
 * @param xmax Maximum value of the grid.
 * @return The bin index.
 */
#define math_bin_index(x, nx, xmin, xmax)                                      \
    (((x) == (xmax)) ? ((nx) - 1)                                              \
                     : floor((nx) * ((x) - (xmin)) / ((xmax) - (xmin))))

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
void math_linspace(real *vec, real a, real b, size_t n);

/**
 * Matrix multiplication.
 *
 * @param a 3x3 matrix.
 * @param b 3x3 matrix.
 * @param c Result of a times b.
 */
#define math_matmul(a, b, c)                                                   \
    do                                                                         \
    {                                                                          \
        (c)[0] = (a)[0] * (b)[0] + (a)[1] * (b)[3] + (a)[2] * (b)[6];          \
        (c)[1] = (a)[0] * (b)[1] + (a)[1] * (b)[4] + (a)[2] * (b)[7];          \
        (c)[2] = (a)[0] * (b)[2] + (a)[1] * (b)[5] + (a)[2] * (b)[8];          \
        (c)[3] = (a)[3] * (b)[0] + (a)[4] * (b)[3] + (a)[5] * (b)[6];          \
        (c)[4] = (a)[3] * (b)[1] + (a)[4] * (b)[4] + (a)[5] * (b)[7];          \
        (c)[5] = (a)[3] * (b)[2] + (a)[4] * (b)[5] + (a)[5] * (b)[8];          \
        (c)[6] = (a)[6] * (b)[0] + (a)[7] * (b)[3] + (a)[8] * (b)[6];          \
        (c)[7] = (a)[6] * (b)[1] + (a)[7] * (b)[4] + (a)[8] * (b)[7];          \
        (c)[8] = (a)[6] * (b)[2] + (a)[7] * (b)[5] + (a)[8] * (b)[8];          \
    } while (0)

/**
 * Matrix vector multiplication.
 *
 * @param a 3x3 matrix.
 * @param b 3D vector.
 * @param c Result of a times b.
 */
#define math_matvecmul(a, b, c)                                                \
    do                                                                         \
    {                                                                          \
        (c)[0] = (b)[0] * (a)[0] + (b)[1] * (a)[1] + (b)[2] * (a)[2];          \
        (c)[1] = (b)[3] * (a)[0] + (b)[4] * (a)[1] + (b)[5] * (a)[2];          \
        (c)[2] = (b)[6] * (a)[0] + (b)[7] * (a)[1] + (b)[8] * (a)[2];          \
    } while (0)

/**
 * Calculate dot product of two 3D vectors.
 *
 * @param a First vector.
 * @param b Second vector.
 * @return The dot product of a and b.
 */
#define math_dot(a, b) ((a)[0] * (b)[0] + (a)[1] * (b)[1] + (a)[2] * (b)[2])

/**
 * Calculate cross product of two 3D vectors.
 *
 * @param a First vector.
 * @param b Second vector.
 * @param c The cross product of a and b.
 */
#define math_cross(a, b, c)                                                    \
    do                                                                         \
    {                                                                          \
        (c)[0] = (a)[1] * (b)[2] - (a)[2] * (b)[1];                            \
        (c)[1] = (a)[2] * (b)[0] - (a)[0] * (b)[2];                            \
        (c)[2] = (a)[0] * (b)[1] - (a)[1] * (b)[0];                            \
    } while (0)

/**
 * Calculate the triple product of three 3D vectors.
 *
 * @param a First vector.
 * @param b Second vector.
 * @param c Third vector.
 * @return The triple product of a x b dot c.
 */
#define math_scalar_triple_product(a, b, c)                                    \
    (-((a)[2] * (b)[1] * (c)[0]) + ((a)[1] * (b)[2] * (c)[0]) +                \
     ((a)[2] * (b)[0] * (c)[1]) - ((a)[0] * (b)[2] * (c)[1]) -                 \
     ((a)[1] * (b)[0] * (c)[2]) + ((a)[0] * (b)[1] * (c)[2]))

/**
 * Calculate norm of a 3D vector.
 *
 * @param a Vector.
 * @return The norm of a.
 */
#define math_norm(a) (sqrt((a)[0] * (a)[0] + (a)[1] * (a)[1] + (a)[2] * (a)[2]))

/**
 * Calculate norm of 3D vector from its components.
 *
 * @param a1 First component.
 * @param a2 Second component.
 * @param a3 Third component.
 * @return The norm of the vector.
 */
#define math_normc(a1, a2, a3) (sqrt((a1) * (a1) + (a2) * (a2) + (a3) * (a3)))

/**
 * Calculate unit vector b from a 3D vector a.
 *
 * @param a Input vector.
 * @param b Output unit vector.
 */
#define math_unit(a, b)                                                        \
    do                                                                         \
    {                                                                          \
        real _n = sqrt((a)[0] * (a)[0] + (a)[1] * (a)[1] + (a)[2] * (a)[2]);   \
        (b)[0] = (a)[0] / _n;                                                  \
        (b)[1] = (a)[1] / _n;                                                  \
        (b)[2] = (a)[2] / _n;                                                  \
    } while (0)

/**
 * Convert cartesian coordinates to cylindrical coordinates.
 *
 * The angle is in radians.
 *
 * @param xyz Input cartesian coordinates [x, y, z].
 * @param rpz Output cylindrical coordinates [r, phi, z].
 */
#define math_xyz2rpz(xyz, rpz)                                                 \
    do                                                                         \
    {                                                                          \
        (rpz)[0] = sqrt((xyz)[0] * (xyz)[0] + (xyz)[1] * (xyz)[1]);            \
        (rpz)[1] = atan2((xyz)[1], (xyz)[0]);                                  \
        (rpz)[2] = (xyz)[2];                                                   \
    } while (0)

/**
 * Convert cylindrical coordinates to cartesian coordinates.
 *
 * The angle is in radians.
 *
 * @param rpz Input cylindrical coordinates [r, phi, z].
 * @param xyz Output cartesian coordinates [x, y, z].
 */
#define math_rpz2xyz(rpz, xyz)                                                 \
    do                                                                         \
    {                                                                          \
        (xyz)[0] = (rpz)[0] * cos((rpz)[1]);                                   \
        (xyz)[1] = (rpz)[0] * sin((rpz)[1]);                                   \
        (xyz)[2] = (rpz)[2];                                                   \
    } while (0)

/**
 * Transform vector from cylindrical to cartesian basis.
 *
 * @param vrpz Input vector in cylindrical coordinates [vr, vphi, vz].
 * @param vxyz Output vector in cartesian coordinates [vx, vy, vz].
 * @param phi Toroidal angle [rad].
 */
#define math_vec_rpz2xyz(vrpz, vxyz, phi)                                      \
    do                                                                         \
    {                                                                          \
        (vxyz)[0] = (vrpz)[0] * cos((phi)) - (vrpz)[1] * sin((phi));           \
        (vxyz)[1] = (vrpz)[0] * sin((phi)) + (vrpz)[1] * cos((phi));           \
        (vxyz)[2] = (vrpz)[2];                                                 \
    } while (0)

/**
 * Transform vector from cartesian to cylindrical basis.
 *
 * @param vxyz Input vector in cartesian coordinates [vx, vy, vz].
 * @param vrpz Output vector in cylindrical coordinates [vr, vphi, vz].
 * @param phi Toroidal angle [rad].
 */
#define math_vec_xyz2rpz(vxyz, vrpz, phi)                                      \
    do                                                                         \
    {                                                                          \
        (vrpz)[0] = (vxyz)[0] * cos((phi)) + (vxyz)[1] * sin((phi));           \
        (vrpz)[1] = -(vxyz)[0] * sin((phi)) + (vxyz)[1] * cos((phi));          \
        (vrpz)[2] = (vxyz)[2];                                                 \
    } while (0)

/**
 * Calculate determinant of a 3x3 matrix (element wise).
 *
 * @param x11 First element of first row.
 * @param x12 Second element of first row.
 * @param x13 Third element of first row.
 * @param x21 First element of second row.
 * @param x22 Second element of second row.
 * @param x23 Third element of second row.
 * @param x31 First element of third row.
 * @param x32 Second element of third row.
 * @param x33 Third element of third row.
 * @return Determinant of matrix.
 */
#define math_determinant3x3(x11, x12, x13, x21, x22, x23, x31, x32, x33)       \
    ((x11) * ((x22) * (x33) - (x23) * (x32)) -                                 \
     (x12) * ((x21) * (x33) - (x23) * (x31)) +                                 \
     (x13) * ((x21) * (x32) - (x22) * (x31)))

/**
 * Compute the modulus of two real numbers.
 *
 * @param x The dividend.
 * @param y The divisor.
 * @return The modulus (remainder) of x and y.
 */
#define fmod(x, y) ((x) - (y) * floor((x) / (y)))

/**
 * Return absolute value of integer.
 *
 * @note This function exists because 'abs' did not work for AMD GPUs.
 * @param x Integer.
 * @return Absolute value of integer.
 */
DECLARE_TARGET
static inline int math_iabs(int x) { return x < 0 ? -x : x; }
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD
/**
 * Convert field vector and Jacobian from cylindrical to cartesian coordinates.
 *
 * @param a_daxyz Output Jacobian and vector (in cartesian coordinates). Has
 *        format [Ax, Ay, Az, dAx/dx, dAx/dy, dAx/dz, dAy/dx, dAy/dy, dAy/dz,
 *        dAz/dx, dAz/dy, dAz/dz].
 * @param a_darpz Input Jacobian and vector (in cylindrical coordinates). Has
 *        format [Ar, Aphi, Az, dAr/dr, dAr/dphi, dAr/dz, dAphi/dr, dAphi/dphi,
 *        dAphi/dz, dAz/dr, dAz/dphi, dAz/dz].
 * @param r Radial coordinate.
 * @param phi Angular coordinate [rad].
 */
void math_jac_rpz2xyz(
    real a_daxyz[12], const real a_darpz[12], real r, real phi);

GPU_DECLARE_TARGET_SIMD
/**
 * Convert field vector and Jacobian from cartesian to cylindrical coordinates.
 *
 * @param a_darpz Output Jacobian and vector (in cylindrical coordinates). Has
 *        format [Ar, Aphi, Az, dAr/dr, dAr/dphi, dAr/dz, dAphi/dr, dAphi/dphi,
 *        dAphi/dz, dAz/dr, dAz/dphi, dAz/dz].
 * @param a_daxyz Input Jacobian and vector (in cartesian coordinates). Has
 *        format [Ax, Ay, Az, dAx/dx, dAx/dy, dAx/dz, dAy/dx, dAy/dy, dAy/dz,
 *        dAz/dx, dAz/dy, dAz/dz].
 * @param r Radial coordinate.
 * @param phi Angular coordinate [rad].
 */
void math_jac_xyz2rpz(
    real a_darpz[12], const real a_daxyz[12], real r, real phi);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(n, xv, yv)
/**
 * Check if coordinates are within polygon.
 *
 * This function checks if the given coordinates are within a 2D polygon using
 * a modified axis crossing method [1]. The edge cases where the point is
 * exactly at the vertex or on the edge of the polygon are not handled and can
 * return either 1 or 0.
 *
 * [1] D.G. Alciatore, R. Miranda. A Winding Number and Point-in-Polygon
 *     Algorithm. Technical report, Colorado State University, 1995.
 *     http://www.engr.colostate.edu/~dga/dga/papers/point_in_polygon.pdf
 *
 * @param n Number of vertices in the polygon.
 * @param xv Vertex x coordinates.
 * @param yv Vertex y coordinates.
 * @param x Tested x coordinate.
 * @param y Tested y coordinate.
 * @return 1 if inside, 0 otherwise.
 */
int math_point_in_polygon(
    size_t n, const real xv[n], const real yv[n], real x, real y);

GPU_DECLARE_TARGET_SIMD_UNIFORM(gamma)
/**
 * Check if a given angle is between two angles.
 *
 * This function checks if a point, defined by gamma, on a unit circle is on
 * the (shorter) arc defined by two angles alpha and beta.
 *
 * @param alpha The first angle [rad].
 * @param beta The second angle [rad].
 * @param gamma The tested angle assumed to be normalized between (0, 2*pi)
 *        [rad].
 *
 * @return -1 if not between, 0 <= k <= 1 otherwise where k is defined such that
 *          gamma = mod(alpha + k * (beta - alpha), 2*pi).
 */
real math_crossed_plane(real alpha, real beta, real gamma);

/**
 * Evaluate vector operations.
 *
 * This function is used to evaluate the vector operations in the code for
 * testing.
 *
 * @param a First vector.
 * @param b Second vector.
 * @param c Third vector.
 * @param dot Result of a dot b.
 * @param cross Result of a x b.
 * @param triple Result of a x b dot c.
 * @param det Result of determinant of matrix [a, b, c].
 * @param norm Result of |a|.
 * @param normc Result of |a| calculated from a[0], a[1], a[2] explicitly.
 * @param unit Result of unit vector of a.
 */
void math_test_eval_vector_operations(
    const real a[3], const real b[3], const real c[3], real dot[1],
    real cross[3], real triple[1], real det[1], real norm[1], real normc[1],
    real unit[3]);

/**
 * Evaluate vector coordinate transformations.
 *
 * This function is used to evaluate the vector transformations in the code for
 * testing.
 *
 * @param xyz Input xyz position.
 * @param rpz Input rpz position.
 * @param vxyz Input vector in cartesian basis at position rpz.
 * @param vrpz Input vector in cylindrical basis at position rpz.
 * @param xyz_out Argument rpz converted to xyz.
 * @param rpz_out Argument xyz converted to rpz.
 * @param vxyz_out Argument vrpz converted to cartesian basis.
 * @param vrpz_out Argument vxyz converted to cylindrical basis.
 */
void math_test_eval_vector_transformations(
    const real xyz[3], const real rpz[3], const real vxyz[3],
    const real vrpz[3], real xyz_out[3], real rpz_out[3], real vxyz_out[3],
    real vrpz_out[3]);

/**
 * Find the bin index on a uniform grid.
 *
 * This function is for testing.
 *
 * @param x Value to find the bin index for.
 * @param nx Number of bins.
 * @param xmin Minimum value of the grid.
 * @param xmax Maximum value of the grid.
 * @param bin_index The bin index.
 */
void math_test_eval_bin_index(
    const size_t nx, const real xmin, const real xmax, const real x,
    size_t bin_index[1]);

/**
 * Evaluate modulus of two real numbers.
 *
 * This function is for testing.
 *
 * @param a Dividend.
 * @param b Divisor.
 * @param out Result.
 */
void math_test_eval_fmod(real a, real b, real out[1]);

/**
 * Evaluate absolute value of an integer.
 *
 * This function is for testing.
 *
 * @param a Integer.
 * @param b Result.
 */
void math_test_eval_iabs(int a, int b[1]);

#endif
