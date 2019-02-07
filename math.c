/**
 * @file math.c
 * @brief Mathematical utility functions
 */
#define _XOPEN_SOURCE 500
#include <math.h>
#include <stdlib.h>
#include "ascot5.h"
#include "math.h"

double math_simpson_helper(double (*f)(double), double a, double b, double eps,
                           double S, double fa, double fb, double fc,
                           int bottom);

/**
 * @brief Convert a Jacobian from cylindrical to cartesian coordinates
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
void math_jac_rpz2xyz(real* rpz, real* xyz, real r, real phi) {
    // Temporary variables
    real c = cos(phi);
    real s = sin(phi);
    real temp[3];

    xyz[0] = rpz[0] * c - rpz[4] * s;
    xyz[4] = rpz[0] * s + rpz[4] * c;
    xyz[8] = rpz[8];

    // Step 1: Vector [dBr/dx dBr/dy dBr/dz]
    temp[0] = rpz[1] * c - rpz[5] * s;
    temp[1] = rpz[2] * c - rpz[6] * s;
    temp[2] = rpz[3] * c - rpz[7] * s;

    // Step 2: Gradient
    xyz[1] = temp[0] * c - temp[1] * s / r + (rpz[0]*s + rpz[4]*c) * s / r;
    xyz[2] = temp[0] * s + temp[1] * c / r - (rpz[0]*s + rpz[4]*c) * c / r;
    xyz[3] = temp[2];

    // Step 1: Vector [dBphi/dx dBphi/dy dBphi/dz]
    temp[0] = rpz[1] * s + rpz[5] * c;
    temp[1] = rpz[2] * s + rpz[6] * c;
    temp[2] = rpz[3] * s + rpz[7] * c;

    // Step 2: Gradient
    xyz[5] = temp[0] * c - temp[1] * s / r + (rpz[0]*c - rpz[4]*s) * s / r;
    xyz[6] = temp[0] * s + temp[1] * c / r - (rpz[0]*c - rpz[4]*s) * c / r;
    xyz[7] = temp[2];

    // Step 1: Vector [dBz/dx dBz/dy dBz/dz]
    temp[0] = rpz[9];
    temp[1] = rpz[10];
    temp[2] = rpz[11];

    // Step 2: Gradient
    xyz[9]  = temp[0] * c - temp[1] * s / r;
    xyz[10] = temp[0] * s + temp[1] * c / r;
    xyz[11] = temp[2];
}

/**
 * @brief Convert a Jacobian from cartesian to cylindrical coordinates
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
 * @param xyz output xyz coordinates in a 3-length array
 * @param rpz input rpz coordinates in a 3-length array
 * @param r   r coordinate of the gradient
 * @param phi phi coordinate of the gradient
 */
void math_jac_xyz2rpz(real* xyz, real* rpz, real r, real phi) {
    // Temporary variables
    real c = cos(phi);
    real s = sin(phi);
    real temp[3];

    rpz[0] =  xyz[0] * c + xyz[4] * s;
    rpz[4] = -xyz[0] * s + xyz[4] * c;
    rpz[8] =  xyz[8];

    // Step 1: Vector [dBx/dr dBx/dphi dBx/dz]
    temp[0] = xyz[1] * c + xyz[5] * s;
    temp[1] = xyz[2] * c + xyz[6] * s;
    temp[2] = xyz[3] * c + xyz[7] * s;

    // Step 2: Gradient
    rpz[1] =  temp[0] * c + temp[1] * s;
    rpz[2] = -temp[0] * s * r + temp[1] * c * r
        - xyz[0] * s + xyz[4] * c;
    rpz[3] =  temp[2];

    // Step 1: Vector [dBy/dr dBy/dphi dBy/dz]
    temp[0] = -xyz[1] * s + xyz[5] * c;
    temp[1] = -xyz[2] * s + xyz[6] * c;
    temp[2] = -xyz[3] * s + xyz[7] * c;

    // Step 2: Gradient
    rpz[5] =  temp[0] * c + temp[1] * s;
    rpz[6] = -temp[0] * s * r + temp[1] * c * r
        - xyz[0] * c - xyz[4] * s;
    rpz[7] =  temp[2];

    // Step 1: Vector [dBz/dr dBz/dphi dBz/dz]
    temp[0] = xyz[9];
    temp[1] = xyz[10];
    temp[2] = xyz[11];

    // Step 2: Gradient
    rpz[9]  =  temp[0] * c + temp[1] * s;
    rpz[10] = -temp[0] * s * r + temp[1] * c * r;
    rpz[11] =  temp[2];
}

/**
 * @brief Matrix multiplication
 *
 * This function performs matrix multiplication. The matrices are given as
 * arrays, the indexing being row-wise.
 *
 * @param matA input array representing a d1 x d2 matrix
 * @param matB input array representing a d2 x d3 matrix
 * @param d1   input scalar dimension (columns in matA and matC)
 * @param d2   input scalar dimension (rows in matA and columns in matB)
 * @param d3   input scalar dimension (rows in matB and matC)
 * @param matC output array representing a d1 x d3 matrix
 */
void math_matmul(real* matA, real* matB, int d1, int d2, int d3, real* matC) {
    real sum;
    for (int i = 0; i < d1; i=i+1) {
        for (int j = 0; j < d3; j=j+1) {
            sum = 0.0;
            for (int k = 0; k < d2; k=k+1){
                sum = sum + matA[k * d1 + i]*matB[j * d2 + k];
            }
            matC[i * d3 + j] = sum;
        }
    }
}

/**
 * @brief Generate normally distributed random numbers
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
 *       independent variate is Y = v2 * sqrt(-2*log(s) / s)
 * @todo Try other random number generators such as those in Intel MKL
 */
real math_normal_rand(void) {
    real X;
    real v1, v2, s;
    do {
        v1 = drand48() * 2 - 1;
        v2 = drand48() * 2 - 1;
        s = v1*v1 + v2*v2;
    } while(s >= 1);

    X = v1 * sqrt(-2*log(s) / s);

    return X;
}

/**
 * @brief Calculate a^p where both a and p are integers (p >= 0)
 *
 * @param a argument
 * @param p power
 *
 * @return a^b
 */
int math_ipow(int a, int p) {
    int pow = 1;
    for(int i = 0; i < p; i++) {
        pow *= a;
    }
    return pow;
}

/**
 * @brief Adaptive Simpsons rule for integral
 *
 * This function uses recursive splitting of the interval
 * until either maximum number of intervals is reached or
 * given relative error tolerance is obtained
 *
 * @param f: pointer to the integrand function (of one variable)
 * @param a: lower limit for the interval
 * @param b: upper limit for the interval
 * @param eps: absolute tolerance
 * @param maxDepth: maximum number of splits for the interval
 *
 * @bug Probably does not work with SIMD
 */
double math_simpson(double (*f)(double), double a, double b, double eps) {
    double c = (a + b)/2, h = b - a;
    double fa = f(a), fb = f(b), fc = f(c);
    double S = (h/6)*(fa + 4*fc + fb);
    return math_simpson_helper(f, a, b, eps, S, fa, fb, fc,
                               math_maxSimpsonDepth);
}

/**
 * @brief Generate linearly spaced vector
 *
 * Generates linearly space vector with n elements whose end points are a and b.
 * If n = 1, return b.
 *
 * @param vec n length array where result is stored
 * @param a start point
 * @param b end point
 * @param n number of elements
 */
void math_linspace(real* vec, real a, real b, int n) {
    if(n == 1) {
        vec[0] = b;
    } else {
        real d = (b-a)/(n-1);
        int i;
        for(i = 0; i < n; i++) {
            vec[i] = a+i*d;
        }
    }
}

/**
 * @brief Find unique numbers and their frequency in given array
 *
 * All argument arrays are assumed to have the same length as the array holding
 * the input numbers. If not all numbers are unique, the remaining entries in
 * unique and count arguments are zero.
 *
 * @param in input array of length n
 * @param unique array where unique numbers are stored, must be of length n
 * @param count count[i] is how many times unique[i] appears in array in,
 *        must be of length n
 * @param n length of the argument arrays
 */
void math_uniquecount(int* in, int* unique, int* count, int n) {

    for(int i=0; i<n; i++) {
        unique[i] = 0;
        count[i] = 0;
    }

    int n_unique = 0;
    for(int i=0; i<n; i++) {

        int test = in[i];
        int isunique = 1;
        for(int j=0; j<n_unique; j++) {
            if(test == unique[j]) {
                isunique = 0;
                count[j] += 1;
            }
        }

        if(isunique) {
            unique[n_unique] = test;
            count[n_unique] = 1;
            n_unique++;
        }
    }
}

/**
 * @brief Helper routine for "math_simpson"
 */
double math_simpson_helper(double (*f)(double), double a, double b, double eps,
                           double S, double fa, double fb, double fc,
                           int bottom) {
    double c = (a + b)/2, h = b - a;
    double d = (a + c)/2, e = (c + b)/2;
    double fd = f(d), fe = f(e);
    double Sleft = (h/12)*(fa + 4*fd + fc);
    double Sright = (h/12)*(fc + 4*fe + fb);
    double S2 = Sleft + Sright;
    if (bottom <= 0 || fabs(S2 - S) <= eps*fabs(S)) {
        return  S2 + (S2 - S)/15;
    }
    return math_simpson_helper(f, a, c, eps, Sleft,  fa, fc, fd, bottom-1)
        +math_simpson_helper(f, c, b, eps, Sright, fc, fb, fe, bottom-1);
}

int rcomp(const void* a, const void* b) {
    real a_val = *((real*) a);
    real b_val = *((real*) b);
    real c_val = *((real*) b + 1);

    if (a_val >= b_val && a_val < c_val) {
        return 0;
    }
    return a_val > b_val ? 1 : -1;
}

real* math_rsearch(const real key, const real* base, int num) {
    return (real*) bsearch(&key, base, num-1, sizeof(real), rcomp);
}

/**
 * @brief Check if coordinates are within polygon
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
 * @param r r coordinate
 * @param z z coordinate
 * @param rv array of r points for the polygon
 * @param zv array of z points for the polygon
 * @param n number of points in the polygon
 */
int math_point_in_polygon(real r, real z, real* rv, real* zv, int n) {
    int hits = 0;

    int i;
    for(i = 0; i < n - 1; i++) {
        real z1 = zv[i] - z;
        real z2 = zv[i+1] - z;
        real r1 = rv[i] - r;
        real r2 = rv[i+1] - r;
        if(z1 * z2 < 0) {
            real ri = r1 + (z1*(r2-r1)) / (z1-z2);
            if(ri > 0) {
                hits++;
            }
        }
    }
    return hits % 2;
}

