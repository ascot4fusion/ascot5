/**
 * Implements math.h.
 */
#include "mathlib.h"
#include "defines.h"
#include <math.h>
#include <stdlib.h>

/**
 * Maximum recursion depth for the simpson integral rule.
 */
#define math_maxSimpsonDepth 20

int rcomp(const void *a, const void *b);
double math_simpson_helper(
    double (*f)(double), double a, double b, double eps, double S, double fa,
    double fb, double fc, int bottom);

real fmod(real x, real y) { return x - y * floor(x / y); }

void math_jac_rpz2xyz(real *rpz, real *xyz, real r, real phi)
{
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
    xyz[1] = temp[0] * c - temp[1] * s / r + (rpz[0] * s + rpz[4] * c) * s / r;
    xyz[2] = temp[0] * s + temp[1] * c / r - (rpz[0] * s + rpz[4] * c) * c / r;
    xyz[3] = temp[2];

    // Step 1: Vector [dBphi/dx dBphi/dy dBphi/dz]
    temp[0] = rpz[1] * s + rpz[5] * c;
    temp[1] = rpz[2] * s + rpz[6] * c;
    temp[2] = rpz[3] * s + rpz[7] * c;

    // Step 2: Gradient
    xyz[5] = temp[0] * c - temp[1] * s / r + (rpz[0] * c - rpz[4] * s) * s / r;
    xyz[6] = temp[0] * s + temp[1] * c / r - (rpz[0] * c - rpz[4] * s) * c / r;
    xyz[7] = temp[2];

    // Step 1: Vector [dBz/dx dBz/dy dBz/dz]
    temp[0] = rpz[9];
    temp[1] = rpz[10];
    temp[2] = rpz[11];

    // Step 2: Gradient
    xyz[9] = temp[0] * c - temp[1] * s / r;
    xyz[10] = temp[0] * s + temp[1] * c / r;
    xyz[11] = temp[2];
}

void math_jac_xyz2rpz(real *xyz, real *rpz, real r, real phi)
{
    // Temporary variables
    real c = cos(phi);
    real s = sin(phi);
    real temp[3];

    // rpz[0] =  xyz[0] * c + xyz[4] * s;
    // rpz[4] = -xyz[0] * s + xyz[4] * c;
    // rpz[8] =  xyz[8];

    // Step 1: Vector [dBx/dr dBx/dphi dBx/dz]
    temp[0] = xyz[1] * c + xyz[5] * s;
    temp[1] = xyz[2] * c + xyz[6] * s;
    temp[2] = xyz[3] * c + xyz[7] * s;

    // Step 2: Gradient
    rpz[1] = temp[0] * c + temp[1] * s;
    rpz[2] = -temp[0] * s * r + temp[1] * c * r - xyz[0] * s + xyz[4] * c;
    rpz[3] = temp[2];

    // Step 1: Vector [dBy/dr dBy/dphi dBy/dz]
    temp[0] = -xyz[1] * s + xyz[5] * c;
    temp[1] = -xyz[2] * s + xyz[6] * c;
    temp[2] = -xyz[3] * s + xyz[7] * c;

    // Step 2: Gradient
    rpz[5] = temp[0] * c + temp[1] * s;
    rpz[6] = -temp[0] * s * r + temp[1] * c * r - xyz[0] * c - xyz[4] * s;
    rpz[7] = temp[2];

    // Step 1: Vector [dBz/dr dBz/dphi dBz/dz]
    temp[0] = xyz[9];
    temp[1] = xyz[10];
    temp[2] = xyz[11];

    // Step 2: Gradient
    rpz[9] = temp[0] * c + temp[1] * s;
    rpz[10] = -temp[0] * s * r + temp[1] * c * r;
    rpz[11] = temp[2];
}

void math_matmul(real *matA, real *matB, int d1, int d2, int d3, real *matC)
{
    real sum;
    for (int i = 0; i < d1; i = i + 1)
    {
        for (int j = 0; j < d3; j = j + 1)
        {
            sum = 0.0;
            for (int k = 0; k < d2; k = k + 1)
            {
                sum = sum + matA[k * d1 + i] * matB[j * d2 + k];
            }
            matC[i * d3 + j] = sum;
        }
    }
}

real math_normal_rand(void)
{
    real X;
    real v1, v2, s;
    do
    {
        v1 = drand48() * 2 - 1;
        v2 = drand48() * 2 - 1;
        s = v1 * v1 + v2 * v2;
    } while (s >= 1);

    X = v1 * sqrt(-2 * log(s) / s);

    return X;
}

int math_ipow(int a, int p)
{
    int pow = 1;
    for (int i = 0; i < p; i++)
    {
        pow *= a;
    }
    return pow;
}

double math_simpson(double (*f)(double), double a, double b, double eps)
{
    double c = (a + b) / 2, h = b - a;
    double fa = f(a), fb = f(b), fc = f(c);
    double S = (h / 6) * (fa + 4 * fc + fb);
    return math_simpson_helper(
        f, a, b, eps, S, fa, fb, fc, math_maxSimpsonDepth);
}

void math_linspace(real *vec, real a, real b, int n)
{
    if (n == 1)
    {
        vec[0] = b;
    }
    else
    {
        real d = (b - a) / (n - 1);
        int i;
        for (i = 0; i < n; i++)
        {
            vec[i] = a + i * d;
        }
    }
}

void math_uniquecount(int *in, int *unique, int *count, int n)
{

    for (int i = 0; i < n; i++)
    {
        unique[i] = 0;
        count[i] = 0;
    }

    int n_unique = 0;
    for (int i = 0; i < n; i++)
    {

        int test = in[i];
        int isunique = 1;
        for (int j = 0; j < n_unique; j++)
        {
            if (test == unique[j])
            {
                isunique = 0;
                count[j] += 1;
            }
        }

        if (isunique)
        {
            unique[n_unique] = test;
            count[n_unique] = 1;
            n_unique++;
        }
    }
}

int math_point_on_plane(real q[3], real t1[3], real t2[3], real t3[3])
{
    real x = q[0], y = q[1], z = q[2];
    real x1 = t1[0], y1 = t1[1], z1 = t1[2];
    real x2 = t2[0], y2 = t2[1], z2 = t2[2];
    real x3 = t3[0], y3 = t3[1], z3 = t3[2];

    int val = 0;
    if (math_determinant3x3(
            x - x1, y - y1, z - z1, x2 - x1, y2 - y1, z2 - z1, x3 - x1, y3 - y1,
            z3 - z1) != 0.0)
    {
        return val;
    }
    if (math_determinant3x3(
            x - x1, y - y1, z - z1, x - x2, y - y2, z - z2, x - x3, y - y3,
            z - z3) != 0.0)
    {
        return val;
    }
    return val;
}

void math_barycentric_coords_triangle(
    real AP[3], real AB[3], real AC[3], real n[3], real *s, real *t)
{
    real n0[3], area;
    math_unit(n, n0);
    area = math_scalar_triple_product(AB, AC, n0);
    *s = math_scalar_triple_product(AP, AC, n0) / area;
    *t = math_scalar_triple_product(AB, AP, n0) / area;
}

/**
 * Helper comparison routine for "math_rsearch".
 *
 * This function checks if a key value is between two consecutive values in a
 * real array. The array to be searched has to be in ascending sorted order.
 *
 * @param a First comparison value (key).
 * @param b First value in array.
 *
 * @return 0 if key is between the consecutive values, 1 if key is greater than
 * the first value, and -1 if key is smaller than the first value.
 */
int rcomp(const void *a, const void *b)
{
    real a_val = *((real *)a);
    real b_val = *((real *)b);
    real c_val = *((real *)b + 1);

    if (a_val >= b_val && a_val < c_val)
    {
        return 0;
    }
    return a_val > b_val ? 1 : -1;
}

real *math_rsearch(const real key, const real *base, int num)
{
    return (real *)bsearch(&key, base, num - 1, sizeof(real), rcomp);
}

int math_point_in_polygon(real r, real z, real *rv, real *zv, int n)
{
    int hits = 0;

    int i;
    for (i = 0; i < n - 1; i++)
    {
        real z1 = zv[i] - z;
        real z2 = zv[i + 1] - z;
        real r1 = rv[i] - r;
        real r2 = rv[i + 1] - r;
        if (z1 * z2 < 0)
        {
            real ri = r1 + (z1 * (r2 - r1)) / (z1 - z2);
            if (ri > 0)
            {
                hits++;
            }
        }
    }
    return hits % 2;
}

/**
 * Helper routine for "math_simpson".
 */
double math_simpson_helper(
    double (*f)(double), double a, double b, double eps, double S, double fa,
    double fb, double fc, int bottom)
{
    double c = (a + b) / 2, h = b - a;
    double d = (a + c) / 2, e = (c + b) / 2;
    double fd = f(d), fe = f(e);
    double Sleft = (h / 12) * (fa + 4 * fd + fc);
    double Sright = (h / 12) * (fc + 4 * fe + fb);
    double S2 = Sleft + Sright;
    if (bottom <= 0 || fabs(S2 - S) <= eps * fabs(S))
    {
        return S2 + (S2 - S) / 15;
    }
    return math_simpson_helper(f, a, c, eps, Sleft, fa, fc, fd, bottom - 1) +
           math_simpson_helper(f, c, b, eps, Sright, fc, fb, fe, bottom - 1);
}
