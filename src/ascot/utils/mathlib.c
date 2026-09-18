/**
 * Implements math.h.
 */
#include "mathlib.h"
#include "consts.h"
#include "defines.h"
#include <math.h>
#include <stdlib.h>

void math_jac_rpz2xyz(
    real a_daxyz[12], const real a_darpz[12], real r, real phi)
{
    const real c = cos(phi);
    const real s = sin(phi);
    const real ir = 1.0 / r;
    const real cc = c * c;
    const real ss = s * s;
    const real sc = s * c;

    a_daxyz[0] = a_darpz[0] * c - a_darpz[1] * s;
    a_daxyz[1] = a_darpz[0] * s + a_darpz[1] * c;
    a_daxyz[2] = a_darpz[2];
    a_daxyz[3] = a_darpz[3] * cc - a_darpz[6] * sc +
                 ir * (-a_darpz[4] * sc + a_darpz[7] * ss + a_darpz[0] * ss +
                       a_darpz[1] * sc);
    a_daxyz[4] = a_darpz[3] * sc - a_darpz[6] * ss +
                 ir * (a_darpz[4] * cc - a_darpz[7] * sc - a_darpz[0] * sc -
                       a_darpz[1] * cc);
    a_daxyz[5] = a_darpz[5] * c - a_darpz[8] * s;
    a_daxyz[6] = a_darpz[3] * sc + a_darpz[6] * cc +
                 ir * (-a_darpz[4] * ss - a_darpz[7] * sc - a_darpz[0] * sc +
                       a_darpz[1] * ss);
    a_daxyz[7] = a_darpz[3] * ss + a_darpz[6] * sc +
                 ir * (a_darpz[4] * sc + a_darpz[7] * cc + a_darpz[0] * cc -
                       a_darpz[1] * sc);
    a_daxyz[8] = a_darpz[5] * s + a_darpz[8] * c;
    a_daxyz[9] = a_darpz[9] * c - a_darpz[10] * s * ir;
    a_daxyz[10] = a_darpz[9] * s + a_darpz[10] * c * ir;
    a_daxyz[11] = a_darpz[11];
}

void math_jac_xyz2rpz(
    real a_darpz[12], const real a_daxyz[12], real r, real phi)
{
    const real c = cos(phi);
    const real s = sin(phi);
    const real cc = c * c;
    const real ss = s * s;
    const real sc = s * c;

    a_darpz[0] = a_daxyz[0] * c + a_daxyz[1] * s;
    a_darpz[1] = -a_daxyz[0] * s + a_daxyz[1] * c;
    a_darpz[2] = a_daxyz[2];
    a_darpz[3] =
        a_daxyz[3] * cc + (a_daxyz[4] + a_daxyz[6]) * sc + a_daxyz[7] * ss;
    a_darpz[6] =
        -a_daxyz[3] * sc + a_daxyz[6] * cc - a_daxyz[4] * ss + a_daxyz[7] * sc;
    a_darpz[9] = a_daxyz[9] * c + a_daxyz[10] * s;
    a_darpz[4] = r * (-a_daxyz[3] * sc + a_daxyz[4] * cc - a_daxyz[6] * ss +
                      a_daxyz[7] * sc) -
                 a_daxyz[0] * s + a_daxyz[1] * c;
    a_darpz[7] = r * (a_daxyz[3] * ss - a_daxyz[4] * sc - a_daxyz[6] * sc +
                      a_daxyz[7] * cc) -
                 a_daxyz[0] * c - a_daxyz[1] * s;
    a_darpz[10] = r * (-a_daxyz[9] * s + a_daxyz[10] * c);
    a_darpz[5] = a_daxyz[5] * c + a_daxyz[8] * s;
    a_darpz[8] = -a_daxyz[5] * s + a_daxyz[8] * c;
    a_darpz[11] = a_daxyz[11];
}

void math_linspace(real *vec, real a, real b, size_t n)
{
    if (n == 1)
    {
        vec[0] = b;
    }
    else
    {
        real d = (b - a) / (n - 1);
        for (size_t i = 0; i < n; i++)
            vec[i] = a + i * d;
    }
}

int math_point_in_polygon(
    size_t n, const real xv[n], const real yv[n], real x, real y)
{
    int winding = 0;

    for (size_t i = 0; i < n; ++i)
    {
        size_t j = (i + 1) % n;

        real y1 = yv[i] - y;
        real y2 = yv[j] - y;
        real x1 = xv[i] - x;
        real x2 = xv[j] - x;

        if (y1 <= 0 && y2 > 0)
        {
            real xi = x1 - y1 * (x2 - x1) / (y2 - y1);
            winding += (xi > 0);
        }
        else if (y1 > 0 && y2 <= 0)
        {
            real xi = x1 - y1 * (x2 - x1) / (y2 - y1);
            winding -= (xi > 0);
        }
    }
    return winding != 0;
}

real math_crossed_plane(real alpha, real beta, real gamma)
{
    alpha = fmod(alpha, CONST_2PI);
    alpha += (alpha < 0) * CONST_2PI;
    beta = fmod(beta, CONST_2PI);
    beta += (beta < 0) * CONST_2PI;

    int crosszero = fabs(alpha - beta) >= CONST_PI;
    alpha = crosszero ? fmod(alpha + CONST_PI, CONST_2PI) : alpha;
    beta = crosszero ? fmod(beta + CONST_PI, CONST_2PI) : beta;
    gamma = crosszero ? fmod(gamma + CONST_PI, CONST_2PI) : gamma;

    int betaissmaller = alpha <= beta;
    real smaller = betaissmaller ? alpha : beta;
    real larger = betaissmaller ? beta : alpha;

    real distance = gamma - smaller;
    real arclength = larger - smaller;

    int isinside = distance <= arclength;
    real nonzeroarc = arclength + (arclength == 0.0);
    real normalized_distance =
        betaissmaller ? distance / nonzeroarc : 1 - distance / nonzeroarc;

    return isinside * normalized_distance + (1.0 - isinside) * (-1.0);
}

void math_test_eval_vector_operations(
    const real a[3], const real b[3], const real c[3], real dot[1],
    real cross[3], real triple[1], real det[1], real norm[1], real normc[1],
    real unit[3])
{
    dot[0] = math_dot(a, b);
    math_cross(a, b, cross);
    triple[0] = math_scalar_triple_product(a, b, c);
    det[0] = math_determinant3x3(
        a[0], b[0], c[0], a[1], b[1], c[1], a[2], b[2], c[2]);
    norm[0] = math_norm(a);
    normc[0] = math_normc(a[0], a[1], a[2]);
    math_unit(a, unit);
}

void math_test_eval_vector_transformations(
    const real xyz[3], const real rpz[3], const real vxyz[3],
    const real vrpz[3], real xyz_out[3], real rpz_out[3], real vxyz_out[3],
    real vrpz_out[3])
{
    math_rpz2xyz(rpz, xyz_out);
    math_xyz2rpz(xyz, rpz_out);
    math_vec_rpz2xyz(vrpz, vxyz_out, rpz[1]);
    math_vec_xyz2rpz(vxyz, vrpz_out, rpz[1]);
}

void math_test_eval_bin_index(
    size_t nx, real xmin, real xmax, real x, size_t bin_index[1])
{
    bin_index[0] = math_bin_index(x, nx, xmin, xmax);
}

void math_test_eval_fmod(real a, real b, real out[1]) { out[0] = fmod(a, b); }

void math_test_eval_iabs(int a, int b[1]) { b[0] = math_iabs(a); }
