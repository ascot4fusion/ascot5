/**
 * Implements biosaw.h.
 */
#include "biosaw.h"
#include "ascot5.h"
#include "consts.h"
#include "math.h"

void biosaw_calc_b(
    int n, int coil_n, real *x, real *y, real *z, real *coil_x, real *coil_y,
    real *coil_z, real *bx, real *by, real *bz)
{

#pragma omp parallel for
    for (int ix = 0; ix < n; ix++)
    {
        bx[ix] = 0;
        by[ix] = 0;
        bz[ix] = 0;

        real p1[3], p2[3];
        p2[0] = coil_x[0];
        p2[1] = coil_y[0];
        p2[2] = coil_z[0];

        for (int i = 1; i < coil_n; i++)
        {
            math_copy(p1, p2);

            p2[0] = coil_x[i];
            p2[1] = coil_y[i];
            p2[2] = coil_z[i];

            real p1p2[3];
            p1p2[0] = p2[0] - p1[0];
            p1p2[1] = p2[1] - p1[1];
            p1p2[2] = p2[2] - p1[2];

            real p1x[3];
            p1x[0] = x[ix] - p1[0];
            p1x[1] = y[ix] - p1[1];
            p1x[2] = z[ix] - p1[2];

            real p2x[3];
            p2x[0] = x[ix] - p2[0];
            p2x[1] = y[ix] - p2[1];
            p2x[2] = z[ix] - p2[2];

            real d1 = math_norm(p1x);
            real d2 = math_norm(p2x);
            real l = math_norm(p1p2);
            real s = math_dot(p1p2, p1x) / math_dot(p1p2, p1p2);
            real h = s * l;

            real xs[3];
            xs[0] = p1[0] + s * p1p2[0] - x[ix];
            xs[1] = p1[1] + s * p1p2[1] - y[ix];
            xs[2] = p1[2] + s * p1p2[2] - z[ix];

            real d = math_norm(xs);

            real abs_B =
                CONST_MU0 / (4 * CONST_PI) * ((l - h) / d2 + h / d1) / d;

            real dir_B[3];
            math_cross(p1x, p2x, dir_B);

            bx[ix] += abs_B * dir_B[0] / math_norm(dir_B);
            by[ix] += abs_B * dir_B[1] / math_norm(dir_B);
            bz[ix] += abs_B * dir_B[2] / math_norm(dir_B);
        }
    }
}
