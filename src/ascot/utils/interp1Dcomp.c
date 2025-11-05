/**
 * Cubic 1D spline interpolation in compact form (see interp.h).
 */
#include "defines.h"
#include "interp.h"
#include "utils/mathlib.h"
#include <math.h>
#include <stdlib.h>

int Spline1D_init(
    Spline1D *spline, size_t nx, BoundaryCondition xbc, real xlim[2],
    real f[nx])
{
    int err = 0;
    spline->c =
        (real *)xmalloc(&err, nx * NSIZE_COMP2D * sizeof(real));
    if (err)
        return 1;

    /* With the periodic boundary condition there's an "extra" grid interval
       that connects the last and first data points. */
    spline->dx = (xlim[1] - xlim[0]) / (nx - 1 + (xbc == PERIODICBC));

    /* Cubic spline along x, using f values to get fxx */
    err = solve_compact_cubic_spline(nx, spline->c, f, xbc);
    if (err)
    {
        free(spline->c);
        return 1;
    }
    for (size_t ix = 0; ix < nx; ix++)
        spline->c[ix * 2 + 1] /= (spline->dx * spline->dx);

    spline->nx = nx;
    spline->xbc = xbc;
    spline->xlim[0] = xlim[0];
    spline->xlim[1] = xlim[1];
    return 0;
}

int Spline1D_eval_f(real f[1], Spline1D *spline, real x)
{
    /* Make sure periodic coordinates are within [min, max] region. */
    x = spline->xbc == PERIODICBC
            ? fmod(x - spline->xlim[0], spline->xlim[1] - spline->xlim[0]) +
                  spline->xlim[0]
            : x;
    x = spline->xbc == PERIODICBC
            ? x + (x < spline->xlim[0]) * (spline->xlim[1] - spline->xlim[0])
            : x;

    /* Abscissa indices. The last term is needed to include the last point in
       the interpolation range. */
    size_t ix = (x - spline->xlim[0]) / spline->dx - 1 * (x == spline->xlim[1]);

    /* Normalized coordinates inside the cell */
    real dx = (x - (spline->xlim[0] + ix * spline->dx)) / spline->dx;

    /* Helper variables */
    real dx3 = dx * (dx * dx - 1.0);
    real dxi = 1.0 - dx;
    real dxi3 = dxi * (dxi * dxi - 1.0);
    real xg2 = spline->dx * spline->dx;

    /* Indices for locating correct coefficients */
    size_t xstep = NSIZE_COMP1D;
    size_t n = ix * NSIZE_COMP1D;

    xstep = (spline->xbc == PERIODICBC && ix == spline->nx - 1)
                ? -(spline->nx - 1) * xstep
                : xstep;

    int err = spline->xbc == NATURALBC &&
              !(x >= spline->xlim[0] && x <= spline->xlim[1]);

    n = err ? 0 : n;
    xstep = err ? 0 : xstep;

    real c00 = spline->c[n + 0];
    real c01 = spline->c[n + 1];
    real c10 = spline->c[n + xstep + 0];
    real c11 = spline->c[n + xstep + 1];

    f[0] = dxi * c00 + dx * c10 + (xg2 / 6) * (dxi3 * c01 + dx3 * c11);
    return err;
}

int Spline1D_eval_f_df(real f_df[3], Spline1D *spline, real x)
{

    /* Make sure periodic coordinates are within [min, max] region. */
    x = spline->xbc == PERIODICBC
            ? fmod(x - spline->xlim[0], spline->xlim[1] - spline->xlim[0]) +
                  spline->xlim[0]
            : x;
    x = spline->xbc == PERIODICBC
            ? x + (x < spline->xlim[0]) * (spline->xlim[1] - spline->xlim[0])
            : x;

    /* Abscissa indices. The last term is needed to include the last point in
       the interpolation range. */
    size_t ix = (x - spline->xlim[0]) / spline->dx - 1 * (x == spline->xlim[1]);

    /* Normalized coordinates inside the cell */
    real dx = (x - (spline->xlim[0] + ix * spline->dx)) / spline->dx;

    /* Helper variables */
    real dx3 = dx * (dx * dx - 1.0);
    real dxi = 1.0 - dx;
    real dxi3 = dxi * (dxi * dxi - 1.0);
    real xg = spline->dx;
    real xg2 = spline->dx * spline->dx;
    real dx3dx = 3 * dx * dx - 1;
    real dxi3dx = -3 * dxi * dxi + 1;
    real xgi = 1.0 / spline->dx;

    /* Indices for locating correct coefficients */
    size_t xstep = NSIZE_COMP1D;
    size_t n = ix * NSIZE_COMP1D;

    xstep = (spline->xbc == PERIODICBC && ix == spline->nx - 1)
                ? -(spline->nx - 1) * xstep
                : xstep;

    int err = spline->xbc == NATURALBC &&
              !(x >= spline->xlim[0] && x <= spline->xlim[1]);

    n = err ? 0 : n;
    xstep = err ? 0 : xstep;

    real c00 = spline->c[n + 0];
    real c01 = spline->c[n + 1];
    real c10 = spline->c[n + xstep + 0];
    real c11 = spline->c[n + xstep + 1];

    f_df[0] = dxi * c00 + dx * c10 + (xg2 / 6) * (dxi3 * c01 + dx3 * c11);
    f_df[1] = xgi * (c10 - c00) + (xg / 6) * (dx3dx * c11 + dxi3dx * c01);
    f_df[2] = dxi * c01 + dx * c11;
    return err;
}
