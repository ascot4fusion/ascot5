/**
 * Cubic 2D spline interpolation in compact form (see interp.h).
 */
#include "defines.h"
#include "interp.h"
#include "utils/mathlib.h"
#include <math.h>
#include <stdlib.h>

int Spline2D_init(
    Spline2D *spline, size_t nx, size_t ny, BoundaryCondition xbc,
    BoundaryCondition ybc, real xlim[2], real ylim[2], real f[nx * ny])
{
    int err = 0;
    real *fx = (real *)xmalloc(&err, nx * sizeof(real));
    if (err)
        return 1;

    real *fy = (real *)xmalloc(&err, ny * sizeof(real));
    if (err)
    {
        free(fx);
        return 1;
    }
    real *cx = (real *)xmalloc(&err, nx * NSIZE_COMP1D * sizeof(real));
    if (err)
    {
        free(fx);
        free(fy);
        return 1;
    }
    real *cy = (real *)xmalloc(&err, ny * NSIZE_COMP1D * sizeof(real));
    if (err)
    {
        free(fx);
        free(fy);
        free(cx);
        return 1;
    }
    spline->c =
        (real *)xmalloc(&err, ny * nx * NSIZE_COMP2D * sizeof(real));
    if (err)
    {
        free(fx);
        free(fy);
        free(cx);
        free(cy);
        return 1;
    }

    /* With the periodic boundary condition there's an "extra" grid interval
       that connects the last and first data points. */
    spline->dx = (xlim[1] - xlim[0]) / (nx - 1 + (xbc == PERIODICBC));
    spline->dy = (ylim[1] - ylim[0]) / (ny - 1 + (ybc == PERIODICBC));

    /* Cubic spline along x for each y, using f values to get fyy */
    for (size_t ix = 0; ix < nx; ix++)
    {
        for (size_t iy = 0; iy < ny; iy++)
            fy[iy] = f[ix * ny + iy];

        err += solve_compact_cubic_spline(ny, cy, fy, ybc);
        for (size_t iy = 0; iy < ny; iy++)
        {
            size_t idx = (ix * ny + iy) * NSIZE_COMP2D;
            spline->c[idx] = cy[iy * 2]; // f
            spline->c[idx + 2] =
                cy[iy * 2 + 1] / (spline->dy * spline->dy); // fyy
        }
    }

    /* Two cubic splines along y for each x, using f and fyy to get fxx and
       fxy */
    for (size_t iy = 0; iy < ny; iy++)
    {
        for (size_t ix = 0; ix < nx; ix++)
            fx[ix] = f[ix * ny + iy];

        err += solve_compact_cubic_spline(nx, cx, fx, xbc);
        for (size_t ix = 0; ix < nx; ix++)
        {
            size_t idx = (ix * ny + iy) * NSIZE_COMP2D;
            spline->c[idx + 1] =
                cx[ix * 2 + 1] / (spline->dx * spline->dx); // fxx

            fx[ix] = spline->c[idx + 2];
        }
        err += solve_compact_cubic_spline(nx, cx, fx, xbc);
        for (size_t ix = 0; ix < nx; ix++)
        {
            size_t idx = (ix * ny + iy) * NSIZE_COMP2D;
            spline->c[idx + 3] =
                cx[ix * 2 + 1] / (spline->dx * spline->dx); // fxxyy
        }
    }

    free(fx);
    free(fy);
    free(cx);
    free(cy);
    if (err)
    {
        free(spline->c);
        return 1;
    }

    spline->nx = nx;
    spline->ny = ny;
    spline->xbc = xbc;
    spline->ybc = ybc;
    spline->xlim[0] = xlim[0];
    spline->xlim[1] = xlim[1];
    spline->ylim[0] = ylim[0];
    spline->ylim[1] = ylim[1];

    return 0;
}

int Spline2D_eval_f(real f[1], Spline2D *spline, real x, real y)
{

    /* Make sure periodic coordinates are within [min, max] region. */
    x = spline->xbc == PERIODICBC
            ? fmod(x - spline->xlim[0], spline->xlim[1] - spline->xlim[0]) +
                  spline->xlim[0]
            : x;
    x = spline->xbc == PERIODICBC
            ? x + (x < spline->xlim[0]) * (spline->xlim[1] - spline->xlim[0])
            : x;

    y = spline->ybc == PERIODICBC
            ? fmod(y - spline->ylim[0], spline->ylim[1] - spline->ylim[0]) +
                  spline->ylim[0]
            : y;
    y = spline->ybc == PERIODICBC
            ? y + (y < spline->ylim[0]) * (spline->ylim[1] - spline->ylim[0])
            : y;

    /* Abscissa indices. The last term is needed to include the last point in
       the interpolation range. */
    size_t ix = (x - spline->xlim[0]) / spline->dx - 1 * (x == spline->xlim[1]);
    size_t iy = (y - spline->ylim[0]) / spline->dy - 1 * (y == spline->ylim[1]);

    /* Normalized coordinates inside the cell */
    real dx = (x - (spline->xlim[0] + ix * spline->dx)) / spline->dx;
    real dy = (y - (spline->ylim[0] + iy * spline->dy)) / spline->dy;

    /* Helper variables */
    real dx3 = dx * (dx * dx - 1.0), dy3 = dy * (dy * dy - 1.0);
    real dxi = 1.0 - dx, dyi = 1.0 - dy;
    real dxi3 = dxi * (dxi * dxi - 1.0), dyi3 = dyi * (dyi * dyi - 1.0);
    real xg2 = spline->dx * spline->dx, yg2 = spline->dy * spline->dy;

    /* Indices for locating correct c */
    size_t ystep = NSIZE_COMP2D;
    size_t xstep = spline->ny * NSIZE_COMP2D;
    size_t n = ix * spline->ny * NSIZE_COMP2D + iy * NSIZE_COMP2D;

    xstep = (spline->xbc == PERIODICBC && ix == spline->nx - 1)
                ? -(spline->nx - 1) * xstep
                : xstep;
    ystep = (spline->ybc == PERIODICBC && iy == spline->ny - 1)
                ? -(spline->ny - 1) * ystep
                : ystep;

    int err = (spline->xbc == NATURALBC &&
               !(x >= spline->xlim[0] && x <= spline->xlim[1])) ||
              (spline->ybc == NATURALBC &&
               !(y >= spline->ylim[0] && y <= spline->ylim[1]));

    n = err ? 0 : n;
    xstep = err ? 0 : xstep;
    ystep = err ? 0 : ystep;

    real c000 = spline->c[n + 0];
    real c001 = spline->c[n + 1];
    real c002 = spline->c[n + 2];
    real c003 = spline->c[n + 3];

    real c010 = spline->c[n + xstep + 0];
    real c011 = spline->c[n + xstep + 1];
    real c012 = spline->c[n + xstep + 2];
    real c013 = spline->c[n + xstep + 3];

    real c100 = spline->c[n + ystep + 0];
    real c101 = spline->c[n + ystep + 1];
    real c102 = spline->c[n + ystep + 2];
    real c103 = spline->c[n + ystep + 3];

    real c110 = spline->c[n + ystep + xstep + 0];
    real c111 = spline->c[n + ystep + xstep + 1];
    real c112 = spline->c[n + ystep + xstep + 2];
    real c113 = spline->c[n + ystep + xstep + 3];

    f[0] = (dxi * (dyi * c000 + dy * c100) + dx * (dyi * c010 + dy * c110)) +
           (xg2 / 6) * (dxi3 * (dyi * c001 + dy * c101) +
                        dx3 * (dyi * c011 + dy * c111)) +
           (yg2 / 6) * (dxi * (dyi3 * c002 + dy3 * c102) +
                        dx * (dyi3 * c012 + dy3 * c112)) +
           (xg2 * yg2 / 36) * (dxi3 * (dyi3 * c003 + dy3 * c103) +
                               dx3 * (dyi3 * c013 + dy3 * c113));
    return err;
}

int Spline2D_eval_f_df(real f_df[6], Spline2D *spline, real x, real y)
{

    /* Make sure periodic coordinates are within [min, max] region. */
    x = spline->xbc == PERIODICBC
            ? fmod(x - spline->xlim[0], spline->xlim[1] - spline->xlim[0]) +
                  spline->xlim[0]
            : x;
    x = spline->xbc == PERIODICBC
            ? x + (x < spline->xlim[0]) * (spline->xlim[1] - spline->xlim[0])
            : x;

    y = spline->ybc == PERIODICBC
            ? fmod(y - spline->ylim[0], spline->ylim[1] - spline->ylim[0]) +
                  spline->ylim[0]
            : y;
    y = spline->ybc == PERIODICBC
            ? y + (y < spline->ylim[0]) * (spline->ylim[1] - spline->ylim[0])
            : y;

    int err = (spline->xbc == NATURALBC &&
               !(x >= spline->xlim[0] && x <= spline->xlim[1])) ||
              (spline->ybc == NATURALBC &&
               !(y >= spline->ylim[0] && y <= spline->ylim[1]));

    /* Abscissa indices. The last term is needed to include the last point in
       the interpolation range. */
    size_t ix = (x - spline->xlim[0]) / spline->dx - 1 * (x == spline->xlim[1]);
    size_t iy = (y - spline->ylim[0]) / spline->dy - 1 * (y == spline->ylim[1]);

    /* Normalize coordinates inside the cell */
    x = (x - spline->xlim[0] - ix * spline->dx) / spline->dx;
    y = (y - spline->ylim[0] - iy * spline->dy) / spline->dy;

    /* Helper variables */
    real dx3 = x * (x * x - 1.0), dy3 = y * (y * y - 1.0);
    real dxi = 1.0 - x, dyi = 1.0 - y;
    real dxi3 = dxi * (dxi * dxi - 1.0), dyi3 = dyi * (dyi * dyi - 1.0);
    real dx = spline->dx, dy = spline->dy;
    real dx2 = dx * dx, dy2 = dy * dy;
    real dx3dx = 3 * x * x - 1, dy3dy = 3 * y * y - 1;
    real dxi3dx = -3 * dxi * dxi + 1, dyi3dy = -3 * dyi * dyi + 1;
    real xgi = 1.0 / spline->dx, ygi = 1.0 / spline->dy;

    /* Indices for locating correct c */
    size_t ystep = NSIZE_COMP2D;
    size_t xstep = spline->ny * NSIZE_COMP2D;
    size_t n = ix * spline->ny * NSIZE_COMP2D + iy * NSIZE_COMP2D;

    xstep = (spline->xbc == PERIODICBC && ix == spline->nx - 1)
                ? -(spline->nx - 1) * xstep
                : xstep;
    ystep = (spline->ybc == PERIODICBC && iy == spline->ny - 1)
                ? -(spline->ny - 1) * ystep
                : ystep;

    n = err ? 0 : n;
    xstep = err ? 0 : xstep;
    ystep = err ? 0 : ystep;

    real c000 = spline->c[n + 0];
    real c001 = spline->c[n + 1];
    real c002 = spline->c[n + 2];
    real c003 = spline->c[n + 3];

    real c010 = spline->c[n + xstep + 0];
    real c011 = spline->c[n + xstep + 1];
    real c012 = spline->c[n + xstep + 2];
    real c013 = spline->c[n + xstep + 3];

    real c100 = spline->c[n + ystep + 0];
    real c101 = spline->c[n + ystep + 1];
    real c102 = spline->c[n + ystep + 2];
    real c103 = spline->c[n + ystep + 3];

    real c110 = spline->c[n + ystep + xstep + 0];
    real c111 = spline->c[n + ystep + xstep + 1];
    real c112 = spline->c[n + ystep + xstep + 2];
    real c113 = spline->c[n + ystep + xstep + 3];

    f_df[0] = (dxi * (dyi * c000 + y * c100) + x * (dyi * c010 + y * c110)) +
              (dx2 / 6) * (dxi3 * (dyi * c001 + y * c101) +
                           dx3 * (dyi * c011 + y * c111)) +
              (dy2 / 6) * (dxi * (dyi3 * c002 + dy3 * c102) +
                           x * (dyi3 * c012 + dy3 * c112)) +
              (dx2 * dy2 / 36) * (dxi3 * (dyi3 * c003 + dy3 * c103) +
                                  dx3 * (dyi3 * c013 + dy3 * c113));

    f_df[1] = xgi * (-(dyi * c000 + y * c100) + (dyi * c010 + y * c110)) +
              (dx / 6) * (dxi3dx * (dyi * c001 + y * c101) +
                          dx3dx * (dyi * c011 + y * c111)) +
              (xgi * dy2 / 6) *
                  (-(dyi3 * c002 + dy3 * c102) + (dyi3 * c012 + dy3 * c112)) +
              (dx * dy2 / 36) * (dxi3dx * (dyi3 * c003 + dy3 * c103) +
                                 dx3dx * (dyi3 * c013 + dy3 * c113));

    f_df[2] = ygi * (dxi * (-c000 + c100) + x * (-c010 + c110)) +
              (dx2 * ygi / 6) * (dxi3 * (-c001 + c101) + dx3 * (-c011 + c111)) +
              (dy / 6) * (dxi * (dyi3dy * c002 + dy3dy * c102) +
                          x * (dyi3dy * c012 + dy3dy * c112)) +
              (dx2 * dy / 36) * (dxi3 * (dyi3dy * c003 + dy3dy * c103) +
                                 dx3 * (dyi3dy * c013 + dy3dy * c113));

    f_df[3] = (dxi * (dyi * c001 + y * c101) + x * (dyi * c011 + y * c111)) +
              (dy2 / 6) * (dxi * (dyi3 * c003 + dy3 * c103) +
                           x * (dyi3 * c013 + dy3 * c113));

    f_df[4] =
        (dxi * (dyi * c002 + y * c102) + x * (dyi * c012 + y * c112)) +
        dx2 / 6 *
            (dxi3 * (dyi * c003 + y * c103) + dx3 * (dyi * c013 + y * c113));

    f_df[5] =
        xgi * ygi * (c000 - c100 - c010 + c110) +
        (dx / 6 * ygi) * (dxi3dx * (-c001 + c101) + dx3dx * (-c011 + c111)) +
        (xgi / 6 * dy) *
            (-(dyi3dy * c002 + dy3dy * c102) + (dyi3dy * c012 + dy3dy * c112)) +
        (dx * dy / 36) * (dxi3dx * (dyi3dy * c003 + dy3dy * c103) +
                          dx3dx * (dyi3dy * c013 + dy3dy * c113));
    return err;
}
