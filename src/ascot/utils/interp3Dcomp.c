/**
 * Cubic 3D spline interpolation in compact form (see interp.h).
 */
#include "consts.h"
#include "defines.h"
#include "interp.h"
#include "utils/mathlib.h"
#include <math.h>
#include <stdlib.h>

int Spline3D_init(
    Spline3D *spline, size_t nx, size_t ny, size_t nz, BoundaryCondition xbc,
    BoundaryCondition ybc, BoundaryCondition zbc, real xlim[2], real ylim[2],
    real zlim[2], real f[nx * ny * nz])
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
    real *fz = (real *)xmalloc(&err, nz * sizeof(real));
    if (err)
    {
        free(fx);
        free(fy);
        return 1;
    }
    real *cx = (real *)xmalloc(&err, nx * NSIZE_COMP1D * sizeof(real));
    if (err)
    {
        free(fx);
        free(fy);
        free(fz);
        return 1;
    }
    real *cy = (real *)xmalloc(&err, ny * NSIZE_COMP1D * sizeof(real));
    if (err)
    {
        free(fx);
        free(fy);
        free(fz);
        free(cx);
        return 1;
    }
    real *cz = (real *)xmalloc(&err, nz * NSIZE_COMP1D * sizeof(real));
    if (err)
    {
        free(fx);
        free(fy);
        free(fz);
        free(cx);
        free(cy);
        return 1;
    }
    spline->c =
        (real *)xmalloc(&err, nz * ny * nx * NSIZE_COMP3D * sizeof(real));
    if (err)
    {
        free(fx);
        free(fy);
        free(fz);
        free(cx);
        free(cy);
        free(cz);
        return 1;
    }

    /* With the periodic boundary condition there's an "extra" grid interval
       that connects the last and first data points. */
    spline->dx = (xlim[1] - xlim[0]) / (nx - 1 + (xbc == PERIODICBC));
    spline->dy = (ylim[1] - ylim[0]) / (ny - 1 + (ybc == PERIODICBC));
    spline->dz = (zlim[1] - zlim[0]) / (nz - 1 + (zbc == PERIODICBC));

    /* Bicubic spline surfaces over yz-grid for each x */
    for (size_t ix = 0; ix < nx; ix++)
    {
        /* Cubic spline along z for each y, using f values to get fzz */
        for (size_t iy = 0; iy < ny; iy++)
        {
            for (size_t iz = 0; iz < nz; iz++)
            {
                fz[iz] = f[ix * ny * nz + iy * nz + iz];
            }
            err += solve_compact_cubic_spline(nz, cz, fz, zbc);
            for (size_t iz = 0; iz < nz; iz++)
            {
                size_t idx = (ix * ny * nz + iy * nz + iz) * NSIZE_COMP3D;
                spline->c[idx] = cz[iz * 2]; // f
                spline->c[idx + 3] =
                    cz[iz * 2 + 1] / (spline->dz * spline->dz); // fzz
            }
        }

        /* Two cubic splines along z for each y, one using f values to
           get fyy, and the other using fzz values to get fyz */
        for (size_t iz = 0; iz < nz; iz++)
        {
            for (size_t iy = 0; iy < ny; iy++)
            {
                fy[iy] = f[ix * ny * nz + iy * nz + iz];
            }
            err += solve_compact_cubic_spline(ny, cy, fy, ybc);
            for (size_t iy = 0; iy < ny; iy++)
            {
                size_t idx = (ix * ny * nz + iy * nz + iz) * NSIZE_COMP3D;
                spline->c[idx + 2] =
                    cy[iy * 2 + 1] / (spline->dy * spline->dy); // fyy
            }
            for (size_t iy = 0; iy < ny; iy++)
            {
                size_t idx = (ix * ny * nz + iy * nz + iz) * NSIZE_COMP3D;
                fy[iy] = spline->c[idx + 3];
            }
            err += solve_compact_cubic_spline(ny, cy, fy, ybc);
            for (size_t iy = 0; iy < ny; iy++)
            {
                size_t idx = (ix * ny * nz + iy * nz + iz) * NSIZE_COMP3D;
                spline->c[idx + 6] =
                    cy[iy * 2 + 1] / (spline->dy * spline->dy); // fyyzz
            }
        }
    }

    /* Four cubic splines along x for each yz-pair, one using f values to get
       fxx, one using fzz to get fxz, one using fyy to get fxy, and one
       using fyz to get fxyz */
    for (size_t iy = 0; iy < ny; iy++)
    {
        for (size_t iz = 0; iz < nz; iz++)
        {
            for (size_t ix = 0; ix < nx; ix++)
            {
                fx[ix] = f[ix * ny * nz + iy * nz + iz];
            }
            err += solve_compact_cubic_spline(nx, cx, fx, xbc);
            for (size_t ix = 0; ix < nx; ix++)
            {
                size_t idx = (ix * ny * nz + iy * nz + iz) * NSIZE_COMP3D;
                spline->c[idx + 1] =
                    cx[ix * 2 + 1] / (spline->dx * spline->dx); // fxx
            }
            for (size_t ix = 0; ix < nx; ix++)
            {
                size_t idx = (ix * ny * nz + iy * nz + iz) * NSIZE_COMP3D;
                fx[ix] = spline->c[idx + 3];
            }
            err += solve_compact_cubic_spline(nx, cx, fx, xbc);
            for (size_t ix = 0; ix < nx; ix++)
            {
                size_t idx = (ix * ny * nz + iy * nz + iz) * NSIZE_COMP3D;
                spline->c[idx + 5] =
                    cx[ix * 2 + 1] / (spline->dx * spline->dx); // fxxzz
            }
            for (size_t ix = 0; ix < nx; ix++)
            {
                size_t idx = (ix * ny * nz + iy * nz + iz) * NSIZE_COMP3D;
                fx[ix] = spline->c[idx + 2];
            }
            err += solve_compact_cubic_spline(nx, cx, fx, xbc);
            for (size_t ix = 0; ix < nx; ix++)
            {
                size_t idx = (ix * ny * nz + iy * nz + iz) * NSIZE_COMP3D;
                spline->c[idx + 4] =
                    cx[ix * 2 + 1] / (spline->dx * spline->dx); // fxxyy
            }
            for (size_t ix = 0; ix < nx; ix++)
            {
                size_t idx = (ix * ny * nz + iy * nz + iz) * NSIZE_COMP3D;
                fx[ix] = spline->c[idx + 6];
            }
            err += solve_compact_cubic_spline(nx, cx, fx, xbc);
            for (size_t ix = 0; ix < nx; ix++)
            {
                size_t idx = (ix * ny * nz + iy * nz + iz) * NSIZE_COMP3D;
                spline->c[idx + 7] =
                    cx[ix * 2 + 1] / (spline->dx * spline->dx); // fxxyyzz
            }
        }
    }

    free(fx);
    free(fy);
    free(fz);
    free(cx);
    free(cy);
    free(cz);
    if (err) {
        free(spline->c);
        return 1;
    }

    spline->nx = nx;
    spline->ny = ny;
    spline->nz = nz;
    spline->xbc = xbc;
    spline->ybc = ybc;
    spline->zbc = zbc;
    spline->xlim[0] = xlim[0];
    spline->xlim[1] = xlim[1];
    spline->ylim[0] = ylim[0];
    spline->ylim[1] = ylim[1];
    spline->zlim[0] = zlim[0];
    spline->zlim[1] = zlim[1];
    return 0;
}

int Spline3D_eval_f(real f[1], Spline3D *spline, real x, real y, real z)
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
    z = spline->zbc == PERIODICBC
            ? fmod(z - spline->zlim[0], spline->zlim[1] - spline->zlim[0]) +
                  spline->zlim[0]
            : z;
    z = spline->zbc == PERIODICBC
            ? z + (z < spline->zlim[0]) * (spline->zlim[1] - spline->zlim[0])
            : z;

    /* Abscissa indices. The last term is needed to include the last point in
       the interpolation range. */
    size_t ix = (x - spline->xlim[0]) / spline->dx - 1 * (x == spline->xlim[1]);
    size_t iy = (y - spline->ylim[0]) / spline->dy - 1 * (y == spline->ylim[1]);
    size_t iz = (z - spline->zlim[0]) / spline->dz - 1 * (z == spline->zlim[1]);

    /* Normalized coordinates inside the cell */
    real dx = (x - (spline->xlim[0] + ix * spline->dx)) / spline->dx;
    real dy = (y - (spline->ylim[0] + iy * spline->dy)) / spline->dy;
    real dz = (z - (spline->zlim[0] + iz * spline->dz)) / spline->dz;

    /* Helper variables */
    real dx3 = dx * (dx * dx - 1.0), dy3 = dy * (dy * dy - 1.0),
         dz3 = dz * (dz * dz - 1.0);
    real dxi = 1.0 - dx, dyi = 1.0 - dy, dzi = 1.0 - dz;
    real dxi3 = dxi * (dxi * dxi - 1.0), dyi3 = dyi * (dyi * dyi - 1.0),
         dzi3 = dzi * (dzi * dzi - 1.0);
    real xg2 = spline->dx * spline->dx, yg2 = spline->dy * spline->dy,
         zg2 = spline->dz * spline->dz;

    /* Indices for locating correct c */
    size_t zstep = NSIZE_COMP3D;
    size_t ystep = spline->nz * NSIZE_COMP3D;
    size_t xstep = spline->nz * spline->ny * NSIZE_COMP3D;
    size_t n =
        (ix * spline->nz * spline->ny + iy * spline->nz + iz) * NSIZE_COMP3D;

    xstep = (spline->xbc == PERIODICBC && ix == spline->nx - 1)
                ? -(spline->nx - 1) * xstep
                : xstep;
    ystep = (spline->ybc == PERIODICBC && iy == spline->ny - 1)
                ? -(spline->ny - 1) * ystep
                : ystep;
    zstep = (spline->zbc == PERIODICBC && iz == spline->nz - 1)
                ? -(spline->nz - 1) * zstep
                : zstep;

    int err = (spline->xbc == NATURALBC &&
               !(x >= spline->xlim[0] && x <= spline->xlim[1])) ||
              (spline->ybc == NATURALBC &&
               !(y >= spline->ylim[0] && y <= spline->ylim[1])) ||
              (spline->zbc == NATURALBC &&
               !(z >= spline->zlim[0] && z <= spline->zlim[1]));

    n = err ? 0 : n;
    xstep = err ? 0 : xstep;
    ystep = err ? 0 : ystep;
    zstep = err ? 0 : zstep;

    real c0000 = spline->c[n + 0];
    real c0001 = spline->c[n + 1];
    real c0002 = spline->c[n + 2];
    real c0003 = spline->c[n + 3];
    real c0004 = spline->c[n + 4];
    real c0005 = spline->c[n + 5];
    real c0006 = spline->c[n + 6];
    real c0007 = spline->c[n + 7];

    real c0010 = spline->c[n + xstep + 0];
    real c0011 = spline->c[n + xstep + 1];
    real c0012 = spline->c[n + xstep + 2];
    real c0013 = spline->c[n + xstep + 3];
    real c0014 = spline->c[n + xstep + 4];
    real c0015 = spline->c[n + xstep + 5];
    real c0016 = spline->c[n + xstep + 6];
    real c0017 = spline->c[n + xstep + 7];

    real c0100 = spline->c[n + ystep + 0];
    real c0101 = spline->c[n + ystep + 1];
    real c0102 = spline->c[n + ystep + 2];
    real c0103 = spline->c[n + ystep + 3];
    real c0104 = spline->c[n + ystep + 4];
    real c0105 = spline->c[n + ystep + 5];
    real c0106 = spline->c[n + ystep + 6];
    real c0107 = spline->c[n + ystep + 7];

    real c1000 = spline->c[n + zstep + 0];
    real c1001 = spline->c[n + zstep + 1];
    real c1002 = spline->c[n + zstep + 2];
    real c1003 = spline->c[n + zstep + 3];
    real c1004 = spline->c[n + zstep + 4];
    real c1005 = spline->c[n + zstep + 5];
    real c1006 = spline->c[n + zstep + 6];
    real c1007 = spline->c[n + zstep + 7];

    real c0110 = spline->c[n + ystep + xstep + 0];
    real c0111 = spline->c[n + ystep + xstep + 1];
    real c0112 = spline->c[n + ystep + xstep + 2];
    real c0113 = spline->c[n + ystep + xstep + 3];
    real c0114 = spline->c[n + ystep + xstep + 4];
    real c0115 = spline->c[n + ystep + xstep + 5];
    real c0116 = spline->c[n + ystep + xstep + 6];
    real c0117 = spline->c[n + ystep + xstep + 7];

    real c1010 = spline->c[n + zstep + xstep + 0];
    real c1011 = spline->c[n + zstep + xstep + 1];
    real c1012 = spline->c[n + zstep + xstep + 2];
    real c1013 = spline->c[n + zstep + xstep + 3];
    real c1014 = spline->c[n + zstep + xstep + 4];
    real c1015 = spline->c[n + zstep + xstep + 5];
    real c1016 = spline->c[n + zstep + xstep + 6];
    real c1017 = spline->c[n + zstep + xstep + 7];

    real c1100 = spline->c[n + zstep + ystep + 0];
    real c1101 = spline->c[n + zstep + ystep + 1];
    real c1102 = spline->c[n + zstep + ystep + 2];
    real c1103 = spline->c[n + zstep + ystep + 3];
    real c1104 = spline->c[n + zstep + ystep + 4];
    real c1105 = spline->c[n + zstep + ystep + 5];
    real c1106 = spline->c[n + zstep + ystep + 6];
    real c1107 = spline->c[n + zstep + ystep + 7];

    real c1110 = spline->c[n + zstep + ystep + xstep + 0];
    real c1111 = spline->c[n + zstep + ystep + xstep + 1];
    real c1112 = spline->c[n + zstep + ystep + xstep + 2];
    real c1113 = spline->c[n + zstep + ystep + xstep + 3];
    real c1114 = spline->c[n + zstep + ystep + xstep + 4];
    real c1115 = spline->c[n + zstep + ystep + xstep + 5];
    real c1116 = spline->c[n + zstep + ystep + xstep + 6];
    real c1117 = spline->c[n + zstep + ystep + xstep + 7];

    f[0] = (dzi * (dxi * (dyi * c0000 + dy * c0100) +
                   dx * (dyi * c0010 + dy * c0110)) +
            dz * (dxi * (dyi * c1000 + dy * c1100) +
                  dx * (dyi * c1010 + dy * c1110))) +
           xg2 / 6 *
               (dzi * (dxi3 * (dyi * c0001 + dy * c0101) +
                       dx3 * (dyi * c0011 + dy * c0111)) +
                dz * (dxi3 * (dyi * c1001 + dy * c1101) +
                      dx3 * (dyi * c1011 + dy * c1111))) +
           yg2 / 6 *
               (dzi * (dxi * (dyi3 * c0002 + dy3 * c0102) +
                       dx * (dyi3 * c0012 + dy3 * c0112)) +
                dz * (dxi * (dyi3 * c1002 + dy3 * c1102) +
                      dx * (dyi3 * c1012 + dy3 * c1112))) +
           zg2 / 6 *
               (dzi3 * (dxi * (dyi * c0003 + dy * c0103) +
                        dx * (dyi * c0013 + dy * c0113)) +
                dz3 * (dxi * (dyi * c1003 + dy * c1103) +
                       dx * (dyi * c1013 + dy * c1113))) +
           xg2 * yg2 / 36 *
               (dzi * (dxi3 * (dyi3 * c0004 + dy3 * c0104) +
                       dx3 * (dyi3 * c0014 + dy3 * c0114)) +
                dz * (dxi3 * (dyi3 * c1004 + dy3 * c1104) +
                      dx3 * (dyi3 * c1014 + dy3 * c1114))) +
           xg2 * zg2 / 36 *
               (dzi3 * (dxi3 * (dyi * c0005 + dy * c0105) +
                        dx3 * (dyi * c0015 + dy * c0115)) +
                dz3 * (dxi3 * (dyi * c1005 + dy * c1105) +
                       dx3 * (dyi * c1015 + dy * c1115))) +
           yg2 * zg2 / 36 *
               (dzi3 * (dxi * (dyi3 * c0006 + dy3 * c0106) +
                        dx * (dyi3 * c0016 + dy3 * c0116)) +
                dz3 * (dxi * (dyi3 * c1006 + dy3 * c1106) +
                       dx * (dyi3 * c1016 + dy3 * c1116))) +
           xg2 * yg2 * zg2 / 216 *
               (dzi3 * (dxi3 * (dyi3 * c0007 + dy3 * c0107) +
                        dx3 * (dyi3 * c0017 + dy3 * c0117)) +
                dz3 * (dxi3 * (dyi3 * c1007 + dy3 * c1107) +
                       dx3 * (dyi3 * c1017 + dy3 * c1117)));

    return err;
}

int Spline3D_eval_f_df(real f_df[10], Spline3D *spline, real x, real y, real z)
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
    z = spline->zbc == PERIODICBC
            ? fmod(z - spline->zlim[0], spline->zlim[1] - spline->zlim[0]) +
                  spline->zlim[0]
            : z;
    z = spline->zbc == PERIODICBC
            ? z + (z < spline->zlim[0]) * (spline->zlim[1] - spline->zlim[0])
            : z;

    /* Abscissa indices. The last term is needed to include the last point in
       the interpolation range. */
    size_t ix = (x - spline->xlim[0]) / spline->dx - 1 * (x == spline->xlim[1]);
    size_t iy = (y - spline->ylim[0]) / spline->dy - 1 * (y == spline->ylim[1]);
    size_t iz = (z - spline->zlim[0]) / spline->dz - 1 * (z == spline->zlim[1]);

    /* Normalized coordinates inside the cell */
    real dx = (x - (spline->xlim[0] + ix * spline->dx)) / spline->dx;
    real dy = (y - (spline->ylim[0] + iy * spline->dy)) / spline->dy;
    real dz = (z - (spline->zlim[0] + iz * spline->dz)) / spline->dz;

    /* Helper variables */
    real dx3 = dx * (dx * dx - 1.0), dy3 = dy * (dy * dy - 1.0),
         dz3 = dz * (dz * dz - 1.0);
    real dxi = 1.0 - dx, dyi = 1.0 - dy, dzi = 1.0 - dz;
    real dxi3 = dxi * (dxi * dxi - 1.0), dyi3 = dyi * (dyi * dyi - 1.0),
         dzi3 = dzi * (dzi * dzi - 1.0);
    real xg = spline->dx, yg = spline->dy, zg = spline->dz;
    real xg2 = xg * xg, yg2 = yg * yg, zg2 = zg * zg;
    real dx3dx = 3 * dx * dx - 1, dy3dy = 3 * dy * dy - 1,
         dz3dz = 3 * dz * dz - 1;
    real dxi3dx = -3 * dxi * dxi + 1, dyi3dy = -3 * dyi * dyi + 1,
         dzi3dz = -3 * dzi * dzi + 1;
    real xgi = 1.0 / spline->dx, ygi = 1.0 / spline->dy, zgi = 1.0 / spline->dz;

    /* Indices for locating correct c */
    size_t zstep = NSIZE_COMP3D;
    size_t ystep = spline->nz * NSIZE_COMP3D;
    size_t xstep = spline->nz * spline->ny * NSIZE_COMP3D;
    size_t n =
        (ix * spline->nz * spline->ny + iy * spline->nz + iz) * NSIZE_COMP3D;

    xstep = (spline->xbc == PERIODICBC && ix == spline->nx - 1)
                ? -(spline->nx - 1) * xstep
                : xstep;
    ystep = (spline->ybc == PERIODICBC && iy == spline->ny - 1)
                ? -(spline->ny - 1) * ystep
                : ystep;
    zstep = (spline->zbc == PERIODICBC && iz == spline->nz - 1)
                ? -(spline->nz - 1) * zstep
                : zstep;

    int err = (spline->xbc == NATURALBC &&
               !(x >= spline->xlim[0] && x <= spline->xlim[1])) ||
              (spline->ybc == NATURALBC &&
               !(y >= spline->ylim[0] && y <= spline->ylim[1])) ||
              (spline->zbc == NATURALBC &&
               !(z >= spline->zlim[0] && z <= spline->zlim[1]));

    n = err ? 0 : n;
    xstep = err ? 0 : xstep;
    ystep = err ? 0 : ystep;
    zstep = err ? 0 : zstep;

    real c0000 = spline->c[n + 0];
    real c0001 = spline->c[n + 1];
    real c0002 = spline->c[n + 2];
    real c0003 = spline->c[n + 3];
    real c0004 = spline->c[n + 4];
    real c0005 = spline->c[n + 5];
    real c0006 = spline->c[n + 6];
    real c0007 = spline->c[n + 7];

    real c0010 = spline->c[n + xstep + 0];
    real c0011 = spline->c[n + xstep + 1];
    real c0012 = spline->c[n + xstep + 2];
    real c0013 = spline->c[n + xstep + 3];
    real c0014 = spline->c[n + xstep + 4];
    real c0015 = spline->c[n + xstep + 5];
    real c0016 = spline->c[n + xstep + 6];
    real c0017 = spline->c[n + xstep + 7];

    real c0100 = spline->c[n + ystep + 0];
    real c0101 = spline->c[n + ystep + 1];
    real c0102 = spline->c[n + ystep + 2];
    real c0103 = spline->c[n + ystep + 3];
    real c0104 = spline->c[n + ystep + 4];
    real c0105 = spline->c[n + ystep + 5];
    real c0106 = spline->c[n + ystep + 6];
    real c0107 = spline->c[n + ystep + 7];

    real c1000 = spline->c[n + zstep + 0];
    real c1001 = spline->c[n + zstep + 1];
    real c1002 = spline->c[n + zstep + 2];
    real c1003 = spline->c[n + zstep + 3];
    real c1004 = spline->c[n + zstep + 4];
    real c1005 = spline->c[n + zstep + 5];
    real c1006 = spline->c[n + zstep + 6];
    real c1007 = spline->c[n + zstep + 7];

    real c0110 = spline->c[n + ystep + xstep + 0];
    real c0111 = spline->c[n + ystep + xstep + 1];
    real c0112 = spline->c[n + ystep + xstep + 2];
    real c0113 = spline->c[n + ystep + xstep + 3];
    real c0114 = spline->c[n + ystep + xstep + 4];
    real c0115 = spline->c[n + ystep + xstep + 5];
    real c0116 = spline->c[n + ystep + xstep + 6];
    real c0117 = spline->c[n + ystep + xstep + 7];

    real c1010 = spline->c[n + zstep + xstep + 0];
    real c1011 = spline->c[n + zstep + xstep + 1];
    real c1012 = spline->c[n + zstep + xstep + 2];
    real c1013 = spline->c[n + zstep + xstep + 3];
    real c1014 = spline->c[n + zstep + xstep + 4];
    real c1015 = spline->c[n + zstep + xstep + 5];
    real c1016 = spline->c[n + zstep + xstep + 6];
    real c1017 = spline->c[n + zstep + xstep + 7];

    real c1100 = spline->c[n + zstep + ystep + 0];
    real c1101 = spline->c[n + zstep + ystep + 1];
    real c1102 = spline->c[n + zstep + ystep + 2];
    real c1103 = spline->c[n + zstep + ystep + 3];
    real c1104 = spline->c[n + zstep + ystep + 4];
    real c1105 = spline->c[n + zstep + ystep + 5];
    real c1106 = spline->c[n + zstep + ystep + 6];
    real c1107 = spline->c[n + zstep + ystep + 7];

    real c1110 = spline->c[n + zstep + ystep + xstep + 0];
    real c1111 = spline->c[n + zstep + ystep + xstep + 1];
    real c1112 = spline->c[n + zstep + ystep + xstep + 2];
    real c1113 = spline->c[n + zstep + ystep + xstep + 3];
    real c1114 = spline->c[n + zstep + ystep + xstep + 4];
    real c1115 = spline->c[n + zstep + ystep + xstep + 5];
    real c1116 = spline->c[n + zstep + ystep + xstep + 6];
    real c1117 = spline->c[n + zstep + ystep + xstep + 7];

    f_df[0] = (dzi * (dxi * (dyi * c0000 + dy * c0100) +
                      dx * (dyi * c0010 + dy * c0110)) +
               dz * (dxi * (dyi * c1000 + dy * c1100) +
                     dx * (dyi * c1010 + dy * c1110))) +
              xg2 / 6 *
                  (dzi * (dxi3 * (dyi * c0001 + dy * c0101) +
                          dx3 * (dyi * c0011 + dy * c0111)) +
                   dz * (dxi3 * (dyi * c1001 + dy * c1101) +
                         dx3 * (dyi * c1011 + dy * c1111))) +
              yg2 / 6 *
                  (dzi * (dxi * (dyi3 * c0002 + dy3 * c0102) +
                          dx * (dyi3 * c0012 + dy3 * c0112)) +
                   dz * (dxi * (dyi3 * c1002 + dy3 * c1102) +
                         dx * (dyi3 * c1012 + dy3 * c1112))) +
              zg2 / 6 *
                  (dzi3 * (dxi * (dyi * c0003 + dy * c0103) +
                           dx * (dyi * c0013 + dy * c0113)) +
                   dz3 * (dxi * (dyi * c1003 + dy * c1103) +
                          dx * (dyi * c1013 + dy * c1113))) +
              xg2 * yg2 / 36 *
                  (dzi * (dxi3 * (dyi3 * c0004 + dy3 * c0104) +
                          dx3 * (dyi3 * c0014 + dy3 * c0114)) +
                   dz * (dxi3 * (dyi3 * c1004 + dy3 * c1104) +
                         dx3 * (dyi3 * c1014 + dy3 * c1114))) +
              xg2 * zg2 / 36 *
                  (dzi3 * (dxi3 * (dyi * c0005 + dy * c0105) +
                           dx3 * (dyi * c0015 + dy * c0115)) +
                   dz3 * (dxi3 * (dyi * c1005 + dy * c1105) +
                          dx3 * (dyi * c1015 + dy * c1115))) +
              yg2 * zg2 / 36 *
                  (dzi3 * (dxi * (dyi3 * c0006 + dy3 * c0106) +
                           dx * (dyi3 * c0016 + dy3 * c0116)) +
                   dz3 * (dxi * (dyi3 * c1006 + dy3 * c1106) +
                          dx * (dyi3 * c1016 + dy3 * c1116))) +
              xg2 * yg2 * zg2 / 216 *
                  (dzi3 * (dxi3 * (dyi3 * c0007 + dy3 * c0107) +
                           dx3 * (dyi3 * c0017 + dy3 * c0117)) +
                   dz3 * (dxi3 * (dyi3 * c1007 + dy3 * c1107) +
                          dx3 * (dyi3 * c1017 + dy3 * c1117)));

    f_df[1] =
        xgi *
            (dzi * (-(dyi * c0000 + dy * c0100) + (dyi * c0010 + dy * c0110)) +
             dz * (-(dyi * c1000 + dy * c1100) + (dyi * c1010 + dy * c1110))) +
        xg / 6 *
            (dzi * (dxi3dx * (dyi * c0001 + dy * c0101) +
                    dx3dx * (dyi * c0011 + dy * c0111)) +
             dz * (dxi3dx * (dyi * c1001 + dy * c1101) +
                   dx3dx * (dyi * c1011 + dy * c1111))) +
        xgi * yg2 / 6 *
            (dzi * (-(dyi3 * c0002 + dy3 * c0102) +
                    (dyi3 * c0012 + dy3 * c0112)) +
             dz * (-(dyi3 * c1002 + dy3 * c1102) +
                   (dyi3 * c1012 + dy3 * c1112))) +
        xgi * zg2 / 6 *
            (dzi3 * (-(dyi * c0003 + dy * c0103) + (dyi * c0013 + dy * c0113)) +
             dz3 * (-(dyi * c1003 + dy * c1103) + (dyi * c1013 + dy * c1113))) +
        xg * yg2 / 36 *
            (dzi * (dxi3dx * (dyi3 * c0004 + dy3 * c0104) +
                    dx3dx * (dyi3 * c0014 + dy3 * c0114)) +
             dz * (dxi3dx * (dyi3 * c1004 + dy3 * c1104) +
                   dx3dx * (dyi3 * c1014 + dy3 * c1114))) +
        xg * zg2 / 36 *
            (dzi3 * (dxi3dx * (dyi * c0005 + dy * c0105) +
                     dx3dx * (dyi * c0015 + dy * c0115)) +
             dz3 * (dxi3dx * (dyi * c1005 + dy * c1105) +
                    dx3dx * (dyi * c1015 + dy * c1115))) +
        xgi * yg2 * zg2 / 36 *
            (dzi3 * (-(dyi3 * c0006 + dy3 * c0106) +
                     (dyi3 * c0016 + dy3 * c0116)) +
             dz3 * (-(dyi3 * c1006 + dy3 * c1106) +
                    (dyi3 * c1016 + dy3 * c1116))) +
        xg * yg2 * zg2 / 216 *
            (dzi3 * (dxi3dx * (dyi3 * c0007 + dy3 * c0107) +
                     dx3dx * (dyi3 * c0017 + dy3 * c0117)) +
             dz3 * (dxi3dx * (dyi3 * c1007 + dy3 * c1107) +
                    dx3dx * (dyi3 * c1017 + dy3 * c1117)));

    f_df[2] = ygi * (dzi * (dxi * (-c0000 + c0100) + dx * (-c0010 + c0110)) +
                     dz * (dxi * (-c1000 + c1100) + dx * (-c1010 + c1110))) +
              ygi * xg2 / 6 *
                  (dzi * (dxi3 * (-c0001 + c0101) + dx3 * (-c0011 + c0111)) +
                   dz * (dxi3 * (-c1001 + c1101) + dx3 * (-c1011 + c1111))) +
              yg / 6 *
                  (dzi * (dxi * (dyi3dy * c0002 + dy3dy * c0102) +
                          dx * (dyi3dy * c0012 + dy3dy * c0112)) +
                   dz * (dxi * (dyi3dy * c1002 + dy3dy * c1102) +
                         dx * (dyi3dy * c1012 + dy3dy * c1112))) +
              ygi * zg2 / 6 *
                  (dzi3 * (dxi * (-c0003 + c0103) + dx * (-c0013 + c0113)) +
                   dz3 * (dxi * (-c1003 + c1103) + dx * (-c1013 + c1113))) +
              xg2 * yg / 36 *
                  (dzi * (dxi3 * (dyi3dy * c0004 + dy3dy * c0104) +
                          dx3 * (dyi3dy * c0014 + dy3dy * c0114)) +
                   dz * (dxi3 * (dyi3dy * c1004 + dy3dy * c1104) +
                         dx3 * (dyi3dy * c1014 + dy3dy * c1114))) +
              ygi * xg2 * zg2 / 36 *
                  (dzi3 * (dxi3 * (-c0005 + c0105) + dx3 * (-c0015 + c0115)) +
                   dz3 * (dxi3 * (-c1005 + c1105) + dx3 * (-c1015 + c1115))) +
              yg * zg2 / 36 *
                  (dzi3 * (dxi * (dyi3dy * c0006 + dy3dy * c0106) +
                           dx * (dyi3dy * c0016 + dy3dy * c0116)) +
                   dz3 * (dxi * (dyi3dy * c1006 + dy3dy * c1106) +
                          dx * (dyi3dy * c1016 + dy3dy * c1116))) +
              xg2 * yg * zg2 / 216 *
                  (dzi3 * (dxi3 * (dyi3dy * c0007 + dy3dy * c0107) +
                           dx3 * (dyi3dy * c0017 + dy3dy * c0117)) +
                   dz3 * (dxi3 * (dyi3dy * c1007 + dy3dy * c1107) +
                          dx3 * (dyi3dy * c1017 + dy3dy * c1117)));

    f_df[3] = zgi * (-(dxi * (dyi * c0000 + dy * c0100) +
                       dx * (dyi * c0010 + dy * c0110)) +
                     (dxi * (dyi * c1000 + dy * c1100) +
                      dx * (dyi * c1010 + dy * c1110))) +
              xg2 * zgi / 6 *
                  (-(dxi3 * (dyi * c0001 + dy * c0101) +
                     dx3 * (dyi * c0011 + dy * c0111)) +
                   (dxi3 * (dyi * c1001 + dy * c1101) +
                    dx3 * (dyi * c1011 + dy * c1111))) +
              yg2 * zgi / 6 *
                  (-(dxi * (dyi3 * c0002 + dy3 * c0102) +
                     dx * (dyi3 * c0012 + dy3 * c0112)) +
                   (dxi * (dyi3 * c1002 + dy3 * c1102) +
                    dx * (dyi3 * c1012 + dy3 * c1112))) +
              zg / 6 *
                  (dzi3dz * (dxi * (dyi * c0003 + dy * c0103) +
                             dx * (dyi * c0013 + dy * c0113)) +
                   dz3dz * (dxi * (dyi * c1003 + dy * c1103) +
                            dx * (dyi * c1013 + dy * c1113))) +
              xg2 * yg2 * zgi / 36 *
                  (-(dxi3 * (dyi3 * c0004 + dy3 * c0104) +
                     dx3 * (dyi3 * c0014 + dy3 * c0114)) +
                   (dxi3 * (dyi3 * c1004 + dy3 * c1104) +
                    dx3 * (dyi3 * c1014 + dy3 * c1114))) +
              xg2 * zg / 36 *
                  (dzi3dz * (dxi3 * (dyi * c0005 + dy * c0105) +
                             dx3 * (dyi * c0015 + dy * c0115)) +
                   dz3dz * (dxi3 * (dyi * c1005 + dy * c1105) +
                            dx3 * (dyi * c1015 + dy * c1115))) +
              yg2 * zg / 36 *
                  (dzi3dz * (dxi * (dyi3 * c0006 + dy3 * c0106) +
                             dx * (dyi3 * c0016 + dy3 * c0116)) +
                   dz3dz * (dxi * (dyi3 * c1006 + dy3 * c1106) +
                            dx * (dyi3 * c1016 + dy3 * c1116))) +
              xg2 * yg2 * zg / 216 *
                  (dzi3dz * (dxi3 * (dyi3 * c0007 + dy3 * c0107) +
                             dx3 * (dyi3 * c0017 + dy3 * c0117)) +
                   dz3dz * (dxi3 * (dyi3 * c1007 + dy3 * c1107) +
                            dx3 * (dyi3 * c1017 + dy3 * c1117)));

    /* d2f/dx2 */
    f_df[4] = (dzi * (dxi * (dyi * c0001 + dy * c0101) +
                      dx * (dyi * c0011 + dy * c0111)) +
               dz * (dxi * (dyi * c1001 + dy * c1101) +
                     dx * (dyi * c1011 + dy * c1111))) +
              yg2 / 6 *
                  (dzi * (dxi * (dyi3 * c0004 + dy3 * c0104) +
                          dx * (dyi3 * c0014 + dy3 * c0114)) +
                   dz * (dxi * (dyi3 * c1004 + dy3 * c1104) +
                         dx * (dyi3 * c1014 + dy3 * c1114))) +
              zg2 / 6 *
                  (dzi3 * (dxi * (dyi * c0005 + dy * c0105) +
                           dx * (dyi * c0015 + dy * c0115)) +
                   dz3 * (dxi * (dyi * c1005 + dy * c1105) +
                          dx * (dyi * c1015 + dy * c1115))) +
              yg2 * zg2 / 36 *
                  (dzi3 * (dxi * (dyi3 * c0007 + dy3 * c0107) +
                           dx * (dyi3 * c0017 + dy3 * c0117)) +
                   dz3 * (dxi * (dyi3 * c1007 + dy3 * c1107) +
                          dx * (dyi3 * c1017 + dy3 * c1117)));

    f_df[5] = (dzi * (dxi * (dyi * c0002 + dy * c0102) +
                      dx * (dyi * c0012 + dy * c0112)) +
               dz * (dxi * (dyi * c1002 + dy * c1102) +
                     dx * (dyi * c1012 + dy * c1112))) +
              xg2 / 6 *
                  (dzi * (dxi3 * (dyi * c0004 + dy * c0104) +
                          dx3 * (dyi * c0014 + dy * c0114)) +
                   dz * (dxi3 * (dyi * c1004 + dy * c1104) +
                         dx3 * (dyi * c1014 + dy * c1114))) +
              zg2 / 6 *
                  (dzi3 * (dxi * (dyi * c0006 + dy * c0106) +
                           dx * (dyi * c0016 + dy * c0116)) +
                   dz3 * (dxi * (dyi * c1006 + dy * c1106) +
                          dx * (dyi * c1016 + dy * c1116))) +
              xg2 * zg2 / 36 *
                  (dzi3 * (dxi3 * (dyi * c0007 + dy * c0107) +
                           dx3 * (dyi * c0017 + dy * c0117)) +
                   dz3 * (dxi3 * (dyi * c1007 + dy * c1107) +
                          dx3 * (dyi * c1017 + dy * c1117)));

    f_df[6] = (dzi * (dxi * (dyi * c0003 + dy * c0103) +
                      dx * (dyi * c0013 + dy * c0113)) +
               dz * (dxi * (dyi * c1003 + dy * c1103) +
                     dx * (dyi * c1013 + dy * c1113))) +
              xg2 / 6 *
                  (dzi * (dxi3 * (dyi * c0005 + dy * c0105) +
                          dx3 * (dyi * c0015 + dy * c0115)) +
                   dz * (dxi3 * (dyi * c1005 + dy * c1105) +
                         dx3 * (dyi * c1015 + dy * c1115))) +
              yg2 / 6 *
                  (dzi * (dxi * (dyi3 * c0006 + dy3 * c0106) +
                          dx * (dyi3 * c0016 + dy3 * c0116)) +
                   dz * (dxi * (dyi3 * c1006 + dy3 * c1106) +
                         dx * (dyi3 * c1016 + dy3 * c1116))) +
              xg2 * yg2 / 36 *
                  (dzi * (dxi3 * (dyi3 * c0007 + dy3 * c0107) +
                          dx3 * (dyi3 * c0017 + dy3 * c0117)) +
                   dz * (dxi3 * (dyi3 * c1007 + dy3 * c1107) +
                         dx3 * (dyi3 * c1017 + dy3 * c1117)));

    f_df[7] =
        xgi * ygi *
            (dzi * ((c0000 - c0100) - (c0010 - c0110)) +
             dz * ((c1000 - c1100) - (c1010 - c1110))) +
        ygi * xg / 6 *
            (dzi * (dxi3dx * (-c0001 + c0101) + dx3dx * (-c0011 + c0111)) +
             dz * (dxi3dx * (-c1001 + c1101) + dx3dx * (-c1011 + c1111))) +
        xgi * yg / 6 *
            (dzi * (-(dyi3dy * c0002 + dy3dy * c0102) +
                    (dyi3dy * c0012 + dy3dy * c0112)) +
             dz * (-(dyi3dy * c1002 + dy3dy * c1102) +
                   (dyi3dy * c1012 + dy3dy * c1112))) +
        xgi * ygi * zg2 / 6 *
            (dzi3 * ((c0003 - c0103) - (c0013 - c0113)) +
             dz3 * ((c1003 - c1103) - (c1013 - c1113))) +
        xg * yg / 36 *
            (dzi * (dxi3dx * (dyi3dy * c0004 + dy3dy * c0104) +
                    dx3dx * (dyi3dy * c0014 + dy3dy * c0114)) +
             dz * (dxi3dx * (dyi3dy * c1004 + dy3dy * c1104) +
                   dx3dx * (dyi3dy * c1014 + dy3dy * c1114))) +
        ygi * xg * zg2 / 36 *
            (dzi3 * (dxi3dx * (-c0005 + c0105) + dx3dx * (-c0015 + c0115)) +
             dz3 * (dxi3dx * (-c1005 + c1105) + dx3dx * (-c1015 + c1115))) +
        xgi * yg * zg2 / 36 *
            (dzi3 * (-(dyi3dy * c0006 + dy3dy * c0106) +
                     (dyi3dy * c0016 + dy3dy * c0116)) +
             dz3 * (-(dyi3dy * c1006 + dy3dy * c1106) +
                    (dyi3dy * c1016 + dy3dy * c1116))) +
        xg * yg * zg2 / 216 *
            (dzi3 * (dxi3dx * (dyi3dy * c0007 + dy3dy * c0107) +
                     dx3dx * (dyi3dy * c0017 + dy3dy * c0117)) +
             dz3 * (dxi3dx * (dyi3dy * c1007 + dy3dy * c1107) +
                    dx3dx * (dyi3dy * c1017 + dy3dy * c1117)));

    f_df[8] =
        xgi * zgi *
            (((dyi * c0000 + dy * c0100) - (dyi * c0010 + dy * c0110)) -
             ((dyi * c1000 + dy * c1100) - (dyi * c1010 + dy * c1110))) +
        xg * zgi / 6 *
            (-(dxi3dx * (dyi * c0001 + dy * c0101) +
               dx3dx * (dyi * c0011 + dy * c0111)) +
             (dxi3dx * (dyi * c1001 + dy * c1101) +
              dx3dx * (dyi * c1011 + dy * c1111))) +
        xgi * yg2 * zgi / 6 *
            (((dyi3 * c0002 + dy3 * c0102) - (dyi3 * c0012 + dy3 * c0112)) -
             ((dyi3 * c1002 + dy3 * c1102) - (dyi3 * c1012 + dy3 * c1112))) +
        xgi * zg / 6 *
            (dzi3dz *
                 (-(dyi * c0003 + dy * c0103) + (dyi * c0013 + dy * c0113)) +
             dz3dz *
                 (-(dyi * c1003 + dy * c1103) + (dyi * c1013 + dy * c1113))) +
        xg * yg2 * zgi / 36 *
            (-(dxi3dx * (dyi3 * c0004 + dy3 * c0104) +
               dx3dx * (dyi3 * c0014 + dy3 * c0114)) +
             (dxi3dx * (dyi3 * c1004 + dy3 * c1104) +
              dx3dx * (dyi3 * c1014 + dy3 * c1114))) +
        xg * zg / 36 *
            (dzi3dz * (dxi3dx * (dyi * c0005 + dy * c0105) +
                       dx3dx * (dyi * c0015 + dy * c0115)) +
             dz3dz * (dxi3dx * (dyi * c1005 + dy * c1105) +
                      dx3dx * (dyi * c1015 + dy * c1115))) +
        xgi * yg2 * zg / 36 *
            (dzi3dz * (-(dyi3 * c0006 + dy3 * c0106) +
                       (dyi3 * c0016 + dy3 * c0116)) +
             dz3dz * (-(dyi3 * c1006 + dy3 * c1106) +
                      (dyi3 * c1016 + dy3 * c1116))) +
        xg * yg2 * zg / 216 *
            (dzi3dz * (dxi3dx * (dyi3 * c0007 + dy3 * c0107) +
                       dx3dx * (dyi3 * c0017 + dy3 * c0117)) +
             dz3dz * (dxi3dx * (dyi3 * c1007 + dy3 * c1107) +
                      dx3dx * (dyi3 * c1017 + dy3 * c1117)));

    f_df[9] = ygi * zgi *
                  ((dxi * (c0000 - c0100) + dx * (c0010 - c0110)) -
                   (dxi * (c1000 - c1100) + dx * (c1010 - c1110))) +
              ygi * xg2 * zgi / 6 *
                  ((dxi3 * (c0001 - c0101) + dx3 * (c0011 - c0111)) -
                   (dxi3 * (c1001 - c1101) + dx3 * (c1011 - c1111))) +
              yg * zgi / 6 *
                  (-(dxi * (dyi3dy * c0002 + dy3dy * c0102) +
                     dx * (dyi3dy * c0012 + dy3dy * c0112)) +
                   (dxi * (dyi3dy * c1002 + dy3dy * c1102) +
                    dx * (dyi3dy * c1012 + dy3dy * c1112))) +
              ygi * zg / 6 *
                  (dzi3dz * (dxi * (-c0003 + c0103) + dx * (-c0013 + c0113)) +
                   dz3dz * (dxi * (-c1003 + c1103) + dx * (-c1013 + c1113))) +
              xg2 * yg * zgi / 36 *
                  (-(dxi3 * (dyi3dy * c0004 + dy3dy * c0104) +
                     dx3 * (dyi3dy * c0014 + dy3dy * c0114)) +
                   (dxi3 * (dyi3dy * c1004 + dy3dy * c1104) +
                    dx3 * (dyi3dy * c1014 + dy3dy * c1114))) +
              ygi * xg2 * zg / 36 *
                  (dzi3dz * (dxi3 * (-c0005 + c0105) + dx3 * (-c0015 + c0115)) +
                   dz3dz * (dxi3 * (-c1005 + c1105) + dx3 * (-c1015 + c1115))) +
              yg * zg / 36 *
                  (dzi3dz * (dxi * (dyi3dy * c0006 + dy3dy * c0106) +
                             dx * (dyi3dy * c0016 + dy3dy * c0116)) +
                   dz3dz * (dxi * (dyi3dy * c1006 + dy3dy * c1106) +
                            dx * (dyi3dy * c1016 + dy3dy * c1116))) +
              xg2 * yg * zg / 216 *
                  (dzi3dz * (dxi3 * (dyi3dy * c0007 + dy3dy * c0107) +
                             dx3 * (dyi3dy * c0017 + dy3dy * c0117)) +
                   dz3dz * (dxi3 * (dyi3dy * c1007 + dy3dy * c1107) +
                            dx3 * (dyi3dy * c1017 + dy3dy * c1117)));
    return err;
}
