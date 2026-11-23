/**
 * Linear interpolation (see "interp.h").
 */
#include "defines.h"
#include "interp.h"
#include "utils/mathlib.h"
#include <math.h>
#include <stdlib.h>

void Linear1D_init(
    Linear1D *linear, size_t nx, int xbc, real xlim[2], real c[nx])
{
    linear->c = c;
    linear->nx = nx;
    linear->xbc = xbc;
    linear->xlim[0] = xlim[0];
    linear->xlim[1] = xlim[1];
    linear->dx = (xlim[1] - xlim[0]) / (nx - 1 * (xbc == NATURALBC));
}

void Linear2D_init(
    Linear2D *linear, size_t nx, size_t ny, int xbc, int ybc, real xlim[2],
    real ylim[2], real c[nx * ny])
{
    linear->c = c;
    linear->nx = nx;
    linear->ny = ny;
    linear->xbc = xbc;
    linear->ybc = ybc;
    linear->xlim[0] = xlim[0];
    linear->xlim[1] = xlim[1];
    linear->ylim[0] = ylim[0];
    linear->ylim[1] = ylim[1];
    linear->dx = (xlim[1] - xlim[0]) / (nx - 1 * (xbc == NATURALBC));
    linear->dy = (ylim[1] - ylim[0]) / (ny - 1 * (ybc == NATURALBC));
}

void Linear3D_init(
    Linear3D *linear, size_t nx, size_t ny, size_t nz, int xbc, int ybc,
    int zbc, real xlim[2], real ylim[2], real zlim[2], real c[nx * ny * nz])
{
    linear->c = c;
    linear->nx = nx;
    linear->ny = ny;
    linear->nz = nz;
    linear->xbc = xbc;
    linear->ybc = ybc;
    linear->zbc = zbc;
    linear->xlim[0] = xlim[0];
    linear->xlim[1] = xlim[1];
    linear->ylim[0] = ylim[0];
    linear->ylim[1] = ylim[1];
    linear->zlim[0] = zlim[0];
    linear->zlim[1] = zlim[1];
    linear->dx = (xlim[1] - xlim[0]) / (nx - 1 * (xbc == NATURALBC));
    linear->dy = (ylim[1] - ylim[0]) / (ny - 1 * (ybc == NATURALBC));
    linear->dz = (zlim[1] - zlim[0]) / (nz - 1 * (zbc == NATURALBC));
}

int Linear1D_eval_f(real f[1], Linear1D *linear, real x)
{
    /* Make sure periodic coordinates are within [max, min] region. */
    if (linear->xbc == PERIODICBC)
    {
        x = fmod(x - linear->xlim[0], linear->xlim[1] - linear->xlim[0]) +
            linear->xlim[0];
        x = x + (x < linear->xlim[0]) * (linear->xlim[1] - linear->xlim[0]);
    }

    /* Cell index and normalized coordinate within the cell */
    real dx;
    size_t i_x = interp_cell(
        &dx, x, linear->xlim[0], linear->dx, linear->nx, linear->xbc);

    size_t x1 = 1; /* Index jump one x forward */

    int err = 0;

    /* Enforce periodic BC or check that the coordinate is within the grid. */
    if (linear->xbc == PERIODICBC && i_x == linear->nx - 1)
    {
        x1 = -(linear->nx - 1) * x1;
    }
    else if (
        linear->xbc == NATURALBC &&
        !(x >= linear->xlim[0] && x <= linear->xlim[1]))
    {
        err = 1;
    }

    if (!err)
    {
        *f = linear->c[i_x] * (1 - dx) + linear->c[i_x + x1] * dx;
    }
    return err;
}

int Linear2D_eval_f(real f[1], Linear2D *linear, real x, real y)
{
    real c00, c01, c10, c11;
    real c0, c1;

    /* Make sure periodic coordinates are within [max, min] region. */
    if (linear->xbc == PERIODICBC)
    {
        x = fmod(x - linear->xlim[0], linear->xlim[1] - linear->xlim[0]) +
            linear->xlim[0];
        x = x + (x < linear->xlim[0]) * (linear->xlim[1] - linear->xlim[0]);
    }
    if (linear->ybc == PERIODICBC)
    {
        y = fmod(y - linear->ylim[0], linear->ylim[1] - linear->ylim[0]) +
            linear->ylim[0];
        y = y + (y < linear->ylim[0]) * (linear->ylim[1] - linear->ylim[0]);
    }

    /* Cell index and normalized coordinate within the cell */
    real dx;
    size_t i_x = interp_cell(
        &dx, x, linear->xlim[0], linear->dx, linear->nx, linear->xbc);

    /* Same for y */
    real dy;
    size_t i_y = interp_cell(
        &dy, y, linear->ylim[0], linear->dy, linear->ny, linear->ybc);

    size_t n = i_x * linear->ny + i_y; /* Index jump to cell       */
    size_t x1 = linear->ny;            /* Index jump one x forward */
    size_t y1 = 1;                     /* Index jump one y forward */

    int err = 0;

    /* Enforce periodic BC or check that the coordinate is within the grid. */
    if (linear->xbc == PERIODICBC && i_x == linear->nx - 1)
    {
        x1 = -(linear->nx - 1) * x1;
    }
    else if (
        linear->xbc == NATURALBC &&
        !(x >= linear->xlim[0] && x <= linear->xlim[1]))
    {
        err = 1;
    }
    if (linear->ybc == PERIODICBC && i_y == linear->ny - 1)
    {
        y1 = -(linear->ny - 1) * y1;
    }
    else if (
        linear->ybc == NATURALBC &&
        !(y >= linear->ylim[0] && y <= linear->ylim[1]))
    {
        err = 1;
    }

    if (!err)
    {
        /* Values at grid cell corners */
        c00 = linear->c[n];
        c10 = linear->c[n + x1];
        c01 = linear->c[n + y1];
        c11 = linear->c[n + y1 + x1];
        /* Interpolate along x */
        c0 = c00 * (1 - dx) + c10 * dx;
        c1 = c01 * (1 - dx) + c11 * dx;
        /* Finally we interpolate these values along y */
        *f = c0 * (1 - dy) + c1 * dy;
    }

    return err;
}

int Linear3D_eval_f(real f[1], Linear3D *linear, real x, real y, real z)
{
    real c000, c100, c001, c101, c010, c110, c011, c111;
    real c00, c01, c10, c11;
    real c0, c1;
    /* Make sure periodic coordinates are within [max, min] region. */
    if (linear->xbc == PERIODICBC)
    {
        x = fmod(x - linear->xlim[0], linear->xlim[1] - linear->xlim[0]) +
            linear->xlim[0];
        x = x + (x < linear->xlim[0]) * (linear->xlim[1] - linear->xlim[0]);
    }
    if (linear->ybc == PERIODICBC)
    {
        y = fmod(y - linear->ylim[0], linear->ylim[1] - linear->ylim[0]) +
            linear->ylim[0];
        y = y + (y < linear->ylim[0]) * (linear->ylim[1] - linear->ylim[0]);
    }
    if (linear->zbc == PERIODICBC)
    {
        z = fmod(z - linear->zlim[0], linear->zlim[1] - linear->zlim[0]) +
            linear->zlim[0];
        z = z + (z < linear->zlim[0]) * (linear->zlim[1] - linear->zlim[0]);
    }

    real dx;
    size_t i_x = interp_cell(
        &dx, x, linear->xlim[0], linear->dx, linear->nx, linear->xbc);

    real dy;
    size_t i_y = interp_cell(
        &dy, y, linear->ylim[0], linear->dy, linear->ny, linear->ybc);

    real dz;
    size_t i_z = interp_cell(
        &dz, z, linear->zlim[0], linear->dz, linear->nz, linear->zbc);

    /* Index jump to cell */
    size_t n = i_x * linear->ny * linear->nz + i_y * linear->nz + i_z;
    size_t x1 = linear->ny * linear->nz; /* Index jump one x forward */
    size_t y1 = linear->nz;              /* Index jump one y forward */
    size_t z1 = 1;                       /* Index jump one z forward */

    int err = 0;

    /* Enforce periodic BC or check that the coordinate is within the grid. */
    if (linear->xbc == PERIODICBC && i_x == linear->nx - 1)
    {
        x1 = -(linear->nx - 1) * x1;
    }
    else if (
        linear->xbc == NATURALBC &&
        !(x >= linear->xlim[0] && x <= linear->xlim[1]))
    {
        err = 1;
    }
    if (linear->ybc == PERIODICBC && i_y == linear->ny - 1)
    {
        y1 = -(linear->ny - 1) * y1;
    }
    else if (
        linear->ybc == NATURALBC &&
        !(y >= linear->ylim[0] && y <= linear->ylim[1]))
    {
        err = 1;
    }
    if (linear->zbc == PERIODICBC && i_z == linear->nz - 1)
    {
        z1 = -(linear->nz - 1) * z1;
    }
    else if (
        linear->zbc == NATURALBC &&
        !(z >= linear->zlim[0] && z <= linear->zlim[1]))
    {
        err = 1;
    }

    if (!err)
    {
        /* Values at grid cell corners */
        c000 = linear->c[n];
        c100 = linear->c[n + x1];
        c001 = linear->c[n + y1];
        c101 = linear->c[n + y1 + x1];
        c010 = linear->c[n + z1];
        c110 = linear->c[n + z1 + x1];
        c011 = linear->c[n + y1 + z1];
        c111 = linear->c[n + y1 + z1 + x1];
        /* Interpolate along x */
        c00 = c000 * (1 - dx) + c100 * dx;
        c01 = c001 * (1 - dx) + c101 * dx;
        c10 = c010 * (1 - dx) + c110 * dx;
        c11 = c011 * (1 - dx) + c111 * dx;
        /* Interpolate these values along z */
        c0 = c00 * (1 - dz) + c10 * dz;
        c1 = c01 * (1 - dz) + c11 * dz;
        /* Finally we interpolate these values along y */
        *f = c0 * (1 - dy) + c1 * dy;
    }

    return err;
}

int Linear2D_eval_f_df(real f_df[3], Linear2D *linear, real x, real y)
{
    real c00, c01, c10, c11;
    real c0, c1;

    /* Make sure periodic coordinates are within [max, min] region. */
    if (linear->xbc == PERIODICBC)
    {
        x = fmod(x - linear->xlim[0], linear->xlim[1] - linear->xlim[0]) +
            linear->xlim[0];
        x = x + (x < linear->xlim[0]) * (linear->xlim[1] - linear->xlim[0]);
    }
    if (linear->ybc == PERIODICBC)
    {
        y = fmod(y - linear->ylim[0], linear->ylim[1] - linear->ylim[0]) +
            linear->ylim[0];
        y = y + (y < linear->ylim[0]) * (linear->ylim[1] - linear->ylim[0]);
    }

    /* Cell index and normalized coordinate within the cell */
    real dx;
    size_t i_x = interp_cell(
        &dx, x, linear->xlim[0], linear->dx, linear->nx, linear->xbc);

    /* Same for y */
    real dy;
    size_t i_y = interp_cell(
        &dy, y, linear->ylim[0], linear->dy, linear->ny, linear->ybc);

    size_t n = i_x * linear->ny + i_y; /* Index jump to cell       */
    size_t x1 = linear->ny;            /* Index jump one x forward */
    size_t y1 = 1;                     /* Index jump one y forward */

    int err = 0;

    /* Enforce periodic BC or check that the coordinate is within the grid. */
    if (linear->xbc == PERIODICBC && i_x == linear->nx - 1)
    {
        x1 = -(linear->nx - 1) * x1;
    }
    else if (
        linear->xbc == NATURALBC &&
        !(x >= linear->xlim[0] && x <= linear->xlim[1]))
    {
        err = 1;
    }
    if (linear->ybc == PERIODICBC && i_y == linear->ny - 1)
    {
        y1 = -(linear->ny - 1) * y1;
    }
    else if (
        linear->ybc == NATURALBC &&
        !(y >= linear->ylim[0] && y <= linear->ylim[1]))
    {
        err = 1;
    }

    if (err)
    {
        return err;
    }

    /* Values at grid cell corners */
    c00 = linear->c[n];
    c10 = linear->c[n + x1];
    c01 = linear->c[n + y1];
    c11 = linear->c[n + y1 + x1];
    /* Interpolate along x */
    c0 = c00 * (1 - dx) + c10 * dx;
    c1 = c01 * (1 - dx) + c11 * dx;
    /* Finally we interpolate these values along y */
    f_df[0] = c0 * (1 - dy) + c1 * dy;

    real xgi = 1.0 / linear->dx;
    real ygi = 1.0 / linear->dy;

    /* f_x */
    f_df[1] = ((c10 - c00) * (1 - dy) + (c11 - c01) * dy) * xgi;
    /* f_y */
    f_df[2] = (c1 - c0) * ygi;

    return 0;
}

int Linear3D_eval_f_df(real f_df[4], Linear3D *linear, real x, real y, real z)
{
    real c000, c100, c001, c101, c010, c110, c011, c111;
    real c00, c01, c10, c11;
    real c0, c1;
    /* Make sure periodic coordinates are within [max, min] region. */
    if (linear->xbc == PERIODICBC)
    {
        x = fmod(x - linear->xlim[0], linear->xlim[1] - linear->xlim[0]) +
            linear->xlim[0];
        x = x + (x < linear->xlim[0]) * (linear->xlim[1] - linear->xlim[0]);
    }
    if (linear->ybc == PERIODICBC)
    {
        y = fmod(y - linear->ylim[0], linear->ylim[1] - linear->ylim[0]) +
            linear->ylim[0];
        y = y + (y < linear->ylim[0]) * (linear->ylim[1] - linear->ylim[0]);
    }
    if (linear->zbc == PERIODICBC)
    {
        z = fmod(z - linear->zlim[0], linear->zlim[1] - linear->zlim[0]) +
            linear->zlim[0];
        z = z + (z < linear->zlim[0]) * (linear->zlim[1] - linear->zlim[0]);
    }

    real dx;
    size_t i_x = interp_cell(
        &dx, x, linear->xlim[0], linear->dx, linear->nx, linear->xbc);

    real dy;
    size_t i_y = interp_cell(
        &dy, y, linear->ylim[0], linear->dy, linear->ny, linear->ybc);

    real dz;
    size_t i_z = interp_cell(
        &dz, z, linear->zlim[0], linear->dz, linear->nz, linear->zbc);

    /* Index jump to cell */
    size_t n = i_x * linear->ny * linear->nz + i_y * linear->nz + i_z;
    size_t x1 = linear->ny * linear->nz; /* Index jump one x forward */
    size_t y1 = linear->nz;              /* Index jump one y forward */
    size_t z1 = 1;                       /* Index jump one z forward */

    int err = 0;

    /* Enforce periodic BC or check that the coordinate is within the grid. */
    if (linear->xbc == PERIODICBC && i_x == linear->nx - 1)
    {
        x1 = -(linear->nx - 1) * x1;
    }
    else if (
        linear->xbc == NATURALBC &&
        !(x >= linear->xlim[0] && x <= linear->xlim[1]))
    {
        err = 1;
    }
    if (linear->ybc == PERIODICBC && i_y == linear->ny - 1)
    {
        y1 = -(linear->ny - 1) * y1;
    }
    else if (
        linear->ybc == NATURALBC &&
        !(y >= linear->ylim[0] && y <= linear->ylim[1]))
    {
        err = 1;
    }
    if (linear->zbc == PERIODICBC && i_z == linear->nz - 1)
    {
        z1 = -(linear->nz - 1) * z1;
    }
    else if (
        linear->zbc == NATURALBC &&
        !(z >= linear->zlim[0] && z <= linear->zlim[1]))
    {
        err = 1;
    }

    if (err)
    {
        return err;
    }

    /* Values at grid cell corners */
    c000 = linear->c[n];
    c100 = linear->c[n + x1];
    c001 = linear->c[n + y1];
    c101 = linear->c[n + y1 + x1];
    c010 = linear->c[n + z1];
    c110 = linear->c[n + z1 + x1];
    c011 = linear->c[n + y1 + z1];
    c111 = linear->c[n + y1 + z1 + x1];
    /* Interpolate along x */
    c00 = c000 * (1 - dx) + c100 * dx;
    c01 = c001 * (1 - dx) + c101 * dx;
    c10 = c010 * (1 - dx) + c110 * dx;
    c11 = c011 * (1 - dx) + c111 * dx;
    /* Interpolate these values along z */
    c0 = c00 * (1 - dz) + c10 * dz;
    c1 = c01 * (1 - dz) + c11 * dz;
    /* Finally we interpolate these values along y */
    f_df[0] = c0 * (1 - dy) + c1 * dy;

    real xgi = 1.0 / linear->dx;
    real ygi = 1.0 / linear->dy;
    real zgi = 1.0 / linear->dz;

    /* f_x */
    f_df[1] = (((c100 - c000) * (1 - dz) + (c110 - c010) * dz) * (1 - dy) +
               ((c101 - c001) * (1 - dz) + (c111 - c011) * dz) * dy) *
              xgi;
    /* f_y */
    f_df[2] = (c1 - c0) * ygi;
    /* f_z */
    f_df[3] = ((c10 - c00) * (1 - dy) + (c11 - c01) * dy) * zgi;

    return 0;
}
