/**
 * Linear interpolation (see interp.h).
 */
#include "defines.h"
#include "interp.h"
#include "utils/mathlib.h"
#include <math.h>
#include <stdlib.h>

void linint1D_init(
    Linear1D *str, real *c, size_t n_x, int bc_x, real x_min, real x_max)
{
    real x_grid = (x_max - x_min) / (n_x - 1 * (bc_x == NATURALBC));

    str->c = c;
    str->n_x = n_x;
    str->bc_x = bc_x;
    str->x_min = x_min;
    str->x_max = x_max;
    str->x_grid = x_grid;
}

int linint1D_eval_f(real *f, Linear1D *str, real x)
{
    /* Make sure periodic coordinates are within [max, min] region. */
    if (str->bc_x == PERIODICBC)
    {
        x = fmod(x - str->x_min, str->x_max - str->x_min) + str->x_min;
        x = x + (x < str->x_min) * (str->x_max - str->x_min);
    }

    /* index for x variable */
    size_t i_x = (x - str->x_min) / str->x_grid;
    /**< Normalized x coordinate in current cell */
    real dx = (x - (str->x_min + i_x * str->x_grid)) / str->x_grid;

    size_t x1 = 1; /* Index jump one x forward */

    int err = 0;

    /* Enforce periodic BC or check that the coordinate is within the grid. */
    if (str->bc_x == PERIODICBC && i_x == str->n_x - 1)
    {
        x1 = -(str->n_x - 1) * x1;
    }
    else if (str->bc_x == NATURALBC && !(x >= str->x_min && x <= str->x_max))
    {
        err = 1;
    }

    if (!err)
    {
        *f = str->c[i_x] * (1 - dx) + str->c[i_x + x1] * dx;
    }
    return err;
}

void linint2D_init(
    Linear2D *str, real *c, size_t n_x, size_t n_y, int bc_x, int bc_y,
    real x_min, real x_max, real y_min, real y_max)
{
    real x_grid = (x_max - x_min) / (n_x - 1 * (bc_x == NATURALBC));
    real y_grid = (y_max - y_min) / (n_y - 1 * (bc_y == NATURALBC));

    str->c = c;
    str->n_x = n_x;
    str->n_y = n_y;
    str->x_min = x_min;
    str->x_max = x_max;
    str->y_min = y_min;
    str->y_max = y_max;
    str->x_grid = x_grid;
    str->y_grid = y_grid;
}

int linint2D_eval_f(real *f, Linear2D *str, real x, real y)
{
    real c00, c01, c10, c11;
    real c0, c1;

    /* Make sure periodic coordinates are within [max, min] region. */
    if (str->bc_x == PERIODICBC)
    {
        x = fmod(x - str->x_min, str->x_max - str->x_min) + str->x_min;
        x = x + (x < str->x_min) * (str->x_max - str->x_min);
    }
    if (str->bc_y == PERIODICBC)
    {
        y = fmod(y - str->y_min, str->y_max - str->y_min) + str->y_min;
        y = y + (y < str->y_min) * (str->y_max - str->y_min);
    }

    /* Index for x variable */
    size_t i_x = (x - str->x_min) / str->x_grid;
    /* Normalized x coordinate in current cell */
    real dx = (x - (str->x_min + i_x * str->x_grid)) / str->x_grid;

    /* Index for y variable */
    size_t i_y = (y - str->y_min) / str->y_grid;
    /* Normalized y coordinate in current cell */
    real dy = (y - (str->y_min + i_y * str->y_grid)) / str->y_grid;

    size_t n = i_y * str->n_x + i_x; /* Index jump to cell       */
    size_t x1 = 1;                   /* Index jump one x forward */
    size_t y1 = str->n_x;            /* Index jump one y forward */

    int err = 0;

    /* Enforce periodic BC or check that the coordinate is within the grid. */
    if (str->bc_x == PERIODICBC && i_x == str->n_x - 1)
    {
        x1 = -(str->n_x - 1) * x1;
    }
    else if (str->bc_x == NATURALBC && !(x >= str->x_min && x <= str->x_max))
    {
        err = 1;
    }
    if (str->bc_y == PERIODICBC && i_y == str->n_y - 1)
    {
        y1 = -(str->n_y - 1) * y1;
    }
    else if (str->bc_y == NATURALBC && !(y >= str->y_min && y <= str->y_max))
    {
        err = 1;
    }

    if (!err)
    {
        /* Values at grid cell corners */
        c00 = str->c[n];
        c10 = str->c[n + x1];
        c01 = str->c[n + y1];
        c11 = str->c[n + y1 + x1];
        /* Interpolate along x */
        c0 = c00 * (1 - dx) + c10 * dx;
        c1 = c01 * (1 - dx) + c11 * dx;
        /* Finally we interpolate these values along y */
        *f = c0 * (1 - dy) + c1 * dy;
    }

    return err;
}

void linint3D_init(
    Linear3D *str, real *c, size_t n_x, size_t n_y, size_t n_z, int bc_x, int bc_y,
    int bc_z, real x_min, real x_max, real y_min, real y_max, real z_min,
    real z_max)
{
    real x_grid = (x_max - x_min) / (n_x - 1 * (bc_x == NATURALBC));
    real y_grid = (y_max - y_min) / (n_y - 1 * (bc_y == NATURALBC));
    real z_grid = (z_max - z_min) / (n_z - 1 * (bc_z == NATURALBC));

    str->c = c;
    str->n_x = n_x;
    str->n_y = n_y;
    str->n_z = n_z;
    str->bc_x = bc_x;
    str->bc_y = bc_y;
    str->bc_z = bc_z;
    str->x_min = x_min;
    str->x_max = x_max;
    str->y_min = y_min;
    str->y_max = y_max;
    str->z_min = z_min;
    str->z_max = z_max;
    str->x_grid = x_grid;
    str->y_grid = y_grid;
    str->z_grid = z_grid;
}

int linint3D_eval_f(real *f, Linear3D *str, real x, real y, real z)
{
    real c000, c100, c001, c101, c010, c110, c011, c111;
    real c00, c01, c10, c11;
    real c0, c1;
    /* Make sure periodic coordinates are within [max, min] region. */
    if (str->bc_x == PERIODICBC)
    {
        x = fmod(x - str->x_min, str->x_max - str->x_min) + str->x_min;
        x = x + (x < str->x_min) * (str->x_max - str->x_min);
    }
    if (str->bc_y == PERIODICBC)
    {
        y = fmod(y - str->y_min, str->y_max - str->y_min) + str->y_min;
        y = y + (y < str->y_min) * (str->y_max - str->y_min);
    }
    if (str->bc_z == PERIODICBC)
    {
        z = fmod(z - str->z_min, str->z_max - str->z_min) + str->z_min;
        z = z + (z < str->z_min) * (str->z_max - str->z_min);
    }

    /* index for x variable */
    size_t i_x = (x - str->x_min) / str->x_grid;
    /* Normalized x coordinate in current cell */
    real dx = (x - (str->x_min + i_x * str->x_grid)) / str->x_grid;

    /* index for y variable */
    size_t i_y = (y - str->y_min) / str->y_grid;
    /* Normalized y coordinate in current cell */
    real dy = (y - (str->y_min + i_y * str->y_grid)) / str->y_grid;

    /* index for z variable */
    size_t i_z = (z - str->z_min) / str->z_grid;
    /* Normalized z coordinate in current cell */
    real dz = (z - (str->z_min + i_z * str->z_grid)) / str->z_grid;

    /* Index jump to cell */
    size_t n = i_y * str->n_z * str->n_x + i_z * str->n_x + i_x;
    size_t x1 = 1;                   /* Index jump one x forward */
    size_t y1 = str->n_z * str->n_x; /* Index jump one y forward */
    size_t z1 = str->n_x;            /* Index jump one z forward */

    int err = 0;

    /* Enforce periodic BC or check that the coordinate is within the grid. */
    if (str->bc_x == PERIODICBC && i_x == str->n_x - 1)
    {
        x1 = -(str->n_x - 1) * x1;
    }
    else if (str->bc_x == NATURALBC && !(x >= str->x_min && x <= str->x_max))
    {
        err = 1;
    }
    if (str->bc_y == PERIODICBC && i_y == str->n_y - 1)
    {
        y1 = -(str->n_y - 1) * y1;
    }
    else if (str->bc_y == NATURALBC && !(y >= str->y_min && y <= str->y_max))
    {
        err = 1;
    }
    if (str->bc_z == PERIODICBC && i_z == str->n_z - 1)
    {
        z1 = -(str->n_z - 1) * z1;
    }
    else if (str->bc_z == NATURALBC && !(z >= str->z_min && z <= str->z_max))
    {
        err = 1;
    }

    if (!err)
    {
        /* Values at grid cell corners */
        c000 = str->c[n];
        c100 = str->c[n + x1];
        c001 = str->c[n + y1];
        c101 = str->c[n + y1 + x1];
        c010 = str->c[n + z1];
        c110 = str->c[n + z1 + x1];
        c011 = str->c[n + y1 + z1];
        c111 = str->c[n + y1 + z1 + x1];
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
