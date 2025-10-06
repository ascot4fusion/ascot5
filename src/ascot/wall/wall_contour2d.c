/**
 * Implements wall_contour2d.h.
 */
#include "wall_contour2d.h"
#include "defines.h"
#include "wall.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int WallContour2D_init(
    WallContour2D *wall, int n, real r[n], real z[n], int flag[n])
{

    wall->n = n;
    wall->flag = (int *)malloc(n * sizeof(int));
    wall->r = (real *)malloc(n * sizeof(real));
    wall->z = (real *)malloc(n * sizeof(real));
    real rmin = r[0], rmax = r[0];
    real zmin = z[0], zmax = z[0];
    for (int i = 0; i < n; i++)
    {
        rmin = fmin(rmin, r[i]);
        rmax = fmax(rmax, r[i]);
        zmin = fmin(zmin, z[i]);
        zmax = fmax(zmax, z[i]);
        wall->r[i] = r[i];
        wall->z[i] = z[i];
        wall->flag[i] = flag[i];
    }
    return 0;
}

void WallContour2D_free(WallContour2D *wall)
{
    free(wall->r);
    free(wall->z);
}

void WallContour2D_offload(WallContour2D *wall)
{
    (void)wall;
    GPU_MAP_TO_DEVICE(
        wall->r [0:wall->n], wall->z [0:wall->n], wall->flag [0:wall->n])
}

int WallContour2D_inside(real r, real z, WallContour2D *wall)
{
    int hits = 0;
    for (int i = 0; i < wall->n; i++)
    {
        real wr1, wr2, wz1, wz2;
        if (i == wall->n - 1)
        {
            wz1 = wall->z[i] - z;
            wz2 = wall->z[0] - z;
            wr1 = wall->r[i] - r;
            wr2 = wall->r[0] - r;
        }
        else
        {
            wz1 = wall->z[i] - z;
            wr1 = wall->r[i] - r;
            wz2 = wall->z[i + 1] - z;
            wr2 = wall->r[i + 1] - r;
        }
        if (wz1 * wz2 < 0)
        {
            real ri = wr1 + (wz1 * (wr2 - wr1)) / (wz1 - wz2);
            if (ri > 0)
            {
                hits++;
            }
        }
    }
    return hits % 2;
}

int WallContour2D_hit_wall(
    real r1, real z1, real r2, real z2, real *w_coll, WallContour2D *wall)
{
    return WallContour2D_find_intersection(r1, z1, r2, z2, w_coll, wall);
}

int WallContour2D_find_intersection(
    real r1, real z1, real r2, real z2, real *w_coll, WallContour2D *wall)
{
    int tile = 0;
    real t0 = 2.0; // Helper variable to pick the closest intersection
    for (int i = 0; i < wall->n; i++)
    {
        real r3, z3, r4, z4;
        if (i == wall->n - 1)
        {
            r3 = wall->r[i];
            z3 = wall->z[i];
            r4 = wall->r[0];
            z4 = wall->z[0];
        }
        else
        {
            r3 = wall->r[i];
            z3 = wall->z[i];
            r4 = wall->r[i + 1];
            z4 = wall->z[i + 1];
        }

        real div = (r1 - r2) * (z3 - z4) - (z1 - z2) * (r3 - r4);
        real t = ((r1 - r3) * (z3 - z4) - (z1 - z3) * (r3 - r4)) / div;
        real u = ((r1 - r3) * (z1 - z2) - (z1 - z3) * (r1 - r2)) / div;
        if (0 <= t && t <= 1.0 && 0 <= u && u <= 1.0 && t < t0)
        {
            t0 = t;
            tile = i + 1;
        }
    }
    *w_coll = t0;
    return tile;
}
