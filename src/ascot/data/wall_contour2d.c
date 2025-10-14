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
    WallContour2D *wall, size_t n, real r[n], real z[n], int flag[n])
{
    wall->n = n;
    wall->flag = (int *)malloc(n * sizeof(int));
    wall->r = (real *)malloc(n * sizeof(real));
    wall->z = (real *)malloc(n * sizeof(real));
    real rmin = r[0], rmax = r[0];
    real zmin = z[0], zmax = z[0];
    for (size_t i = 0; i < n; i++)
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

size_t WallContour2D_eval_intersection(
    real w_coll[1], real r1, real z1, real r2, real z2, WallContour2D *wall)
{
    size_t tile = 0;
    w_coll[0] = 1.0;
    for (size_t i = 0; i < wall->n; i++)
    {
        real r3, z3, r4, z4;
        size_t next_index = i == wall->n - 1 ? 0 : i + 1;
        r3 = wall->r[i];
        z3 = wall->z[i];
        r4 = wall->r[next_index];
        z4 = wall->r[next_index];

        real div = (r1 - r2) * (z3 - z4) - (z1 - z2) * (r3 - r4);
        real t = ((r1 - r3) * (z3 - z4) - (z1 - z3) * (r3 - r4)) / div;
        real u = ((r1 - r3) * (z1 - z2) - (z1 - z3) * (r1 - r2)) / div;

        int closer_hit =
            0 <= t && t <= 1.0 && 0 <= u && u <= 1.0 && t < w_coll[0];
        tile = closer_hit ? i + 1 : tile;
        w_coll[0] = closer_hit ? t : w_coll[0];
    }
    return tile;
}
