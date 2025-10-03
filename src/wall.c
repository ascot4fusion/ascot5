/**
 * Implements the wall interface (see wall.h).
 */
#include "wall.h"
#include "ascot5.h"
#include "wall/wall_2d.h"
#include "wall/wall_3d.h"
#include <math.h>
#include <stdio.h>

void wall_free(wall_data *wall)
{
    switch (wall->type)
    {
    case wall_contour2d:
        WallContour2D_free(wall->contour2d);
        break;
    case wall_triangular3d:
        WallTriangular3D_free(wall->triangular3d);
        break;
    }
}

void wall_offload(wall_data *wall)
{
    switch (wall->type)
    {
    case wall_contour2d:
        WallContour2D_offload(wall->contour2d);
        break;
    case wall_triangular3d:
        WallTriangular3D_offload(wall->triangular3d);
        break;
    }
}

int wall_hit_wall(
    real r1, real phi1, real z1, real r2, real phi2, real z2, wall_data *wall,
    real *w_coll)
{
    int ret = 0;
    switch (wall->type)
    {
    case wall_contour2d:
        ret = WallContour2D_hit_wall(r1, z1, r2, z2, w_coll, wall->contour2d);
        break;
    case wall_triangular3d:
        ret = WallTriangular3D_hit_wall(
            r1, phi1, z1, r2, phi2, z2, w_coll, wall->triangular3d);
        break;
    }
    return ret;
}

int wall_get_n_elements(wall_data *wall)
{
    int ret = 0;
    switch (wall->type)
    {
    case wall_contour2d:
        ret = wall->contour2d->n;
        break;
    case wall_triangular3d:
        ret = wall->triangular3d->n;
        break;
    }
    return ret;
}

int wall_get_flag(wall_data *wall, int idx)
{
    int flag = 0;
    switch (wall->type)
    {
    case wall_contour2d:
        flag = wall->contour2d->flag[idx];
        break;
    case wall_triangular3d:
        flag = wall->triangular3d->flag[idx];
        break;
    }
    return flag;
}