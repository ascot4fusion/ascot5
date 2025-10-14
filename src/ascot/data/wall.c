/**
 * Implements the wall interface (see wall.h).
 */
#include "wall.h"
#include "defines.h"
#include "wall_contour2d.h"
#include "wall_triangular3d.h"
#include <math.h>
#include <stdio.h>

void Wall_free(Wall *wall)
{
    switch (wall->type)
    {
    case WALL_CONTOUR2D:
        WallContour2D_free(wall->contour2d);
        break;
    case WALL_TRIANGULAR3D:
        WallTriangular3D_free(wall->triangular3d);
        break;
    }
}

void Wall_offload(Wall *wall)
{
    switch (wall->type)
    {
    case WALL_CONTOUR2D:
        WallContour2D_offload(wall->contour2d);
        break;
    case WALL_TRIANGULAR3D:
        WallTriangular3D_offload(wall->triangular3d);
        break;
    }
}

size_t Wall_eval_intersection(
    real w_coll[1], real r1, real phi1, real z1, real r2, real phi2, real z2,
    Wall *wall)
{
    size_t ret = 0;
    switch (wall->type)
    {
    case WALL_CONTOUR2D:
        ret = WallContour2D_eval_intersection(
            w_coll, r1, z1, r2, z2, wall->contour2d);
        break;
    case WALL_TRIANGULAR3D:
        ret = WallTriangular3D_eval_intersection(
            w_coll, r1, phi1, z1, r2, phi2, z2, wall->triangular3d);
        break;
    }
    return ret;
}

size_t Wall_get_n_elements(Wall *wall)
{
    size_t n = 0;
    switch (wall->type)
    {
    case WALL_CONTOUR2D:
        n = wall->contour2d->n;
        break;
    case WALL_TRIANGULAR3D:
        n = wall->triangular3d->n;
        break;
    }
    return n;
}

int Wall_get_flag(size_t idx, Wall *wall)
{
    int flag = 0;
    switch (wall->type)
    {
    case WALL_CONTOUR2D:
        flag = wall->contour2d->flag[idx + 1];
        break;
    case WALL_TRIANGULAR3D:
        flag = wall->triangular3d->flag[idx + 1];
        break;
    }
    return flag;
}