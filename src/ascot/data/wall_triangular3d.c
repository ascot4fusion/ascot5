/**
 * Implements wall_triangular3d.c.
 */
#include "wall_triangular3d.h"
#include "defines.h"
#include "utils/list.h"
#include "utils/mathlib.h"
#include "utils/octree.h"
#include "wall.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int WallTriangular3D_init(
    WallTriangular3D *wall, size_t n, real *vertices, int *flag)
{

    wall->n = n;
    wall->vertices = (real *)malloc(9 * n * sizeof(real));
    real xmin = vertices[0], xmax = vertices[0];
    real ymin = vertices[1], ymax = vertices[1];
    real zmin = vertices[2], zmax = vertices[2];
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            wall->vertices[i * 9 + j * 3 + 0] = vertices[i * 9 + j * 3 + 0];
            wall->vertices[i * 9 + j * 3 + 1] = vertices[i * 9 + j * 3 + 1];
            wall->vertices[i * 9 + j * 3 + 2] = vertices[i * 9 + j * 3 + 2];

            /* Find min & max values of the volume occupied by the wall */
            xmin = fmin(xmin, vertices[i * 9 + j * 3 + 0]);
            xmax = fmax(xmax, vertices[i * 9 + j * 3 + 0]);
            ymin = fmin(ymin, vertices[i * 9 + j * 3 + 1]);
            ymax = fmax(ymax, vertices[i * 9 + j * 3 + 1]);
            zmin = fmin(zmin, vertices[i * 9 + j * 3 + 2]);
            zmax = fmax(zmax, vertices[i * 9 + j * 3 + 2]);
        }
    }

    /* Add a little bit of padding so we don't need to worry about triangles
       clipping the edges */
    wall->xmin = xmin - 0.1;
    wall->xmax = xmax + 0.1;
    wall->ymin = ymin - 0.1;
    wall->ymax = ymax + 0.1;
    wall->zmin = zmin - 0.1;
    wall->zmax = zmax + 0.1;

    /* Depth of the octree in which the triangles are sorted */
    wall->depth = WALL_OCTREE_DEPTH;
    size_t ngrid = 1;
    for (size_t i = 0; i < wall->depth - 1; i++)
    {
        ngrid *= 2;
    }
    wall->ngrid = ngrid;

    wall->dx = (wall->xmax - wall->xmin) / wall->ngrid;
    wall->dy = (wall->ymax - wall->ymin) / wall->ngrid;
    wall->dz = (wall->zmax - wall->zmin) / wall->ngrid;

    wall->depth = WALL_OCTREE_DEPTH;
    WallTriangular3D_init_octree(wall);
    wall->flag = (int *)malloc(n * sizeof(int));
    memcpy(wall->flag, flag, wall->n * sizeof(int));

    return 0;
}

void WallTriangular3D_free(WallTriangular3D *wall)
{
    free(wall->tree_array);
    free(wall->vertices);
}

void WallTriangular3D_offload(WallTriangular3D *wall)
{
    (void)wall;
    GPU_MAP_TO_DEVICE(
        data->wall_flag [0:data->n], data->vertices [0:data->n * 9],
        data->tree_array [0:data->tree_array_size])
}

size_t WallTriangular3D_hit_wall(
    real r1, real phi1, real z1, real r2, real phi2, real z2, real *w_coll,
    WallTriangular3D *wall)
{
    real rpz1[3], rpz2[3];
    rpz1[0] = r1;
    rpz1[1] = phi1;
    rpz1[2] = z1;
    rpz2[0] = r2;
    rpz2[1] = phi2;
    rpz2[2] = z2;

    real q1[3], q2[3];
    math_rpz2xyz(rpz1, q1);
    math_rpz2xyz(rpz2, q2);

    ptrdiff_t ix1 = (int)floor(
        (q1[0] - wall->xmin) / ((wall->xmax - wall->xmin) / (wall->ngrid)));
    ptrdiff_t iy1 = (int)floor(
        (q1[1] - wall->ymin) / ((wall->ymax - wall->ymin) / (wall->ngrid)));
    ptrdiff_t iz1 = (int)floor(
        (q1[2] - wall->zmin) / ((wall->zmax - wall->zmin) / (wall->ngrid)));

    ptrdiff_t ix2 = (int)floor(
        (q2[0] - wall->xmin) / ((wall->xmax - wall->xmin) / (wall->ngrid)));
    ptrdiff_t iy2 = (int)floor(
        (q2[1] - wall->ymin) / ((wall->ymax - wall->ymin) / (wall->ngrid)));
    ptrdiff_t iz2 = (int)floor(
        (q2[2] - wall->zmin) / ((wall->zmax - wall->zmin) / (wall->ngrid)));

    size_t hit_tri = 0;
    real smallest_w = 1.1;

    for (size_t i = 0; i <= llabs(ix2 - ix1); i++)
    {
        for (size_t j = 0; j <= llabs(iy2 - iy1); j++)
        {
            for (size_t k = 0; k <= llabs(iz2 - iz1); k++)
            {
                ptrdiff_t ix = ix1 + i * ((int)copysign(1, ix2 - ix1));
                ptrdiff_t iy = iy1 + j * ((int)copysign(1, iy2 - iy1));
                ptrdiff_t iz = iz1 + k * ((int)copysign(1, iz2 - iz1));

                if (ix >= 0 && ix < wall->ngrid && iy >= 0 &&
                    iy < wall->ngrid && iz >= 0 && iz < wall->ngrid)
                {

                    size_t ilist = wall->tree_array
                                    [ix * wall->ngrid * wall->ngrid +
                                     iy * wall->ngrid + iz];

                    for (size_t l = 0; l < wall->tree_array[ilist]; l++)
                    {
                        size_t itri = wall->tree_array[ilist + l + 1];
                        real w = octree_tri_collision(
                            q1, q2, &wall->vertices[9 * itri],
                            &wall->vertices[9 * itri + 3],
                            &wall->vertices[9 * itri + 6]);
                        if (w >= 0 && w < smallest_w)
                        {
                            smallest_w = w;
                            hit_tri = itri + 1;
                        }
                    }
                }
            }
        }
    }
    *w_coll = smallest_w;
    return hit_tri;
}

void WallTriangular3D_init_octree(WallTriangular3D *wall)
{

    /* Construct the octree and store triangles there */
    octree_node *tree;
    octree_create(
        &tree, wall->xmin, wall->xmax, wall->ymin, wall->ymax, wall->zmin,
        wall->zmax, wall->depth);
    for (size_t i = 0; i < wall->n; i++)
    {
        real t1[3], t2[3], t3[3];
        t1[0] = wall->vertices[i * 9];
        t1[1] = wall->vertices[i * 9 + 1];
        t1[2] = wall->vertices[i * 9 + 2];
        t2[0] = wall->vertices[i * 9 + 3];
        t2[1] = wall->vertices[i * 9 + 4];
        t2[2] = wall->vertices[i * 9 + 5];
        t3[0] = wall->vertices[i * 9 + 6];
        t3[1] = wall->vertices[i * 9 + 7];
        t3[2] = wall->vertices[i * 9 + 8];
        octree_add(tree, t1, t2, t3, i);
    }

    /* Create lists for triangles in each grid square and fill the lists
       by querying the octree in each grid point */
    size_t ncell = wall->ngrid * wall->ngrid * wall->ngrid;
    list_int_node **tri_list =
        (list_int_node **)malloc(ncell * sizeof(list_int_node *));
    for (size_t ix = 0; ix < wall->ngrid; ix++)
    {
        for (size_t iy = 0; iy < wall->ngrid; iy++)
        {
            for (size_t iz = 0; iz < wall->ngrid; iz++)
            {
                real p[3];
                p[0] = wall->xmin + ix * wall->dx + 0.5 * wall->dx;
                p[1] = wall->ymin + iy * wall->dy + 0.5 * wall->dy;
                p[2] = wall->zmin + iz * wall->dz + 0.5 * wall->dz;

                int cell_index =
                    ix * wall->ngrid * wall->ngrid + iy * wall->ngrid + iz;
                tri_list[cell_index] = octree_get(tree, p);
            }
        }
    }

    /* construct an array from the triangle lists */
    size_t list_size = 0;
    for (size_t i = 0; i < ncell; i++)
    {
        list_size += list_int_size(tri_list[i]);
    }
    wall->tree_array = (size_t *)malloc((2 * ncell + list_size) * sizeof(size_t));
    wall->n_tree_array = 2 * ncell + list_size;

    size_t next_empty_list = ncell;
    for (size_t i = 0; i < ncell; i++)
    {
        /* First ncell elements store the position where the actual cell data
         * begins in tree_array */
        wall->tree_array[i] = next_empty_list;

        /* The first data point in the actual cell data is the number of
         * triangles in this cell */
        wall->tree_array[next_empty_list] = list_int_size(tri_list[i]);

        /* Store triangle IDs that are located in this cell */
        for (size_t j = 0; j < wall->tree_array[next_empty_list]; j++)
        {
            wall->tree_array[next_empty_list + j + 1] =
                list_int_get(tri_list[i], j);
        }
        next_empty_list += wall->tree_array[next_empty_list] + 1;
    }
    free(tri_list);
    octree_free(&tree);
}

void WallTriangular3D_init_tree(WallTriangular3D *wall, real *offload_array)
{
    /* create a list for holding the triangle ids in each cell */
    size_t ncell = wall->ngrid * wall->ngrid * wall->ngrid;
    list_int_node **tri_list =
        (list_int_node **)malloc(ncell * sizeof(list_int_node *));
    for (size_t i = 0; i < ncell; i++)
    {
        list_int_create(&tri_list[i]);
    }

    /* iterate through all triangles and cells and fill the lists */
    for (size_t i = 0; i < wall->n; i++)
    {
        real t1[3], t2[3], t3[3];
        t1[0] = offload_array[i * 9];
        t1[1] = offload_array[i * 9 + 1];
        t1[2] = offload_array[i * 9 + 2];
        t2[0] = offload_array[i * 9 + 3];
        t2[1] = offload_array[i * 9 + 4];
        t2[2] = offload_array[i * 9 + 5];
        t3[0] = offload_array[i * 9 + 6];
        t3[1] = offload_array[i * 9 + 7];
        t3[2] = offload_array[i * 9 + 8];

        #pragma omp parallel for
        for (size_t ix = 0; ix < wall->ngrid; ix++)
        {
            for (size_t iy = 0; iy < wall->ngrid; iy++)
            {
                for (size_t iz = 0; iz < wall->ngrid; iz++)
                {
                    real c1[3], c2[3];
                    real epsilon = 1e-6;
                    c1[0] = wall->xmin + ix * wall->dx - epsilon;
                    c2[0] = wall->xmin + (ix + 1) * wall->dx + epsilon;
                    c1[1] = wall->ymin + iy * wall->dy - epsilon;
                    c2[1] = wall->ymin + (iy + 1) * wall->dy + epsilon;
                    c1[2] = wall->zmin + iz * wall->dz - epsilon;
                    c2[2] = wall->zmin + (iz + 1) * wall->dz + epsilon;
                    int result = octree_tri_in_cube(t1, t2, t3, c1, c2);
                    size_t cell_index =
                        ix * wall->ngrid * wall->ngrid + iy * wall->ngrid + iz;
                    if (result > 0)
                    {
                        list_int_add(tri_list[cell_index], i);
                    }
                }
            }
        }
    }

    /* construct an array from the triangle lists */
    size_t list_size = 0;
    for (size_t i = 0; i < ncell; i++)
    {
        list_size += list_int_size(tri_list[i]);
    }

    wall->tree_array = (size_t *)malloc((2 * ncell + list_size) * sizeof(size_t));
    wall->n_tree_array = 2 * ncell + list_size;

    size_t next_empty_list = ncell;
    for (size_t i = 0; i < ncell; i++)
    {
        wall->tree_array[i] = next_empty_list;
        wall->tree_array[next_empty_list] = list_int_size(tri_list[i]);
        for (size_t j = 0; j < wall->tree_array[next_empty_list]; j++)
        {
            wall->tree_array[next_empty_list + j + 1] =
                list_int_get(tri_list[i], j);
        }
        next_empty_list += wall->tree_array[next_empty_list] + 1;
        list_int_free(&tri_list[i]);
    }
    free(tri_list);
}

size_t WallTriangular3D_hit_wall_full(
    real r1, real phi1, real z1, real r2, real phi2, real z2, real *w_coll,
    WallTriangular3D *wall)
{
    real rpz1[3], rpz2[3];
    rpz1[0] = r1;
    rpz1[1] = phi1;
    rpz1[2] = z1;
    rpz2[0] = r2;
    rpz2[1] = phi2;
    rpz2[2] = z2;

    real q1[3], q2[3];
    math_rpz2xyz(rpz1, q1);
    math_rpz2xyz(rpz2, q2);

    size_t hit_tri = 0;
    real smallest_w = 1.1;
    real w;

    for (size_t j = 0; j < wall->n; j++)
    {
        w = octree_tri_collision(
            q1, q2, &wall->vertices[9 * j], &wall->vertices[9 * j + 3],
            &wall->vertices[9 * j + 6]);
        if (w > 0)
        {
            if (w < smallest_w)
            {
                smallest_w = w;
                hit_tri = j + 1;
            }
        }
    }

    *w_coll = smallest_w;
    return hit_tri;
}
