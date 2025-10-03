/**
 * Implements wall_3d.c.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ascot5.h"
#include "math.h"
#include "list.h"
#include "octree.h"
#include "wall_3d.h"


int WallTriangular3D_init(
    WallTriangular3D* wall, int n, real* vertices, int* flag
) {

    wall->n = n;
    wall->vertices = (real*)malloc(9 * n * sizeof(real));
    real xmin = vertices[0], xmax = vertices[0];
    real ymin = vertices[1], ymax = vertices[1];
    real zmin = vertices[2], zmax = vertices[2];
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < 3; j++) {
            wall->vertices[i*9 + j*3 + 0] = vertices[i*9 + j*3 + 0];
            wall->vertices[i*9 + j*3 + 1] = vertices[i*9 + j*3 + 1];
            wall->vertices[i*9 + j*3 + 2] = vertices[i*9 + j*3 + 2];

            /* Find min & max values of the volume occupied by the wall */
            xmin = fmin( xmin, vertices[i*9 + j*3 + 0] );
            xmax = fmax( xmax, vertices[i*9 + j*3 + 0] );
            ymin = fmin( ymin, vertices[i*9 + j*3 + 1] );
            ymax = fmax( ymax, vertices[i*9 + j*3 + 1] );
            zmin = fmin( zmin, vertices[i*9 + j*3 + 2] );
            zmax = fmax( zmax, vertices[i*9 + j*3 + 2] );
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
    int ngrid = 1;
    for(int i = 0; i < wall->depth - 1; i++) {
        ngrid *= 2;
    }
    wall->ngrid = ngrid;

    wall->dx = (wall->xmax - wall->xmin) / wall->ngrid;
    wall->dy = (wall->ymax - wall->ymin) / wall->ngrid;
    wall->dz = (wall->zmax - wall->zmin) / wall->ngrid;

    wall->depth = WALL_OCTREE_DEPTH;
    WallTriangular3D_init_octree(wall);
    wall->flag = (int*)malloc(n * sizeof(int));
    memcpy(wall->flag, flag, wall->n * sizeof(int));

    return 0;
}


void WallTriangular3D_free(WallTriangular3D* wall) {
    free(wall->tree_array);
    free(wall->vertices);
}


void WallTriangular3D_offload(WallTriangular3D* wall) {
    (void)wall;
    GPU_MAP_TO_DEVICE(
        data->wall_flag[0:data->n],
        data->vertices[0:data->n*9],
        data->tree_array[0:data->tree_array_size]
    )
}


int WallTriangular3D_hit_wall(
    real r1, real phi1, real z1, real r2, real phi2, real z2, real* w_coll,
    WallTriangular3D* wall
) {
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

    int ix1 = (int) floor((q1[0] - wall->xmin)
                          / ((wall->xmax - wall->xmin) / (wall->ngrid)));
    int iy1 = (int) floor((q1[1] - wall->ymin)
                          / ((wall->ymax - wall->ymin) / (wall->ngrid)));
    int iz1 = (int) floor((q1[2] - wall->zmin)
                          / ((wall->zmax - wall->zmin) / (wall->ngrid)));

    int ix2 = (int) floor((q2[0] - wall->xmin)
                          / ((wall->xmax - wall->xmin) / (wall->ngrid)));
    int iy2 = (int) floor((q2[1] - wall->ymin)
                          / ((wall->ymax - wall->ymin) / (wall->ngrid)));
    int iz2 = (int) floor((q2[2] - wall->zmin)
                          / ((wall->zmax - wall->zmin) / (wall->ngrid)));

    int hit_tri = 0;
    real smallest_w = 1.1;

    for(int i = 0; i <= abs(ix2-ix1); i++) {
        for(int j = 0; j <= abs(iy2-iy1); j++) {
            for(int k = 0; k <= abs(iz2-iz1); k++) {
                int ix = ix1 + i*((int) copysign(1, ix2-ix1));
                int iy = iy1 + j*((int) copysign(1, iy2-iy1));
                int iz = iz1 + k*((int) copysign(1, iz2-iz1));

                if(ix >= 0 && ix < wall->ngrid && iy >= 0 && iy < wall->ngrid
                   && iz >= 0 && iz < wall->ngrid) {

                    int ilist = wall->tree_array[ix*wall->ngrid*wall->ngrid
                                                  + iy*wall->ngrid + iz];

                    for(int l = 0; l < wall->tree_array[ilist]; l++) {
                        int itri = wall->tree_array[ilist+l+1];
                        real w = WallTriangular3D_tri_collision(
                            q1, q2,
                            &wall->vertices[9*itri],
                            &wall->vertices[9*itri+3],
                            &wall->vertices[9*itri+6]);
                        if(w >= 0 && w < smallest_w) {
                            smallest_w = w;
                            hit_tri = itri+1;
                        }
                    }
                }
            }
        }
    }
    *w_coll = smallest_w;
    return hit_tri;
}


int WallTriangular3D_tri_in_cube(
    real t1[3], real t2[3], real t3[3], real bb1[3], real bb2[3]
) {
    /* check if any point is inside the cube */
    if(bb1[0] <= t1[0] && t1[0] <= bb2[0]
        && bb1[1] <= t1[1] && t1[1] <= bb2[1]
        && bb1[2] <= t1[2] && t1[2] <= bb2[2])
        return 1;
    if(bb1[0] < t2[0] && t2[0] <= bb2[0]
        && bb1[1] <= t2[1] && t2[1] <= bb2[1]
        && bb1[2] <= t2[2] && t2[2] <= bb2[2])
        return 1;
    if(bb1[0] <= t3[0] && t3[0] <= bb2[0]
        && bb1[1] <= t3[1] && t3[1] <= bb2[1]
        && bb1[2] <= t3[2] && t3[2] <= bb2[2])
        return 1;

    /* no such luck; check if any of the cube edges intersects the triangle */
    real c000[3]; c000[0] = bb1[0]; c000[1] = bb1[1]; c000[2] = bb1[2];
    real c100[3]; c100[0] = bb2[0]; c100[1] = bb1[1]; c100[2] = bb1[2];
    real c010[3]; c010[0] = bb1[0]; c010[1] = bb2[1]; c010[2] = bb1[2];
    real c110[3]; c110[0] = bb2[0]; c110[1] = bb2[1]; c110[2] = bb1[2];
    real c001[3]; c001[0] = bb1[0]; c001[1] = bb1[1]; c001[2] = bb2[2];
    real c101[3]; c101[0] = bb2[0]; c101[1] = bb1[1]; c101[2] = bb2[2];
    real c011[3]; c011[0] = bb1[0]; c011[1] = bb2[1]; c011[2] = bb2[2];
    real c111[3]; c111[0] = bb2[0]; c111[1] = bb2[1]; c111[2] = bb2[2];

    if(   WallTriangular3D_tri_collision(c000, c100, t1, t2, t3) >= 0
       || WallTriangular3D_tri_collision(c000, c010, t1, t2, t3) >= 0
       || WallTriangular3D_tri_collision(c110, c010, t1, t2, t3) >= 0
       || WallTriangular3D_tri_collision(c110, c100, t1, t2, t3) >= 0
       || WallTriangular3D_tri_collision(c000, c001, t1, t2, t3) >= 0
       || WallTriangular3D_tri_collision(c010, c011, t1, t2, t3) >= 0
       || WallTriangular3D_tri_collision(c100, c101, t1, t2, t3) >= 0
       || WallTriangular3D_tri_collision(c110, c111, t1, t2, t3) >= 0
       || WallTriangular3D_tri_collision(c001, c101, t1, t2, t3) >= 0
       || WallTriangular3D_tri_collision(c001, c011, t1, t2, t3) >= 0
       || WallTriangular3D_tri_collision(c111, c011, t1, t2, t3) >= 0
       || WallTriangular3D_tri_collision(c111, c101, t1, t2, t3) >= 0)
        return 1;

    /* check for triangle edges intersecting cube quads */
    if(   WallTriangular3D_quad_collision(t1, t2, c000, c100, c110, c010) == 1
       || WallTriangular3D_quad_collision(t1, t2, c000, c010, c011, c001) == 1
       || WallTriangular3D_quad_collision(t1, t2, c000, c100, c101, c001) == 1
       || WallTriangular3D_quad_collision(t1, t2, c010, c110, c111, c011) == 1
       || WallTriangular3D_quad_collision(t1, t2, c100, c101, c111, c110) == 1
       || WallTriangular3D_quad_collision(t1, t2, c001, c101, c111, c011) == 1)
        return 1;
    if(   WallTriangular3D_quad_collision(t3, t2, c000, c100, c110, c010) == 1
       || WallTriangular3D_quad_collision(t3, t2, c000, c010, c011, c001) == 1
       || WallTriangular3D_quad_collision(t3, t2, c000, c100, c101, c001) == 1
       || WallTriangular3D_quad_collision(t3, t2, c010, c110, c111, c011) == 1
       || WallTriangular3D_quad_collision(t3, t2, c100, c101, c111, c110) == 1
       || WallTriangular3D_quad_collision(t3, t2, c001, c101, c111, c011) == 1)
        return 1;
    if(   WallTriangular3D_quad_collision(t1, t3, c000, c100, c110, c010) == 1
       || WallTriangular3D_quad_collision(t1, t3, c000, c010, c011, c001) == 1
       || WallTriangular3D_quad_collision(t1, t3, c000, c100, c101, c001) == 1
       || WallTriangular3D_quad_collision(t1, t3, c010, c110, c111, c011) == 1
       || WallTriangular3D_quad_collision(t1, t3, c100, c101, c111, c110) == 1
       || WallTriangular3D_quad_collision(t1, t3, c001, c101, c111, c011) == 1)
        return 1;
    return 0;
}


double WallTriangular3D_tri_collision(
    real q1[3], real q2[3], real t1[3], real t2[3], real t3[3]
) {
    real q12[3], Q12[3];
    Q12[0] = q2[0] - q1[0];
    Q12[1] = q2[1] - q1[1];
    Q12[2] = q2[2] - q1[2];
    math_unit(Q12, q12);

    real edge12[3];
    edge12[0] = t2[0] - t1[0];
    edge12[1] = t2[1] - t1[1];
    edge12[2] = t2[2] - t1[2];

    real edge13[3];
    edge13[0] = t3[0] - t1[0];
    edge13[1] = t3[1] - t1[1];
    edge13[2] = t3[2] - t1[2];

    real h[3];
    math_cross(q12, edge13, h);
    real det = math_dot(h, edge12);

    /* Check that the triangle has non-zero area */
    real normal[3], area;
    math_cross(edge12, edge13, normal);
    area = math_norm(normal);

    real w = -1.0;
    if( area > WALL_EPSILON ) {
        /* If ray is parallel to the triangle, nudge it a little bit so we don't
           have to handle annoying special cases */
        if( fabs(det) < WALL_EPSILON ) {
            Q12[0] = q2[0] - q1[0] + 2 * WALL_EPSILON * normal[0] / area;
            Q12[1] = q2[1] - q1[1] + 2 * WALL_EPSILON * normal[1] / area;
            Q12[2] = q2[2] - q1[2] + 2 * WALL_EPSILON * normal[2] / area;
            math_unit(Q12, q12);
            math_cross(q12, edge13, h);
            det = math_dot(h, edge12);
        }

        real tq11[3];
        tq11[0] = q1[0] - t1[0];
        tq11[1] = q1[1] - t1[1];
        tq11[2] = q1[2] - t1[2];

        real n[3];
        math_cross(tq11, edge12, n);

        real u = math_dot(h, tq11) / det;
        real v = math_dot(q12, n) / det;

        if( ( u >= 0.0 && u <= 1.0 ) && ( v >= 0.0 && u + v <= 1.0 )  ) {
            w = ( math_dot(n, edge13) / det ) / math_norm(Q12);
            if( w > 1.0 ) {
                w = -1.0;
            }
        }
    }

    return w;

}


int WallTriangular3D_quad_collision(
    real q1[3], real q2[3], real t1[3], real t2[3], real t3[3], real t4[3]
) {
    if(WallTriangular3D_tri_collision(q1, q2, t1, t2, t3) >= 0
       || WallTriangular3D_tri_collision(q1, q2, t1, t3, t4) >= 0)
        return 1;
    else
        return 0;
}


void WallTriangular3D_init_octree(WallTriangular3D* wall) {

    /* Construct the octree and store triangles there */
    octree_node* tree;
    octree_create(&tree, wall->xmin, wall->xmax, wall->ymin, wall->ymax, wall->zmin, wall->zmax,
                  wall->depth);
    int i;
    for(i = 0; i < wall->n; i++) {
        real t1[3], t2[3], t3[3];
        t1[0] = wall->vertices[i*9];
        t1[1] = wall->vertices[i*9+1];
        t1[2] = wall->vertices[i*9+2];
        t2[0] = wall->vertices[i*9+3];
        t2[1] = wall->vertices[i*9+4];
        t2[2] = wall->vertices[i*9+5];
        t3[0] = wall->vertices[i*9+6];
        t3[1] = wall->vertices[i*9+7];
        t3[2] = wall->vertices[i*9+8];
        octree_add(tree, t1, t2, t3, i);
    }

    /* Create lists for triangles in each grid square and fill the lists
       by querying the octree in each grid point */
    int ncell = wall->ngrid*wall->ngrid*wall->ngrid;
    list_int_node** tri_list =
        (list_int_node**) malloc(ncell * sizeof(list_int_node*));
    int ix, iy, iz;
    for(ix = 0; ix < wall->ngrid; ix++) {
        for(iy = 0; iy < wall->ngrid; iy++) {
            for(iz = 0; iz < wall->ngrid; iz++) {
                real p[3];
                p[0] = wall->xmin + ix * wall->dx + 0.5*wall->dx;
                p[1] = wall->ymin + iy * wall->dy + 0.5*wall->dy;
                p[2] = wall->zmin + iz * wall->dz + 0.5*wall->dz;

                int cell_index = ix*wall->ngrid*wall->ngrid+iy*wall->ngrid+iz;
                tri_list[cell_index] = octree_get(tree, p);
            }
        }
    }

    /* construct an array from the triangle lists */
    int list_size = 0;
    for(i = 0; i < ncell; i++) {
        list_size += list_int_size(tri_list[i]);
    }
    wall->tree_array = (int*) malloc((2*ncell + list_size)*sizeof(int));
    wall->tree_array_size = 2*ncell + list_size;

    int next_empty_list = ncell;
    for(i = 0; i < ncell; i++) {
        /* First ncell elements store the position where the actual cell data
         * begins in tree_array */
        wall->tree_array[i] = next_empty_list;

        /* The first data point in the actual cell data is the number of
         * triangles in this cell */
        wall->tree_array[next_empty_list] = list_int_size(tri_list[i]);

        /* Store triangle IDs that are located in this cell */
        for(int j = 0; j < wall->tree_array[next_empty_list]; j++) {
            wall->tree_array[next_empty_list+j+1] = list_int_get(tri_list[i], j);
        }
        next_empty_list += wall->tree_array[next_empty_list] + 1;
    }
    free(tri_list);
    octree_free(&tree);
}


void WallTriangular3D_init_tree(WallTriangular3D* wall, real* offload_array) {
    /* create a list for holding the triangle ids in each cell */
    int ncell = wall->ngrid * wall->ngrid * wall->ngrid;
    list_int_node** tri_list = (list_int_node**) malloc(ncell
                                                    *sizeof(list_int_node*));
    int i;
    for(i = 0; i < ncell; i++) {
        list_int_create(&tri_list[i]);
    }

    /* iterate through all triangles and cells and fill the lists */
    for(i = 0; i < wall->n; i++) {
        real t1[3], t2[3], t3[3];
        t1[0] = offload_array[i*9];
        t1[1] = offload_array[i*9+1];
        t1[2] = offload_array[i*9+2];
        t2[0] = offload_array[i*9+3];
        t2[1] = offload_array[i*9+4];
        t2[2] = offload_array[i*9+5];
        t3[0] = offload_array[i*9+6];
        t3[1] = offload_array[i*9+7];
        t3[2] = offload_array[i*9+8];

        int ix, iy, iz;
        #pragma omp parallel for private(ix, iy, iz)
        for(ix = 0; ix < wall->ngrid; ix++) {
            for(iy = 0; iy < wall->ngrid; iy++) {
                for(iz = 0; iz < wall->ngrid; iz++) {
                    real c1[3], c2[3];
                    real epsilon = 1e-6;
                    c1[0] = wall->xmin + ix * wall->dx - epsilon;
                    c2[0] = wall->xmin + (ix+1) * wall->dx + epsilon;
                    c1[1] = wall->ymin + iy * wall->dy - epsilon;
                    c2[1] = wall->ymin + (iy+1) * wall->dy + epsilon;
                    c1[2] = wall->zmin + iz * wall->dz - epsilon;
                    c2[2] = wall->zmin + (iz+1) * wall->dz + epsilon;
                    int result = WallTriangular3D_tri_in_cube(t1, t2, t3, c1, c2);
                    int cell_index = ix*wall->ngrid*wall->ngrid+iy*wall->ngrid+iz;
                    if(result > 0) {
                        list_int_add(tri_list[cell_index], i);
                    }
                }
            }
        }
    }

    /* construct an array from the triangle lists */
    int list_size = 0;
    for(i = 0; i < ncell; i++) {
        list_size += list_int_size(tri_list[i]);
    }

    wall->tree_array = (int*) malloc((2*ncell + list_size)*sizeof(int));
    wall->tree_array_size = 2*ncell + list_size;

    int next_empty_list = ncell;
    for(i = 0; i < ncell; i++) {
        wall->tree_array[i] = next_empty_list;
        wall->tree_array[next_empty_list] = list_int_size(tri_list[i]);
        int j;
        for(j = 0; j < wall->tree_array[next_empty_list]; j++) {
            wall->tree_array[next_empty_list+j+1] = list_int_get(tri_list[i], j);
        }
        next_empty_list += wall->tree_array[next_empty_list] + 1;
        list_int_free(&tri_list[i]);
    }
    free(tri_list);
}


int WallTriangular3D_hit_wall_full(
    real r1, real phi1, real z1, real r2, real phi2, real z2, real* w_coll,
    WallTriangular3D* wall
) {
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

    int hit_tri = 0;
    real smallest_w = 1.1;
    real w;
    int j;

    for(j = 0; j < wall->n; j++) {
        w = WallTriangular3D_tri_collision(q1, q2, &wall->vertices[9*j],
                &wall->vertices[9*j+3], &wall->vertices[9*j+6]);
        if(w > 0) {
            if(w < smallest_w) {
                smallest_w = w;
                hit_tri = j+1;
            }
        }
    }

    *w_coll = smallest_w;
    return hit_tri;
}
