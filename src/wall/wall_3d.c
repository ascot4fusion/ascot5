/**
 * @file wall_3d.c
 * @brief 3D wall collision checks
 *
 * 3D wall model where wall consists of small triangles that form a surface
 * mesh. A wall collision happens when a marker intersects a wall triangle's
 * plane. During initialization, the computational model is divided into
 * smaller cell using octree, and wall triangles are divided according to
 * which cell(s) they inhabit. Collision checks are only made with respect to
 * triangles that are in the same cell as the marker.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../ascot5.h"
#include "wall_3d.h"
#include "../math.h"
#include "../list.h"
#include "../octree.h"
#include "../print.h"

/**
 * @brief Initialize 3D wall data and check inputs
 *
 * Before calling this function, the offload struct is expected to hold
 * triangle positions as
 *
 * [x1_1, y1_1, z1_1, x2_1, y2_1, z2_1, x3_1, y3_1, z3_1,... ],
 *
 * where first index is for the triangle vertex and second for the triangle
 * itself.
 *
 * This function fill rest of the offload struct and constructs the octree,
 * while also printing some values as sanity check.

 * The default octree depth is defined by macro WALL_OCTREE_DEPTH in wall_3d.h.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array
 * @param int_offload_array pointer to offload array containing integers
 *
 * @return zero if initialization succeeded
 */
int wall_3d_init_offload(wall_3d_offload_data* offload_data,
                         real** offload_array, int** int_offload_array) {

    /* Find min & max values of the volume occupied by the wall triangles. */
    real xmin = (*offload_array)[0], xmax = (*offload_array)[0];
    real ymin = (*offload_array)[1], ymax = (*offload_array)[1];
    real zmin = (*offload_array)[2], zmax = (*offload_array)[2];
    for(int i=0; i<offload_data->n*3; i++) {
        xmin = fmin( xmin, (*offload_array)[i*3 + 0] );
        xmax = fmax( xmax, (*offload_array)[i*3 + 0] );
        ymin = fmin( ymin, (*offload_array)[i*3 + 1] );
        ymax = fmax( ymax, (*offload_array)[i*3 + 1] );
        zmin = fmin( zmin, (*offload_array)[i*3 + 2] );
        zmax = fmax( zmax, (*offload_array)[i*3 + 2] );
    }

    /* Add a little bit of padding so we don't need to worry about triangles
       clipping the edges */
    offload_data->xmin = xmin - 0.1;
    offload_data->xmax = xmax + 0.1;
    offload_data->ymin = ymin - 0.1;
    offload_data->ymax = ymax + 0.1;
    offload_data->zmin = zmin - 0.1;
    offload_data->zmax = zmax + 0.1;

    /* Depth of the octree in which the triangles are sorted */
    offload_data->depth = WALL_OCTREE_DEPTH;
    int ngrid = 1;
    for(int i = 0; i < offload_data->depth - 1; i++) {
        ngrid *= 2;
    }
    offload_data->ngrid = ngrid;

    offload_data->xgrid = (offload_data->xmax - offload_data->xmin)
                          / offload_data->ngrid;
    offload_data->ygrid = (offload_data->ymax - offload_data->ymin)
                          / offload_data->ngrid;
    offload_data->zgrid = (offload_data->zmax - offload_data->zmin)
                          / offload_data->ngrid;

    print_out(VERBOSE_IO, "\n3D wall model (wall_3D)\n");
    print_out(VERBOSE_IO,
              "Number of wall elements %d\n"
              "Spanning xmin = %2.3f m, xmax = %2.3f m\n"
              "         ymin = %2.3f m, ymax = %2.3f m\n"
              "         zmin = %2.3f m, zmax = %2.3f m\n",
              offload_data->n,
              offload_data->xmin, offload_data->xmax, offload_data->ymin,
              offload_data->ymax, offload_data->zmin, offload_data->zmax);

    /* Store flags temporarily */
    int* flag = (int*) malloc(offload_data->n * sizeof(int));
    memcpy(flag, *int_offload_array, offload_data->n * sizeof(int));
    free(*int_offload_array);

    /* Initialize octree. Note that this function reserves space for flags */
    wall_3d_init_octree(offload_data, *offload_array, int_offload_array);

    /* Copy flags to offload array */
    int idx = offload_data->offload_array_length - offload_data->n;
    memcpy(&(*int_offload_array)[idx], flag, offload_data->n * sizeof(int));
    free(flag);

    return 0;
}

/**
 * @brief Free offload array and reset parameters
 *
 * This function deallocates the offload_array.
 *
 * This function is host only.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array
 * @param int_offload_array pointer to offload array containing integers
 */
void wall_3d_free_offload(wall_3d_offload_data* offload_data,
                          real** offload_array, int** int_offload_array) {
    free(*offload_array);
    *offload_array = NULL;

    free(*int_offload_array);
    *int_offload_array = NULL;
}

/**
 * @brief Initialize wall data struct on target
 *
 * This function copies the wall parameters from the offload struct to the
 * struct on target and sets the wall data pointers to correct offsets in the
 * offload array.
 *
 * @param w pointer to data struct on target
 * @param offload_data pointer to offload data struct
 * @param offload_array offload array
 * @param int_offload_array offload array containing integers
 */
void wall_3d_init(wall_3d_data* w, wall_3d_offload_data* offload_data,
                  real* offload_array, int* int_offload_array) {
    w->n = offload_data->n;
    w->xmin  = offload_data->xmin;
    w->xmax  = offload_data->xmax;
    w->xgrid = offload_data->xgrid;
    w->ymin  = offload_data->ymin;
    w->ymax  = offload_data->ymax;
    w->ygrid = offload_data->ygrid;
    w->zmin  = offload_data->zmin;
    w->zmax  = offload_data->zmax;
    w->zgrid = offload_data->zgrid;
    w->depth = offload_data->depth;
    w->ngrid = offload_data->ngrid;
    w->wall_tris = &offload_array[0];

    w->tree_array_size = offload_data->int_offload_array_length - w->n;
    w->tree_array = &int_offload_array[0];
    w->flag       = &int_offload_array[w->tree_array_size];
}

/**
 * @brief Construct wall octree iteratively
 *
 * Constructs the octree array by iterating through all wall triangles and
 * octree grid to identify triangles belonging to each grid cell.
 *
 * Slow, only for testing purposes.
 *
 * @param w pointer to wall data
 * @param offload_array offload array
 */
void wall_3d_init_tree(wall_3d_data* w, real* offload_array) {
    /* create a list for holding the triangle ids in each cell */
    int ncell = w->ngrid*w->ngrid*w->ngrid;
    list_int_node** tri_list = (list_int_node**) malloc(ncell
                                                    *sizeof(list_int_node*));
    int i;
    for(i = 0; i < ncell; i++) {
        list_int_create(&tri_list[i]);
    }

    /* iterate through all triangles and cells and fill the lists */
    for(i = 0; i < w->n; i++) {
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
        for(ix = 0; ix < w->ngrid; ix++) {
            for(iy = 0; iy < w->ngrid; iy++) {
                for(iz = 0; iz < w->ngrid; iz++) {
                    real c1[3], c2[3];
                    real epsilon = 1e-6;
                    c1[0] = w->xmin + ix * w->xgrid - epsilon;
                    c2[0] = w->xmin + (ix+1) * w->xgrid + epsilon;
                    c1[1] = w->ymin + iy * w->ygrid - epsilon;
                    c2[1] = w->ymin + (iy+1) * w->ygrid + epsilon;
                    c1[2] = w->zmin + iz * w->zgrid - epsilon;
                    c2[2] = w->zmin + (iz+1) * w->zgrid + epsilon;
                    int result = wall_3d_tri_in_cube(t1, t2, t3, c1, c2);
                    int cell_index = ix*w->ngrid*w->ngrid+iy*w->ngrid+iz;
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

    w->tree_array_size = 2*ncell + list_size;
    w->tree_array = (int*) malloc((w->tree_array_size)*sizeof(int));

    int next_empty_list = ncell;
    for(i = 0; i < ncell; i++) {
        w->tree_array[i] = next_empty_list;
        w->tree_array[next_empty_list] = list_int_size(tri_list[i]);
        int j;
        for(j = 0; j < w->tree_array[next_empty_list]; j++) {
            w->tree_array[next_empty_list+j+1] = list_int_get(tri_list[i], j);
        }
        next_empty_list += w->tree_array[next_empty_list] + 1;
        list_int_free(&tri_list[i]);
    }
    free(tri_list);
}

/**
 * @brief Construct wall octree recursively
 *
 * Constructs the octree array by iterating through all wall triangles and
 * placing them into an octree structure
 *
 * Note that this function allocates extra space at the end of the tree_array
 * to store wall element flags.
 *
 * @param w pointer to wall data
 * @param offload_array the offload array
 * @param tree_array pointer to array storing what octree cells contain
 *        which triangles
 */
void wall_3d_init_octree(wall_3d_offload_data* w, real* offload_array,
                         int** tree_array) {
    if(w->n > 1000000){
        print_out(VERBOSE_NORMAL,
                  "Starting to initialize 3D-wall octree with %d triangles.\n",
                  w->n);
    }

    /* Construct the octree and store triangles there */
    octree_node* tree;
    octree_create(&tree, w->xmin, w->xmax, w->ymin, w->ymax, w->zmin, w->zmax,
                  w->depth);
    int i;
    for(i = 0; i < w->n; i++) {
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
        octree_add(tree, t1, t2, t3, i);
        if (i%1000000==0 && i > 0){
            print_out(VERBOSE_NORMAL, "  Adding triangle %10d/%d.\n",i,w->n);
        }
    }

    /* Create lists for triangles in each grid square and fill the lists
       by querying the octree in each grid point */
    int ncell = w->ngrid*w->ngrid*w->ngrid;
    list_int_node** tri_list =
        (list_int_node**) malloc(ncell * sizeof(list_int_node*));
    int ix, iy, iz;
    for(ix = 0; ix < w->ngrid; ix++) {
        for(iy = 0; iy < w->ngrid; iy++) {
            for(iz = 0; iz < w->ngrid; iz++) {
                real p[3];
                p[0] = w->xmin + ix * w->xgrid + 0.5*w->xgrid;
                p[1] = w->ymin + iy * w->ygrid + 0.5*w->ygrid;
                p[2] = w->zmin + iz * w->zgrid + 0.5*w->zgrid;

                int cell_index = ix*w->ngrid*w->ngrid+iy*w->ngrid+iz;
                tri_list[cell_index] = octree_get(tree, p);
            }
        }
    }

    /* construct an array from the triangle lists */
    int list_size = 0;
    for(i = 0; i < ncell; i++) {
        list_size += list_int_size(tri_list[i]);
    }

    w->int_offload_array_length = 2*ncell + list_size + w->n;
    *tree_array = (int*) malloc((w->int_offload_array_length)*sizeof(int));

    int next_empty_list = ncell;
    for(i = 0; i < ncell; i++) {
        /* First ncell elements store the position where the actual cell data
         * begins in tree_array */
        (*tree_array)[i] = next_empty_list;

        /* The first data point in the actual cell data is the number of
         * triangles in this cell */
        (*tree_array)[next_empty_list] = list_int_size(tri_list[i]);

        /* Store triangle IDs that are located in this cell */
        for(int j = 0; j < (*tree_array)[next_empty_list]; j++) {
            (*tree_array)[next_empty_list+j+1] = list_int_get(tri_list[i], j);
        }
        next_empty_list += (*tree_array)[next_empty_list] + 1;
    }
    free(tri_list);
    octree_free(&tree);
}

/**
 * @brief Check if trajectory from (r1, phi1, z1) to (r2, phi2, z2) intersects
 *        the wall using the octree structure
 *
 * @param r1 start point R coordinate [m]
 * @param phi1 start point phi coordinate [rad]
 * @param z1 start point z coordinate [rad]
 * @param r2 end point R coordinate [m]
 * @param phi2 end point phi coordinate [rad]
 * @param z2 end point z coordinate [rad]
 * @param wdata pointer to data struct on target
 * @param w_coll pointer for storing the parameter in P = P1 + w_coll * (P2-P1),
 *        where P is the point where the collision occurred.
 *
 * @return id, which is the first element id if hit, zero otherwise
 */
int wall_3d_hit_wall(real r1, real phi1, real z1, real r2, real phi2,
            real z2, wall_3d_data* wdata, real* w_coll) {
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

    int ix1 = (int) floor((q1[0] - wdata->xmin)
                          / ((wdata->xmax - wdata->xmin) / (wdata->ngrid)));
    int iy1 = (int) floor((q1[1] - wdata->ymin)
                          / ((wdata->ymax - wdata->ymin) / (wdata->ngrid)));
    int iz1 = (int) floor((q1[2] - wdata->zmin)
                          / ((wdata->zmax - wdata->zmin) / (wdata->ngrid)));

    int ix2 = (int) floor((q2[0] - wdata->xmin)
                          / ((wdata->xmax - wdata->xmin) / (wdata->ngrid)));
    int iy2 = (int) floor((q2[1] - wdata->ymin)
                          / ((wdata->ymax - wdata->ymin) / (wdata->ngrid)));
    int iz2 = (int) floor((q2[2] - wdata->zmin)
                          / ((wdata->zmax - wdata->zmin) / (wdata->ngrid)));

    int hit_tri = 0;
    real smallest_w = 1.1;

    for(int i = 0; i <= abs(ix2-ix1); i++) {
        for(int j = 0; j <= abs(iy2-iy1); j++) {
            for(int k = 0; k <= abs(iz2-iz1); k++) {
                int ix = ix1 + i*((int) copysign(1, ix2-ix1));
                int iy = iy1 + j*((int) copysign(1, iy2-iy1));
                int iz = iz1 + k*((int) copysign(1, iz2-iz1));

                if(ix >= 0 && ix < wdata->ngrid && iy >= 0 && iy < wdata->ngrid
                   && iz >= 0 && iz < wdata->ngrid) {

                    int ilist = wdata->tree_array[ix*wdata->ngrid*wdata->ngrid
                                                  + iy*wdata->ngrid + iz];

                    for(int l = 0; l < wdata->tree_array[ilist]; l++) {
                        int itri = wdata->tree_array[ilist+l+1];
                        real w = wall_3d_tri_collision(
                            q1, q2,
                            &wdata->wall_tris[9*itri],
                            &wdata->wall_tris[9*itri+3],
                            &wdata->wall_tris[9*itri+6]);
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

/**
 * @brief Check if trajectory from (r1, phi1, z1) to (r2, phi2, z2) intersects
 *        the wall against all triangles
 *
 * @param r1 start point R coordinate [m]
 * @param phi1 start point phi coordinate [rad]
 * @param z1 start point z coordinate [rad]
 * @param r2 end point R coordinate [m]
 * @param phi2 end point phi coordinate [rad]
 * @param z2 end point z coordinate [rad]
 * @param wdata pointer to data struct on target
 * @param w_coll pointer for storing the parameter in P = P1 + w_coll * (P2-P1),
 *        where P is the point where the collision occurred.
 *
 * @return id is wall element id if hit, zero otherwise*
 */
int wall_3d_hit_wall_full(real r1, real phi1, real z1, real r2, real phi2,
                          real z2, wall_3d_data* wdata, real* w_coll) {
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

    for(j = 0; j < wdata->n; j++) {
        w = wall_3d_tri_collision(q1, q2, &wdata->wall_tris[9*j],
                &wdata->wall_tris[9*j+3], &wdata->wall_tris[9*j+6]);
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

/**
 * @brief Check if any part of a triangle is inside a box
 *
 * @param t1 xyz coordinates of first triangle vertex [m]
 * @param t2 xyz coordinates of second triangle vertex [m]
 * @param t3 xyz coordinates of third triangle vertex [m]
 * @param bb1 bounding box minimum xyz coordinates [m]
 * @param bb2 bounding box maximum xyz coordinates [m]
 *
 * @return zero if not any part of the triangle is within the box
 */
int wall_3d_tri_in_cube(real t1[3], real t2[3], real t3[3], real bb1[3],
                        real bb2[3]) {
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

    if(   wall_3d_tri_collision(c000, c100, t1, t2, t3) >= 0
       || wall_3d_tri_collision(c000, c010, t1, t2, t3) >= 0
       || wall_3d_tri_collision(c110, c010, t1, t2, t3) >= 0
       || wall_3d_tri_collision(c110, c100, t1, t2, t3) >= 0
       || wall_3d_tri_collision(c000, c001, t1, t2, t3) >= 0
       || wall_3d_tri_collision(c010, c011, t1, t2, t3) >= 0
       || wall_3d_tri_collision(c100, c101, t1, t2, t3) >= 0
       || wall_3d_tri_collision(c110, c111, t1, t2, t3) >= 0
       || wall_3d_tri_collision(c001, c101, t1, t2, t3) >= 0
       || wall_3d_tri_collision(c001, c011, t1, t2, t3) >= 0
       || wall_3d_tri_collision(c111, c011, t1, t2, t3) >= 0
       || wall_3d_tri_collision(c111, c101, t1, t2, t3) >= 0)
        return 1;

    /* check for triangle edges intersecting cube quads */
    if(   wall_3d_quad_collision(t1, t2, c000, c100, c110, c010) == 1
       || wall_3d_quad_collision(t1, t2, c000, c010, c011, c001) == 1
       || wall_3d_quad_collision(t1, t2, c000, c100, c101, c001) == 1
       || wall_3d_quad_collision(t1, t2, c010, c110, c111, c011) == 1
       || wall_3d_quad_collision(t1, t2, c100, c101, c111, c110) == 1
       || wall_3d_quad_collision(t1, t2, c001, c101, c111, c011) == 1)
        return 1;
    if(   wall_3d_quad_collision(t3, t2, c000, c100, c110, c010) == 1
       || wall_3d_quad_collision(t3, t2, c000, c010, c011, c001) == 1
       || wall_3d_quad_collision(t3, t2, c000, c100, c101, c001) == 1
       || wall_3d_quad_collision(t3, t2, c010, c110, c111, c011) == 1
       || wall_3d_quad_collision(t3, t2, c100, c101, c111, c110) == 1
       || wall_3d_quad_collision(t3, t2, c001, c101, c111, c011) == 1)
        return 1;
    if(   wall_3d_quad_collision(t1, t3, c000, c100, c110, c010) == 1
       || wall_3d_quad_collision(t1, t3, c000, c010, c011, c001) == 1
       || wall_3d_quad_collision(t1, t3, c000, c100, c101, c001) == 1
       || wall_3d_quad_collision(t1, t3, c010, c110, c111, c011) == 1
       || wall_3d_quad_collision(t1, t3, c100, c101, c111, c110) == 1
       || wall_3d_quad_collision(t1, t3, c001, c101, c111, c011) == 1)
        return 1;
    return 0;
}

/**
 * @brief Check if a line segment intersects a triangle
 *
 * This routine implements the Möller-Trumbore algorithm.
 *
 * @param q1 line segment start point xyz coordinates [m]
 * @param q2 line segment end point xyz coordinates [m]
 * @param t1 xyz coordinates of first triangle vertex [m]
 * @param t2 xyz coordinates of second triangle vertex [m]
 * @param t3 xyz coordinates of third triangle vertex [m]
 *
 * @return A positive number w which is defined so that vector q1 + w*(q2-q1)
 *         is the intersection point. A negative number is returned if no there
 *         is no intersection
 */
double wall_3d_tri_collision(real q1[3], real q2[3], real t1[3], real t2[3],
            real t3[3]) {
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

/**
 * @brief Check if a line segment intersects a quad (assumed planar)
 *
 * @param q1 line segment start point xyz coordinates [m]
 * @param q2 line segment end point xyz coordinates [m]
 * @param t1 xyz coordinates of first quad vertex [m]
 * @param t2 xyz coordinates of second quad vertex [m]
 * @param t3 xyz coordinates of third quad vertex [m]
 * @param t4 xyz coordinates of fourth quad vertex [m]
 *
 * @return Zero if no intersection, positive number otherwise
 */
int wall_3d_quad_collision(real q1[3], real q2[3], real t1[3], real t2[3],
                           real t3[3], real t4[3]) {
    if(wall_3d_tri_collision(q1, q2, t1, t2, t3) >= 0
       || wall_3d_tri_collision(q1, q2, t1, t3, t4) >= 0)
        return 1;
    else
        return 0;
}
