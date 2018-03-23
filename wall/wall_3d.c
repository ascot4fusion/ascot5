/**
 * @file wall_3d.c
 * @brief 3D wall collision checks
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../ascot5.h"
#include "wall_3d.h"
#include "../math.h"
#include "../list.h"
#include "../octree.h"

/**
 * @brief Load 3D wall data and prepare parameters
 *
 * Reads a 3D wall from input.wall_3d (not the same format as ASCOT4!), stores
 * parameters in offload struct and allocates and fills the offload array.
 *
 * The input file contains an integer giving the number of wall triangles
 * followed by a list of coordinates for the corners of each triangle
 * (x1,y1,z1,x2,y2,z2,x3,y3,z3).
 *
 * @todo Error checking
 * @todo Implement ASCOT4 input.wall_3d into ascot4_interface
 */
void wall_3d_init_offload(wall_3d_offload_data* offload_data,
                          real** offload_array) {
    // Dummy function
}

/**
 * @brief Free offload array and reset parameters
 */
void wall_3d_free_offload(wall_3d_offload_data* offload_data,
                          real** offload_array) {
    free(*offload_array);
    *offload_array = NULL;
}

/**
 * @brief Initialize 3D wall data struct on target
 */
void wall_3d_init(wall_3d_data* w, wall_3d_offload_data* offload_data,
                  real* offload_array) {
    w->n = offload_data->n;
    w->xmin = offload_data->xmin;
    w->xmax = offload_data->xmax;
    w->xgrid = offload_data->xgrid;
    w->ymin = offload_data->ymin;
    w->ymax = offload_data->ymax;
    w->ygrid = offload_data->ygrid;
    w->zmin = offload_data->zmin;
    w->zmax = offload_data->zmax;
    w->zgrid = offload_data->zgrid;
    w->depth = offload_data->depth;
    w->ngrid = offload_data->ngrid;
    w->wall_tris = &offload_array[0];
    wall_3d_init_octree(w, offload_array);

    /* old slow method */
    /* wall_3d_init_tree(w, offload_array); */
}

/**
 * @brief Construct wall octree iteratively
 *
 * Constructs the octree array by iterating through all wall triangles and
 * octree grid to identify triangles belonging to each grid cell.
 *
 * Slow, only for testing purposes.
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
 * placing them into an octree structure.
 */
void wall_3d_init_octree(wall_3d_data* w, real* offload_array) {
    /* construct the octree and store triangles there */
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
    }

    /* create lists for triangles in each grid square and fill the lists
       by querying the octree in each grid point */
    int ncell = w->ngrid*w->ngrid*w->ngrid;
    list_int_node** tri_list = (list_int_node**) malloc(ncell
                                                    *sizeof(list_int_node*));
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
    }
    free(tri_list);
} 

/**
 * @brief Check if trajectory from (r1, phi1, z1) to (r2, phi2, z2) intersects
 *        the wall using the octree structure
 */
int wall_3d_hit_wall(real r1, real phi1, real z1, real r2, real phi2,
                     real z2, wall_3d_data* wdata) {
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
                    real w = wall_3d_tri_collision(q1, q2,
                                                   &wdata->wall_tris[9*itri],
                                                   &wdata->wall_tris[9*itri+3],
                                                   &wdata->wall_tris[9*itri+6]);
                    if(w >= 0) {
                        hit_tri = itri+1;
                    }
                }
                }
            }
        }
    }

    return hit_tri;
}

/**
 * @brief Check if trajectory from (r1, phi1, z1) to (r2, phi2, z2) intersects
 *        the wall against all triangles
 */
int wall_3d_hit_wall_full(real r1, real phi1, real z1, real r2, real phi2,
                          real z2, wall_3d_data* wdata) {
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

    real w;
    int j;
    for(j = 0; j < wdata->n; j++) {
        w = wall_3d_tri_collision(q1, q2, &wdata->wall_tris[9*j],
                &wdata->wall_tris[9*j+3], &wdata->wall_tris[9*j+6]);
        if(w >= 0) {
            break;
        }
    }

    if(j == wdata->n) {
        return 0;
    }
    else {
        return 1;
    }
}

/**
 * @brief Check if any part of a triangle is inside a box
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
 */
double wall_3d_tri_collision(real q1[3], real q2[3], real t1[3], real t2[3],
                             real t3[3]) {
    real t2t1[3];
    t2t1[0] = t2[0] - t1[0];
    t2t1[1] = t2[1] - t1[1];
    t2t1[2] = t2[2] - t1[2];

    real t3t1[3];
    t3t1[0] = t3[0] - t1[0];
    t3t1[1] = t3[1] - t1[1];
    t3t1[2] = t3[2] - t1[2];

    real n[3];
    math_cross(t2t1, t3t1, n);

    real t1q1[3];
    t1q1[0] = t1[0]-q1[0];
    t1q1[1] = t1[1]-q1[1];
    t1q1[2] = t1[2]-q1[2];

    real q2q1[3];
    q2q1[0] = q2[0]-q1[0];
    q2q1[1] = q2[1]-q1[1];
    q2q1[2] = q2[2]-q1[2];

    real w = math_dot(t1q1, n) / math_dot(q2q1, n);

    real p[3];
    p[0] = q1[0] + w * (q2[0] - q1[0]) - t1[0];
    p[1] = q1[1] + w * (q2[1] - q1[1]) - t1[1];
    p[2] = q1[2] + w * (q2[2] - q1[2]) - t1[2];

    if(w >= 0 && w <= 1) {
        real v = (math_dot(t3t1,p)
             - math_dot(t3t1,t2t1) * math_dot(t2t1,p) / math_dot(t2t1,t2t1))
            / (math_dot(t3t1,t3t1) - math_dot(t2t1,t3t1)*math_dot(t2t1,t3t1)
                                     / math_dot(t2t1,t2t1));
        real u = (math_dot(t2t1, p) - v * math_dot(t2t1,t3t1))
                 / math_dot(t2t1,t2t1);

        if(v >= 0 && u >= 0 && v+u <= 1) {
            return w;
        } else {
            return -1;
        }
    }
    else {
        return -1;
    }
}

/**
 * @brief Check if a line segment intersects a quad (assumed planar)
 */
int wall_3d_quad_collision(real q1[3], real q2[3], real t1[3], real t2[3],
                           real t3[3], real t4[3]) {
    if(wall_3d_tri_collision(q1, q2, t1, t2, t3) >= 0
       || wall_3d_tri_collision(q1, q2, t1, t3, t4) >= 0)
        return 1;
    else
        return 0;
}
