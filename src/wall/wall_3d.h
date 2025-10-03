/**
 * @file wall_3d.h
 * Header file for wall_3d.c
 */
#ifndef WALL_3D_H
#define WALL_3D_H
#include "ascot5.h"
#include "wall.h"
#include "offload.h"


/** Small value to check if x = 0 (i.e. abs(x) < WALL_EPSILON) */
#define WALL_EPSILON 1e-9


/**
 * @brief Initialize 3D wall data and check inputs

 * The default octree depth is defined by macro WALL_OCTREE_DEPTH in wall_3d.h.
 *
 * @param data pointer to the data struct
 * @param nelements number of elements in the wall polygon
 * @param vertices each triangle's vertices' (x,y,z) coordinates in an array as
 *        [..., x1_i, y1_i, z1_i, x2_i, y2_i, z2_i, x3_i, y3_i, z3_i, ...]
 * @param flag integer label for grouping wall elements together.
 *
 * @return zero if initialization succeeded
 */
int WallTriangular3D_init(
    WallTriangular3D* wall, int n, real* vertices, int* flag
);


/**
 * @brief Free allocated resources
 *
 * @param data pointer to data struct
 */
void WallTriangular3D_free(WallTriangular3D* wall);


/**
 * @brief Offload data to the accelerator.
 *
 * @param data pointer to the data struct
 */
void WallTriangular3D_offload(WallTriangular3D* wall);


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
GPU_DECLARE_TARGET_SIMD_UNIFORM(wall)
int WallTriangular3D_hit_wall(
    real r1, real phi1, real z1, real r2, real phi2, real z2, real* w_coll,
    WallTriangular3D* wall
);
DECLARE_TARGET_END


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
GPU_DECLARE_TARGET_SIMD
double WallTriangular3D_tri_collision(
    real q1[3], real q2[3], real t1[3], real t2[3], real t3[3]
);
DECLARE_TARGET_END


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
DECLARE_TARGET
int WallTriangular3D_tri_in_cube(
    real t1[3], real t2[3], real t3[3], real bb1[3], real bb2[3]
);
DECLARE_TARGET_END


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
DECLARE_TARGET
int WallTriangular3D_quad_collision(
    real q1[3], real q2[3], real t1[3], real t2[3], real t3[3], real t4[3]
);
DECLARE_TARGET_END


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
void WallTriangular3D_init_octree(WallTriangular3D* wall);


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
void WallTriangular3D_init_tree(WallTriangular3D* wall, real* offload_array);


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
int WallTriangular3D_hit_wall_full(
    real r1, real phi1, real z1, real r2, real phi2, real z2, real* w_coll,
    WallTriangular3D* wall
);

#endif
