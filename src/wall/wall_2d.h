/**
 * @file wall_2d.h
 * Header file for wall_2d.c
 */
#ifndef WALL_2D_H
#define WALL_2D_H
#include "ascot5.h"
#include "wall.h"
#include "offload.h"

/**
 * @brief Load 2D wall data and prepare parameters
 *
 * @param data pointer to the data struct
 * @param nelements number of elements in the wall polygon
 * @param r R coordinates for the wall polygon points
 * @param z z coordinates for the wall polygon points
 * @param flag integer label for grouping wall elements together.
 *
 * @return zero to indicate success
 */
int WallContour2D_init(
    WallContour2D* wall, int n, real r[n], real z[n], int flag[n]
);


/**
 * @brief Free allocated resources
 *
 * @param data pointer to the data struct
 */
void WallContour2D_free(WallContour2D* wall);


/**
 * @brief Offload data to the accelerator.
 *
 * @param data pointer to the data struct
 */
void WallContour2D_offload(WallContour2D* wall);


/**
 * @brief Check if coordinates are within 2D polygon wall
 *
 * This function checks if the given coordinates are within the walls defined
 * by a 2D polygon using a modified axis crossing method Origin is moved
 * to the coordinates and the number of wall segments crossing the positive
 * R-axis are calculated. If this is odd, the point is inside the polygon.
 *
 * @param r R coordinate [m]
 * @param z z coordinate [m]
 * @param w 2D wall data structure
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(wall)
int WallContour2D_inside(real r, real z, WallContour2D* wall);
DECLARE_TARGET_END


/**
 * @brief Check if trajectory from (r1, phi1, z1) to (r2, phi2, z2) intersects
 *        the wall
 *
 * @param r1 start point R coordinate [m]
 * @param phi1 start point phi coordinate [rad]
 * @param z1 start point z coordinate [m]
 * @param r2 end point R coordinate [m]
 * @param phi2 end point phi coordinate [rad]
 * @param z2 end point z coordinate [m]
 * @param w pointer to data struct on target
 * @param w_coll pointer for storing the parameter in P = P1 + w_coll * (P2-P1),
 *        where P is the point where the collision occurred.
 *
 * @return wall element ID if hit, zero otherwise
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(wall)
int WallContour2D_hit_wall(
    real r1, real z1, real r2, real z2, real* w_coll, WallContour2D* wall
);


/**
 * @brief Find intersection between the wall element and line segment
 *
 * If there are multiple intersections, the one that is closest to P1
 * is returned.
 *
 * @param r1 R1 coordinate of the line segment [P1,P2] [m]
 * @param z1 z1 coordinate of the line segment [P1,P2] [m]
 * @param r2 R2 coordinate of the line segment [P1,P2] [m]
 * @param z2 z2 coordinate of the line segment [P1,P2] [m]
 * @param w pointer to the wall data
 * @param w_coll pointer for storing the parameter in P = P1 + w_coll * (P2-P1),
 *        where P is the point where the collision occurred.
 *
 * @return int wall element id if hit, zero otherwise
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(wall)
int WallContour2D_find_intersection(
    real r1, real z1, real r2, real z2, real* w_coll, WallContour2D* wall
);
DECLARE_TARGET_END

#endif
