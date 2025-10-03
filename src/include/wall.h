/**
 * Interface through which wall data is initialized and accessed.
 */
#ifndef WALL_H
#define WALL_H

#include "ascot5.h"
#include "offload.h"

/**
 * Wall data types.
 */
typedef enum wall_type
{
    wall_contour2d,    /**< Corresponds to WallContour2D.                     */
    wall_triangular3d, /**< 3D wall model consisting of triangles */
} wall_type;

/**
 * 2D wall data parameters.
 *
 * Note: The start and end point of wall polygon does not have to concide.
 */
typedef struct
{
    int n;     /**< Number of points in the wall polygon      */
    real *r;   /**< R coordinates for the wall polygon points */
    real *z;   /**< z coordinates for the wall polygon points */
    int *flag; /**< Array of wall element flags               */
} WallContour2D;

/**
 * 3D wall data parameters.
 */
typedef struct
{
    int n;               /**< Number of wall triangles.                       */
    int depth;           /**< Depth of the octree.                            */
    int ngrid;           /**< Number of cells computational volume is divided
                              to in each direction. ngrid = 2^(depth-1).      */
    int tree_array_size; /* Number of elements in tree_array.                 */
    int *flag;           /**< Array of wall triangle flags.                   */
    real xmin;           /**< Minimum extend on x-direction [m].              */
    real xmax;           /**< Maximum extend on x-direction [m].              */
    real dx;             /**< Octree cell width in x-direction [m].           */
    real ymin;           /**< Minimum extend on y-direction [m].              */
    real ymax;           /**< Maximum extend on y-direction [m].              */
    real dy;             /**<  Octree cell width in y-direction [m].          */
    real zmin;           /**< Minimum extend on z-direction [m].              */
    real zmax;           /**< Maximum extend on z-direction [m].              */
    real dz;             /**< Octree cell width in z-direction [m].           */
    real *vertices;      /**< Array of wall triangle coordinates.             */

    /**
     * Array storing information what triangles given octree cell stores.
     *
     * First ncell elements store the array position where data for a given cell
     * begins, where the cell index is icell = ix * ngrid**2 + iy * ngrid + iz.
     * The first element in the cell data, i.e. tree_array[tree_array[icell]],
     * is the number of triangles in this cell, ntriangle, and the next
     * ntriangle elements are the triangle indices.
     */
    int *tree_array;
} WallTriangular3D;

/**
 * Wall model simulation data.
 *
 * The intended usage is that only single type of data is used at a time. This
 * is declared using the `type` field.
 */
typedef struct
{
    WallContour2D *contour2d; /**< 2D model or NULL if not active         */
    WallTriangular3D *triangular3d; /**< 3D model or NULL if not active */
    wall_type type; /**< Wall model type wrapped by this struct */
} wall_data;

/**
 * Free allocated resources.
 *
 * @param wall Pointer to the data struct.
 */
void wall_free(wall_data *wall);

/**
 *  Offload data to the accelerator.
 *
 * @param wall pointer to the data struct.
 */
void wall_offload(wall_data *wall);

/**
 * Check if a given directed line segment intersects the wall.
 *
 * This function is intended to be used to check whether a marker collides with
 * the wall. If there is a collision, this function returns an identification
 * number specific to that wall tile. If the marker hits multiple wall elements,
 * only the first one is returned.
 *
 * This is a SIMD function.
 *
 * @param r1 Start point R coordinate [m].
 * @param phi1 Start point phi coordinate [rad].
 * @param z1 Start point z coordinate [rad].
 * @param r2 End point R coordinate [m].
 * @param phi2 End point phi coordinate [rad].
 * @param z2 End point z coordinate [rad].
 * @param wall Pointer to data struct on target.
 * @param w_coll Pointer for storing the parameter in P = P1 + w_coll * (P2-P1),
 *        where P is the point where the collision occurred.
 *
 * @return Wall element id if hit, zero otherwise.
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(wall)
int wall_hit_wall(
    real r1, real phi1, real z1, real r2, real phi2, real z2, wall_data *wall,
    real *w_coll);
DECLARE_TARGET_END

/**
 * Return the number of wall elements.
 *
 * @param wall pointer to wall data struct on target.
 *
 * @return Number of wall elements or zero on failure.
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(wall)
int wall_get_n_elements(wall_data *wall);
DECLARE_TARGET_END

/**
 * Return the flag of a wall element.
 *
 * @param wall Pointer to wall data struct on target.
 * @param idx Wall element index.
 *
 * @return Flag of the wall element.
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(wall)
int wall_get_flag(wall_data *wall, int idx);
DECLARE_TARGET_END

#endif
