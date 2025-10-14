/**
 * @file wall_data.h
 * Wall data types.
 */
#ifndef WALL_DATA_H
#define WALL_DATA_H

/**
 * Wall data types.
 */
typedef enum Wall_type
{
    WALL_CONTOUR2D,    /**< Corresponds to WallContour2D.                     */
    WALL_TRIANGULAR3D, /**< Corresponds to WallTriangular3D.                  */
} Wall_type;

/**
 * Parameters for the 2D wall.
 *
 * Note: The start and end point of wall polygon does not have to concide.
 */
typedef struct
{
    size_t n;  /**< Number of points in the wall polygon.                     */
    real *r;   /**< Wall polygon point R coordinates [m].                     */
    real *z;   /**< Wall polygon point z coordinates [m].                     */
    int *flag; /**< Flags used to label wall elements.                        */
} WallContour2D;

/**
 * 3D wall data parameters.
 */
typedef struct
{
    size_t n;            /**< Number of wall triangles.                       */
    size_t depth;        /**< Depth of the octree.                            */
    size_t ngrid;        /**< Number of cells computational volume is divided
                              to in each direction. ngrid = 2^(depth-1).      */
    size_t n_tree_array; /**< Number of elements in tree_array.               */
    int *flag;           /**< Flags used to label wall elements.              */
    real xmin;           /**< Minimum extend on x-direction [m].              */
    real xmax;           /**< Maximum extend on x-direction [m].              */
    real dx;             /**< Octree cell width in x-direction [m].           */
    real ymin;           /**< Minimum extend on y-direction [m].              */
    real ymax;           /**< Maximum extend on y-direction [m].              */
    real dy;             /**< Octree cell width in y-direction [m].           */
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
    size_t *tree_array;
} WallTriangular3D;

/**
 * Wall model simulation data.
 *
 * The intended usage is that only single type of data is used at a time. This
 * is declared using the `type` field.
 */
typedef struct
{
    WallContour2D *contour2d;       /**< 2D contour.                          */
    WallTriangular3D *triangular3d; /**< 3D triangular mesh.                  */
    Wall_type type;                 /**< Current wall type.                   */
} Wall;

#endif
