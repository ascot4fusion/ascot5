/**
 * Simple octree for storing triangles.
 *
 * In octree each node branches to eight childnodes. Octree is used to
 * repeatedly divide volume to eight boxes of equal size and volume until the
 * volume in the leaf nodes is of desired size. Boxes are aligned to cartesian
 * basis vectors.
 *
 * This particular octree is intended to be used to group 3D triangles so that
 * a triangle is assigned to an octree node whose volume the triangle belongs
 * to. A triangle belongs to a volume if even one of its vertices is in that
 * volume. One triangle can therefore belong to several octree leaf nodes. Only
 * leaf nodes store triangles.
 */
#ifndef OCTREE_H
#define OCTREE_H
#include "defines.h"
#include "parallel.h"
#include "list.h"

/** Small value to check if x = 0 (i.e. abs(x) < WALL_EPSILON) */
#define WALL_EPSILON 1e-9


/**
 * Struct representing single octree node.
 *
 * Stores eight child nodes, bounding box of the volume this node encloses, and
 * a linked list containing IDs of triangles belonging to this node.
 */
typedef struct octree_node {
    struct octree_node* n000; /**< The [xmin, ymin, zmin] child node.         */
    struct octree_node* n100; /**< The [xmax, ymin, zmin] child node.         */
    struct octree_node* n010; /**< The [xmin, ymax, zmin] child node.         */
    struct octree_node* n110; /**< The [xmax, ymax, zmin] child node.         */
    struct octree_node* n001; /**< The [xmin, ymin, zmax] child node.         */
    struct octree_node* n101; /**< The [xmax, ymin, zmax] child node.         */
    struct octree_node* n011; /**< The [xmin, ymax, zmax] child node.         */
    struct octree_node* n111; /**< The [xmax, ymax, zmax] child node.         */
    real bb1[3];              /**< Bounding box xyz minimum limit.            */
    real bb2[3];              /**< Bounding box xyz maximum limit.            */
    list_int_node* list;      /**< Linked list for storing triangle IDs.      */
} octree_node;


/**
 * Create octree of given depth.
 *
 * This function creates recursively a complete octree hierarchy with given
 * number of levels. Each node have bounding box asigned and there is a small
 * overlap between the boxes to avoid numerical artifact where a point is
 * exactly between two boxes but belongs to neither.
 *
 * The linked lists on leaf nodes are initialized.
 *
 * @param node Pointer to parent node from which the octree sprawls.
 * @param x_min Minimum x coordinate of the parent node.
 * @param x_max Maximum x coordinate of the parent node.
 * @param y_min Minimum y coordinate of the parent node.
 * @param y_max Maximum y coordinate of the parent node.
 * @param z_min Minimum z coordinate of the parent node.
 * @param z_max Maximum z coordinate of the parent node.
 * @param depth Levels of octree nodes to be created. If depth=1, this node will
 *        be a leaf node. Final volume will be 1/8^(depth-1) th of the initial
 *        volume.
 */
void octree_create(
    octree_node** node, real x_min, real x_max, real y_min, real y_max,
    real z_min, real z_max, int depth
);


/**
 * Free octree node and all its child nodes.
 *
 * Deallocates node and its child node recursively. Linked lists are also freed.
 *
 * @param node Pointer to octree node.
 */
void octree_free(octree_node** node);


/**
 * Add triangle to the node(s) it belongs to.
 *
 * This function uses recursion to travel the octree and find all leaf nodes
 * the given triangle belongs to.In other words, at each step this function
 * determines the child node(s) the triangle belongs to and calls this function
 * for that/those node(s).
 *
 * @param node Node to which or to which child nodes triangle is added to.
 * @param t1 Triangle first vertex xyz coordinates.
 * @param t2 Triangle second vertex xyz coordinates.
 * @param t3 Triangle third vertex xyz coordinates.
 * @param id Triangle id which is stored in the node(s) the triangle belongs to.
 */
void octree_add(octree_node* node, real t1[3], real t2[3], real t3[3], int id);

/**
 * @brief Get that leaf node's linked list the given coordinate belongs to
 *
 * This function uses recursion to travel through the octree, determining at
 * each step the correct branch to follow next until the leaf node is found.
 * In other words, at each step this function determines the child node the
 * point belongs to and calls this function for that node.
 *
 * The point is assumed to belong to the volume of the node used in the
 * argument.
 *
 * @param node octree node that is traversed
 * @param p xyz coordinates of the point
 *
 * @return linked list of the leaf node given point belongs to
 */
list_int_node* octree_get(octree_node* node, real p[3]);


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
int octree_tri_in_cube(
    real t1[3], real t2[3], real t3[3], real bb1[3], real bb2[3]
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
double octree_tri_collision(
    real q1[3], real q2[3], real t1[3], real t2[3], real t3[3]
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
int octree_quad_collision(
    real q1[3], real q2[3], real t1[3], real t2[3], real t3[3], real t4[3]
);
DECLARE_TARGET_END

#endif
