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
#include "list.h"
#include "parallel.h"

/** Small value to check if x = 0 (i.e. abs(x) < WALL_EPSILON) [m]. */
#define WALL_EPSILON 1e-9

/**
 * Struct representing single octree node.
 *
 * Stores eight child nodes, bounding box of the volume this node encloses, and
 * a linked list containing IDs of triangles belonging to this node.
 */
typedef struct Octree
{
    struct Octree *n000; /**< The [xmin, ymin, zmin] of child nodes [m].      */
    struct Octree *n100; /**< The [xmax, ymin, zmin] of child nodes [m].      */
    struct Octree *n010; /**< The [xmin, ymax, zmin] of child nodes [m].      */
    struct Octree *n110; /**< The [xmax, ymax, zmin] of child nodes [m].      */
    struct Octree *n001; /**< The [xmin, ymin, zmax] of child nodes [m].      */
    struct Octree *n101; /**< The [xmax, ymin, zmax] of child nodes [m].      */
    struct Octree *n011; /**< The [xmin, ymax, zmax] of child nodes [m].      */
    struct Octree *n111; /**< The [xmax, ymax, zmax] of child nodes [m].      */
    float bb1[3];        /**< Bounding box xyz minimum limit [m].             */
    float bb2[3];        /**< Bounding box xyz maximum limit [m].             */
    ListInt *list;       /**< Linked list for storing triangle IDs.           */
} Octree;

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
 * @param x_min Minimum x coordinate of the parent node [m].
 * @param x_max Maximum x coordinate of the parent node [m].
 * @param y_min Minimum y coordinate of the parent node [m].
 * @param y_max Maximum y coordinate of the parent node [m].
 * @param z_min Minimum z coordinate of the parent node [m].
 * @param z_max Maximum z coordinate of the parent node [m].
 * @param depth Levels of octree nodes to be created. If depth=1, this node will
 *        be a leaf node. Final volume will be 1/8^(depth-1) th of the initial
 *        volume.
 * @return 0 on success and 1 if failed to allocate sufficient memory.
 */
int Octree_create(
    Octree **node, float x_min, float x_max, float y_min, float y_max,
    float z_min, float z_max, size_t depth);

/**
 * Free octree node and all its child nodes.
 *
 * Deallocates node and its child node recursively. Linked lists are also freed.
 *
 * @param node Pointer to octree node.
 */
void Octree_free(Octree **node);

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
 * @param id Triangle ID which is stored in the node(s) the triangle belongs to.
 */
void Octree_add(Octree *node, float t1[3], float t2[3], float t3[3], size_t id);

/**
 * Get that leaf node's linked list the given coordinate belongs to.
 *
 * This function uses recursion to travel through the octree, determining at
 * each step the correct branch to follow next until the leaf node is found.
 * In other words, at each step this function determines the child node the
 * point belongs to and calls this function for that node.
 *
 * The point is assumed to belong to the volume of the node used in the
 * argument.
 *
 * @param node Octree node that is traversed.
 * @param p Point xyz coordinates.
 *
 * @return Linked list of the leaf node given point belongs to.
 */
ListInt *Octree_get(Octree *node, real p[3]);

DECLARE_TARGET
/**
 * Check if any part of a triangle is inside a box.
 *
 * @param t1 First triangle vertex xyz coordinates [m].
 * @param t2 Second triangle vertex xyz coordinates [m].
 * @param t3 Third triangle vertex xyz coordinates [m].
 * @param bb1 Bounding box minimum xyz coordinates [m].
 * @param bb2 Bounding box maximum xyz coordinates [m].
 *
 * @return Zero if not any part of the triangle is within the box, positive
 *         number otherwise
 */
int Octree_tri_in_cube(
    float t1[3], float t2[3], float t3[3], float bb1[3], float bb2[3]);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD
/**
 * Check if a line segment intersects a triangle.
 *
 * This routine implements the Möller-Trumbore algorithm.
 *
 * @param q1 Line segment start point xyz coordinates [m].
 * @param q2 Line segment end point xyz coordinates [m].
 * @param t1 First triangle vertex xyz coordinates [m].
 * @param t2 Second triangle vertex xyz coordinates [m].
 * @param t3 Third triangle vertex xyz coordinates [m].
 *
 * @return A positive number w which is defined so that vector q1 + w*(q2-q1)
 *         is the intersection point. A negative number is returned if there
 *         is no intersection.
 */
float Octree_tri_collision(
    real q1[3], real q2[3], float t1[3], float t2[3], float t3[3]);
DECLARE_TARGET_END

/**
 * Check if a line segment intersects a quad (assumed planar).
 *
 * @param q1 Line segment start point xyz coordinates [m].
 * @param q2 Line segment end point xyz coordinates [m].
 * @param t1 First quad vertex xyz coordinates [m].
 * @param t2 Second quad vertex xyz coordinates [m].
 * @param t3 Third quad vertex xyz coordinates [m].
 * @param t3 Fourth quad vertex xyz coordinates [m].
 *
 * @return Zero if no intersection, positive number otherwise.
 */
#define Octree_quad_collision(q1, q2, t1, t2, t3, t4)                          \
    (Octree_tri_collision((q1), (q2), (t1), (t2), (t3)) >= 0 ||                \
     Octree_tri_collision((q1), (q2), (t1), (t3), (t4)) >= 0)

#endif
