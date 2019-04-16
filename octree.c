/**
 * @file octree.c
 * @brief Simple octree for storing triangles
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
 *
 * This module is used for an efficient collision check in 3D wall model.
 * @see wall_3d.c
 */
#include <stdlib.h>
#include "ascot5.h"
#include "octree.h"
#include "list.h"
#include "wall/wall_3d.h"

/**
 * @brief Create octree of given depth
 *
 * This function creates recursively a complete octree hierarchy with given
 * number of levels. Each node have bounding box asigned and there is a small
 * overlap between the boxes to avoid numerical artifact where a point is
 * exactly between two boxes but belongs to neither.
 *
 * The linked lists on leaf nodes are initialized.
 *
 * @param node pointer to parent node from which the octree sprawls
 * @param x_min minimum x coordinate of the parent node
 * @param x_max maximum x coordinate of the parent node
 * @param y_min minimum y coordinate of the parent node
 * @param y_max maximum y coordinate of the parent node
 * @param z_min minimum z coordinate of the parent node
 * @param z_max maximum z coordinate of the parent node
 * @param depth levels of octree nodes to be created. If depth=1, this node will
 *        be a leaf node. Final volume will be 1/8^(depth-1) th of the initial
 *        volume.
 */
void octree_create(octree_node** node, real x_min, real x_max, real y_min,
                   real y_max, real z_min, real z_max, int depth) {
    octree_node* n = (octree_node*) malloc(sizeof(octree_node));
    (*node) = n;
    real epsilon = 1e-6;
    n->bb1[0] = x_min - epsilon;
    n->bb2[0] = x_max + epsilon;
    n->bb1[1] = y_min - epsilon;
    n->bb2[1] = y_max + epsilon;
    n->bb1[2] = z_min - epsilon;
    n->bb2[2] = z_max + epsilon;
    real x = (x_max + x_min) / 2;
    real y = (y_max + y_min) / 2;
    real z = (z_max + z_min) / 2;
    if(depth > 1) {
        octree_create(&(n->n000), x_min, x, y_min, y, z_min, z, depth-1);
        octree_create(&(n->n100), x, x_max, y_min, y, z_min, z, depth-1);
        octree_create(&(n->n010), x_min, x, y, y_max, z_min, z, depth-1);
        octree_create(&(n->n110), x, x_max, y, y_max, z_min, z, depth-1);
        octree_create(&(n->n001), x_min, x, y_min, y, z, z_max, depth-1);
        octree_create(&(n->n101), x, x_max, y_min, y, z, z_max, depth-1);
        octree_create(&(n->n011), x_min, x, y, y_max, z, z_max, depth-1);
        octree_create(&(n->n111), x, x_max, y, y_max, z, z_max, depth-1);
        n->list = NULL;
    }
    else {
        n->n000 = NULL;
        n->n100 = NULL;
        n->n010 = NULL;
        n->n110 = NULL;
        n->n001 = NULL;
        n->n101 = NULL;
        n->n011 = NULL;
        n->n111 = NULL;
        list_int_create(&(n->list));
    }
}

/**
 * @brief Free octree node and all its child nodes
 *
 * Deallocates node and its child node recursively. Linked lists are also freed.
 *
 * @param node pointer to octree node
 */
void octree_free(octree_node** node) {
    if((*node)->n000 != NULL) {
        octree_free( &((*node)->n000) );
        octree_free( &((*node)->n100) );
        octree_free( &((*node)->n010) );
        octree_free( &((*node)->n110) );
        octree_free( &((*node)->n001) );
        octree_free( &((*node)->n101) );
        octree_free( &((*node)->n011) );
        octree_free( &((*node)->n111) );
    }
    else {
        list_int_free( &((*node)->list) );
    }
    free(*node);
    *node = NULL;
}

/**
 * @brief Add triangle to the node(s) it belongs to.
 *
 * This function uses recursion to travel the octree and find all leaf nodes
 * the given triangle belongs to.In other words, at each step this function
 * determines the child node(s) the triangle belongs to and calls this function
 * for that/those node(s).
 *
 * @param node octree node to which or to which child nodes triangle is added to
 * @param t1 triangle first vertex xyz coordinates
 * @param t2 triangle second vertex xyz coordinates
 * @param t3 triangle third vertex xyz coordinates
 * @param id triangle id which is stored in the node(s) the triangle belongs to
 */
void octree_add(octree_node* node, real t1[3], real t2[3], real t3[3], int id) {
    if(node->n000 == NULL) {
        /* leaf node, add the triangle to the list */
        list_int_add(node->list, id);
    }
    else {
        if(wall_3d_tri_in_cube(t1, t2, t3, node->n000->bb1, node->n000->bb2)>0)
            octree_add(node->n000, t1, t2, t3, id);
        if(wall_3d_tri_in_cube(t1, t2, t3, node->n100->bb1, node->n100->bb2)>0)
            octree_add(node->n100, t1, t2, t3, id);
        if(wall_3d_tri_in_cube(t1, t2, t3, node->n010->bb1, node->n010->bb2)>0)
            octree_add(node->n010, t1, t2, t3, id);
        if(wall_3d_tri_in_cube(t1, t2, t3, node->n110->bb1, node->n110->bb2)>0)
            octree_add(node->n110, t1, t2, t3, id);
        if(wall_3d_tri_in_cube(t1, t2, t3, node->n001->bb1, node->n001->bb2)>0)
            octree_add(node->n001, t1, t2, t3, id);
        if(wall_3d_tri_in_cube(t1, t2, t3, node->n101->bb1, node->n101->bb2)>0)
            octree_add(node->n101, t1, t2, t3, id);
        if(wall_3d_tri_in_cube(t1, t2, t3, node->n011->bb1, node->n011->bb2)>0)
            octree_add(node->n011, t1, t2, t3, id);
        if(wall_3d_tri_in_cube(t1, t2, t3, node->n111->bb1, node->n111->bb2)>0)
            octree_add(node->n111, t1, t2, t3, id);
    }
}

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
list_int_node* octree_get(octree_node* node, real p[3]) {
    if(node->n000 == NULL) {
        return node->list;
    }
    else {
        real x = (node->bb1[0] + node->bb2[0]) / 2;
        real y = (node->bb1[1] + node->bb2[1]) / 2;
        real z = (node->bb1[2] + node->bb2[2]) / 2;

        if(p[0] < x && p[1] < y && p[2] < z)
            return octree_get(node->n000, p);
        if(p[0] > x && p[1] < y && p[2] < z)
            return octree_get(node->n100, p);
        if(p[0] < x && p[1] > y && p[2] < z)
            return octree_get(node->n010, p);
        if(p[0] > x && p[1] > y && p[2] < z)
            return octree_get(node->n110, p);
        if(p[0] < x && p[1] < y && p[2] > z)
            return octree_get(node->n001, p);
        if(p[0] > x && p[1] < y && p[2] > z)
            return octree_get(node->n101, p);
        if(p[0] < x && p[1] > y && p[2] > z)
            return octree_get(node->n011, p);
        if(p[0] > x && p[1] > y && p[2] > z)
            return octree_get(node->n111, p);
    }

    /* We should not ever get here but this
       supresses compiler warnings. */
    return NULL;
}
