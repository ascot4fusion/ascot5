/**
 * @file octree.h
 * @brief Header file for octree.c
 */
#ifndef OCTREE_H
#define OCTREE_H

#include "ascot5.h"
#include "list.h"

/**
 * @brief Struct representing single octree node
 *
 * Stores eight child nodes, bounding box of the volume this node encloses, and
 * a linked list containing IDs of triangles belonging to this node.
 */
typedef struct octree_node {
    struct octree_node* n000; /**< [xmin, ymin, zmin] child node        */
    struct octree_node* n100; /**< [xmax, ymin, zmin] child node        */
    struct octree_node* n010; /**< [xmin, ymax, zmin] child node        */
    struct octree_node* n110; /**< [xmax, ymax, zmin] child node        */
    struct octree_node* n001; /**< [xmin, ymin, zmax] child node        */
    struct octree_node* n101; /**< [xmax, ymin, zmax] child node        */
    struct octree_node* n011; /**< [xmin, ymax, zmax] child node        */
    struct octree_node* n111; /**< [xmax, ymax, zmax] child node        */
    real bb1[3];              /**< Bounding box xyz minimum limit       */
    real bb2[3];              /**< Bounding box xyz maximum limit       */
    list_int_node* list;      /**< Linked list for storing triangle IDs */
} octree_node;

#pragma omp declare target
void octree_create(octree_node** node, real x_min, real x_max,
                        real y_min, real y_max, real z_min, real z_max,
                        int depth);
void octree_free(octree_node** node);
void octree_add(octree_node* node, real t1[3], real t2[3], real t3[3], int id);
list_int_node* octree_get(octree_node* node, real p[3]);
#pragma omp end declare target

#endif
