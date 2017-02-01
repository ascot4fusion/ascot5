/**
 * @file octree.h
 * @brief Header file for octree.c
 */
#ifndef OCTREE_H
#define OCTREE_H

#include "ascot5.h"
#include "list.h"

typedef struct octree_node {
    struct octree_node* n000;
    struct octree_node* n100;
    struct octree_node* n010;
    struct octree_node* n110;
    struct octree_node* n001;
    struct octree_node* n101;
    struct octree_node* n011;
    struct octree_node* n111;
    real bb1[3], bb2[3];
    list_int_node* list;
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
