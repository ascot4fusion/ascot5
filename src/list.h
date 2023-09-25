/**
 * @file list.h
 * @brief Header file for list.c
 */
#ifndef LIST_H
#define LIST_H
#include "offload_acc_omp.h"

/**
 * @brief Linked list node that stores int data
 */
typedef struct list_int_node {
    struct list_int_node* next; /**< Next node in chain               */
    int data;                   /**< Data that is stored in this node */
} list_int_node;

#pragma omp declare target
DECLARE_TARGET
void list_int_create(list_int_node** list);
DECLARE_TARGET_END
DECLARE_TARGET
void list_int_free(list_int_node** list);
DECLARE_TARGET_END
DECLARE_TARGET
void list_int_add(list_int_node* list, int data);
DECLARE_TARGET_END
DECLARE_TARGET
int list_int_get(list_int_node* list, int index);
DECLARE_TARGET_END
DECLARE_TARGET
int list_int_size(list_int_node* list);
DECLARE_TARGET_END
#pragma omp end declare target

#endif
