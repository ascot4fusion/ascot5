/**
 * Linked list.
 */
#ifndef LIST_H
#define LIST_H
#include "parallel.h"

/**
 * Linked list node that stores data of type `int`.
 */
typedef struct list_int_node
{
    int data;                   /**< Data that is stored in this node.        */
    struct list_int_node *next; /**< Next node in chain or NULL if last.      */
} list_int_node;

/**
 * Create an empty list.
 *
 * @param list Pointer to the created list.
 */
DECLARE_TARGET
void list_int_create(list_int_node **list);
DECLARE_TARGET_END

/**
 * Deallocate this list and all lists it is linked to.
 *
 * @param list Pointer to the list to be freed.
 */
DECLARE_TARGET
void list_int_free(list_int_node **list);
DECLARE_TARGET_END

/**
 * Add new node to the end of the chain.
 *
 * @param list List node to which new node is linked.
 * @param data Value to be stored in the new node.
 */
DECLARE_TARGET
void list_int_add(list_int_node *list, int data);
DECLARE_TARGET_END

/**
 * Retrieve the data stored in a list node.
 *
 * @param list First node in the list.
 * @param index Node index where data is retrieved (use zero for this node).
 *
 * @return The stored data.
 */
DECLARE_TARGET
int list_int_get(list_int_node *list, int index);
DECLARE_TARGET_END

/**
 * Get number of nodes that come after this node.
 *
 * @param list List node.
 *
 * @return Number of nodes this node is followed by.
 */
DECLARE_TARGET
int list_int_size(list_int_node *list);
DECLARE_TARGET_END

#endif
