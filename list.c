/**
 * @file list.c
 * @brief Simple linked list
 *
 * Linked list where each node stores an int data and a pointer to the next node
 * in chain.
 */
#include <stdlib.h>
#include "list.h"

/**
 * @brief Create an empty list
 *
 * @param list pointer to the created list
 */
void list_int_create(list_int_node** list) {
    (*list) = (list_int_node*) malloc(sizeof(list_int_node));
    (*list)->data = 0;
    (*list)->next = NULL;

}

/**
 * @brief Deallocate this list and all lists it is linked to
 *
 * @param list pointer to the list to be freed
 */
void list_int_free(list_int_node** list) {
    list_int_node* node = (*list);
    int n = list_int_size(*list);
    list_int_node* next_node = (*list)->next;
    int i;
    for(i = 0; i < n; i++) {
        free(node);
        node = next_node;
        next_node = next_node->next;
    }
    free(node);
    (*list) = NULL;
}

/**
 * @brief Add new node to the end of the chain
 *
 * @param list list node to which new node is linked
 * @param data int value to be stored in the new node
 */
void list_int_add(list_int_node* list, int data) {
    list_int_node* node = list;
    while(node->next != NULL) {
        node = node->next;
    }

    list_int_node* new_node = (list_int_node*) malloc(sizeof(list_int_node));
    new_node->next = NULL;

    node->next = new_node;
    node->data = data;
}

/**
 * @brief Retrieve the data stored in a list node
 *
 * @param list list node
 * @param index node index where data is retrieved, zero refers to current node
 *
 * @return the stored data
 */
int list_int_get(list_int_node* list, int index) {
    int i = 0;
    list_int_node* node = list;
    while(i < index) {
        node = node->next;
        i++;
    }
    return node->data;
}

/**
 * @brief Get list size
 *
 * @param list list node
 *
 * @return number of nodes this node is linked to
 */
int list_int_size(list_int_node* list) {
    int i = 0;
    list_int_node* node = list;
    while(node->next != NULL) {
        node = node->next;
        i++;
    }
    return i;
}
