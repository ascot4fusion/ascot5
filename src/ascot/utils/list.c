/**
 * Implements list.h.
 */
#include "list.h"
#include <stdlib.h>

void list_int_create(list_int_node **list)
{
    (*list) = (list_int_node *)malloc(sizeof(list_int_node));
    (*list)->data = 0;
    (*list)->next = NULL;
}

void list_int_free(list_int_node **list)
{
    list_int_node *node = (*list);
    int n = list_int_size(*list);
    list_int_node *next_node = (*list)->next;
    int i;
    for (i = 0; i < n; i++)
    {
        free(node);
        node = next_node;
        next_node = next_node->next;
    }
    free(node);
    (*list) = NULL;
}

void list_int_add(list_int_node *list, int data)
{
    list_int_node *node = list;
    while (node->next != NULL)
    {
        node = node->next;
    }

    list_int_node *new_node = (list_int_node *)malloc(sizeof(list_int_node));
    new_node->next = NULL;

    node->next = new_node;
    node->data = data;
}

int list_int_get(list_int_node *list, int index)
{
    int i = 0;
    list_int_node *node = list;
    while (i < index)
    {
        node = node->next;
        i++;
    }
    return node->data;
}

int list_int_size(list_int_node *list)
{
    int i = 0;
    list_int_node *node = list;
    while (node->next != NULL)
    {
        node = node->next;
        i++;
    }
    return i;
}
