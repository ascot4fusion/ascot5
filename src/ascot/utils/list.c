/**
 * Implements "list.h".
 */
#include "list.h"
#include <stdlib.h>

void list_int_create(ListInt **list)
{
    (*list) = (ListInt *)malloc(sizeof(ListInt));
    (*list)->data = 0;
    (*list)->next = NULL;
}

void list_int_free(ListInt **list)
{
    ListInt *node = (*list);
    int n = list_int_size(*list);
    ListInt *next_node = (*list)->next;
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

void list_int_add(ListInt *list, int data)
{
    ListInt *node = list;
    while (node->next != NULL)
    {
        node = node->next;
    }

    ListInt *new_node = (ListInt *)malloc(sizeof(ListInt));
    new_node->next = NULL;

    node->next = new_node;
    node->data = data;
}

int list_int_get(ListInt *list, int index)
{
    int i = 0;
    ListInt *node = list;
    while (i < index)
    {
        node = node->next;
        i++;
    }
    return node->data;
}

int list_int_size(ListInt *list)
{
    int i = 0;
    ListInt *node = list;
    while (node->next != NULL)
    {
        node = node->next;
        i++;
    }
    return i;
}
