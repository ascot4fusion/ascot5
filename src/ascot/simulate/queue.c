#include "simulate.h"
#include "data/bfield.h"
#include "data/marker.h"
#include <stdlib.h>

/**
 * Iterate through marker array and return index of the finished marker.
 *
 * This function iterates through the markers that are currently being simulated
 * and stops when it finds the first marker that is finished.
 */
size_t MarkerQueue_cycle(
    size_t *next_in_queue, MarkerQueue *q, size_t nmrk, size_t start,
    size_t ids[nmrk], int running[nmrk]) {
    size_t idx;
    for (idx = start; idx < nmrk; idx++)
    {
        int marker_finished = ids[idx] > 0 && !running[idx];
        int vector_initial_fill = ids[idx] == 0 && q->next < q->n;
        if (marker_finished)
        {
            #pragma omp critical
            q->finished++;
            break;
        }
        if (vector_initial_fill)
            break;
    }

    *next_in_queue = q->n;
    int vector_has_empty_slot = idx < nmrk;
    if(vector_has_empty_slot) {
        size_t next;
        #pragma omp critical
        next = q->next++;
        int empty_queue = next >= q->n;
        if(!empty_queue)
            *next_in_queue = next;
    }

    return idx;
}
