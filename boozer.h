/**
 * @file boozer.h
 * @brief Header file for boozer.c
 */
#ifndef BOOZER_H
#define BOOZER_H

#include "ascot5.h"
#include "error.h"

/**
 * @brief Boozer parameters that will be offloaded to target
 */
typedef struct {
    int offload_array_length; /**< Number of elements in offload_array        */
} boozer_offload_data;

/**
 * @brief Boozer parameters on the target
 */
typedef struct {
} boozer_data;

int boozer_init_offload(boozer_offload_data* offload_data,
                        real** offload_array);

#pragma omp declare target
void boozer_init(boozer_data* boozerdata, boozer_offload_data* offload_data,
                 real* offload_array);

#pragma omp declare simd uniform(boozerdata)
a5err boozer_eval(boozer_data* boozerdata);

#pragma omp end declare target

#endif
