/**
 * @file offload.h
 * @brief Header file for offload.h
 */
#ifndef OFFLOAD_H
#define OFFLOAD_H

typedef struct {
    real* offload_array;
    size_t offload_array_length;
    size_t unpack_pos;
} offload_package;

void offload_init_offload(offload_package* o);
void offload_free_offload(offload_package* o);

void offload_pack(offload_package* o, real* pack_array, size_t pack_length);

#pragma omp declare target
void offload_init(offload_package* o, size_t offload_array_length,
                  real* offload_array);
real* offload_unpack(offload_package* o, size_t pack_length);
#pragma omp end declare target

#endif
