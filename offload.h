/**
 * @file offload.h
 * @brief Header file for offload.h
 */
#ifndef OFFLOAD_H
#define OFFLOAD_H

typedef struct {
    size_t offload_array_length;
    size_t unpack_pos;
} offload_package;

void offload_init_offload(offload_package* o, real** offload_array);
void offload_free_offload(offload_package* o, real** offload_array);

void offload_pack(offload_package* o, real** offload_array, real* pack_array,
                  size_t pack_length);

#pragma omp declare target
real* offload_unpack(offload_package* o, real* offload_array,
                     size_t pack_length);
#pragma omp end declare target

#endif
