/**
 * @file offload.h
 * @brief Header file for offload.h
 */
#ifndef OFFLOAD_H
#define OFFLOAD_H

/**
 * @brief Struct to keep track of the offload array length and unpack status.
 */
typedef struct {
    size_t offload_array_length; /**< Total length of the common offload
                                      array                                 */
    size_t unpack_pos;           /**< Position of the beginning of the next
                                      data pack that will be unpacked       */
    size_t int_offload_array_length; /**< Total length of the common int
                                          offload array                     */
    size_t int_unpack_pos;       /**< Position of the beginning of the next
                                      data pack that will be unpacked in int
                                      array                                 */
} offload_package;

void offload_init_offload(offload_package* o, real** pack_array,
                          int** int_array);
void offload_free_offload(offload_package* o, real** pack_array,
                          int** int_array);

void offload_pack(offload_package* o, real** offload_array, real* pack_array,
                  size_t pack_length, int** int_offload_array,
                  int* int_pack_array, size_t int_pack_length);

#pragma omp declare target
void offload_unpack(offload_package* o, real* pack_array,
                    size_t pack_length, int* int_offload_array,
                    size_t int_pack_length, real** ptr, int** intptr);
#pragma omp end declare target

#endif
