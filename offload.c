/**
 * @file offload.c
 * @brief Offload functions
 *
 * An offload package that is used to collect all offload arrays into a single
 * one so that offloading can be done nice and clean. offload_pack() takes an
 * offload data array and appends it to the common offload array.
 * offload_unpack() returns a pointer at the start of the common offload array.
 * offload_unpack() must be called in the same order as offload_pack() was
 * called when data was packaged.
 */
#include <stdlib.h>
#include <string.h>
#include "ascot5.h"
#include "offload.h"

/**
 * @brief Initialize offload package
 *
 * Sets offload_package.pack_array_length and offload_package.unpack_pos to
 * zero and assigns NULL to offload_array pointer.
 *
 * @param o uninitialized offload package
 * @param offload_array pointer to packarray
 */
void offload_init_offload(offload_package* o, real** offload_array) {
    *offload_array = NULL;
    o->offload_array_length = 0;
    o->unpack_pos = 0;
}

/**
 * @brief Free offload array and set offload_package to clean state
 *
 * @param o offload package
 * @param offload_array pointer to offload_array
 */
void offload_free_offload(offload_package* o, real** offload_array) {
    free(*offload_array);
    *offload_array = NULL;
    o->offload_array_length = 0;
    o->unpack_pos = 0;
}

/**
 * @brief Pack an offload array to package array
 *
 * Reallocates the offload array so that the new data from the package array can
 * be appended. Once the data has been appended, the old offload array is freed.
 * The new length of the offload array is stored to offload_package struct.
 *
 * @param o offload package
 * @param offload_array pointer to offload array where data will be packed
 * @param pack_array pack array containing the data to be appended
 * @param pack_length length of the offload_array
 */
void offload_pack(offload_package* o, real** offload_array, real* pack_array,
                  size_t pack_length) {

    size_t new_length = o->pack_array_length + pack_length;
    real* new_array = (real*) malloc(new_length * sizeof(real));

    if(o->pack_array_length > 0) {
        memcpy(new_array, *offload_array, o->pack_array_length*sizeof(real));
    }

    memcpy(new_array+o->pack_array_length, pack_array,
           pack_length*sizeof(real));

    free(*offload_array);

    *offload_array = new_array;
    o->pack_array_length = new_length;
}

/**
 * @brief Unpack offload array from the package
 *
 * The data is not actually moved anywhere. Unpacking means that pointers to the
 * start of packed offload arrays are returned. Unpacking must be done in the
 * same order as data was packed.
 *
 * @param o offload package
 * @param offload_array offload array
 * @param pack_length length of the data that is unpacked
 */
real* offload_unpack(offload_package* o, real* offload_array,
                     size_t pack_length) {
    real* ptr = offload_array + o->unpack_pos;

    o->unpack_pos += pack_length;

    return ptr;
}
