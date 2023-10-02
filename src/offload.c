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
 *
 * Since input data may contain both floats and integers, we store those
 * separately (i.e. there is one common offload array for floats and one common
 * offload array for integers).
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
 * @param int_offload_array pointer to int packarray
 */
void offload_init_offload(offload_package* o, real** offload_array,
                          int** int_offload_array) {
    *offload_array     = NULL;
    *int_offload_array = NULL;
    o->offload_array_length     = 0;
    o->int_offload_array_length = 0;
    o->unpack_pos     = 0;
    o->int_unpack_pos = 0;
}


/**
 * @brief Free offload array and set offload_package to clean state
 *
 * @param o offload package
 * @param offload_array pointer to offload_array
 * @param int_offload_array pointer to int packarray
 */
void offload_free_offload(offload_package* o, real** offload_array,
                          int** int_offload_array) {
    free(*offload_array);
    free(*int_offload_array);
    *offload_array     = NULL;
    *int_offload_array = NULL;
    o->offload_array_length     = 0;
    o->int_offload_array_length = 0;
    o->unpack_pos     = 0;
    o->int_unpack_pos = 0;
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
 * @param int_offload_array pointer to offload array where integer data will be
 *        packed
 * @param int_pack_array pack array containing the integer data to be appended
 * @param int_pack_length length of the int_offload_array
 */
void offload_pack(offload_package* o, real** offload_array, real* pack_array,
                  size_t pack_length, int** int_offload_array,
                  int* int_pack_array, size_t int_pack_length) {

    /* Float array */
    if( pack_length > 0 ) {
        size_t new_length = o->offload_array_length + pack_length;
        real* new_array = (real*) malloc(new_length * sizeof(real));

        if(o->offload_array_length > 0) {
            memcpy(new_array, *offload_array,
                   o->offload_array_length*sizeof(real));
        }

        memcpy(new_array+o->offload_array_length, pack_array,
               pack_length*sizeof(real));

        free(*offload_array);

        *offload_array = new_array;
        o->offload_array_length = new_length;
    }

    /* Int array */
    if( int_pack_length > 0 ) {
        size_t int_new_length = o->int_offload_array_length + int_pack_length;
        int* int_new_array = (int*) malloc(int_new_length * sizeof(int));

        if(o->int_offload_array_length > 0) {
            memcpy(int_new_array, *int_offload_array,
                   o->int_offload_array_length*sizeof(int));
        }

        memcpy(int_new_array+o->int_offload_array_length, int_pack_array,
               int_pack_length*sizeof(int));

        free(*int_offload_array);

        *int_offload_array = int_new_array;
        o->int_offload_array_length = int_new_length;
    }
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
 * @param int_offload_array int offload array
 * @param int_pack_length
 * @param ptr pointer where the unpacked data begins
 * @param intptr pointer where the unpacked int data begins
 */
void offload_unpack(offload_package* o, real* offload_array,
                    size_t pack_length, int* int_offload_array,
                    size_t int_pack_length, real** ptr, int** intptr) {
    *ptr = offload_array + o->unpack_pos;
    o->unpack_pos += pack_length;

    *intptr = int_offload_array + o->int_unpack_pos;
    o->int_unpack_pos += int_pack_length;
}
