/**
 * @file boozer.c
 * @brief Module for transforming between cylindrical and Boozer coordinates.
 */
#include "ascot5.h"
#include "error.h"
#include "boozer.h"

/**
 * @brief Load Boozer data and prepare parameters for offload.
 *
 * This function fills the boozer offload struct with parameters and allocates
 * and fills the offload array. Sets offload array length in the offload struct.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 *
 * @return zero if initialization succeeded.
 *
 * @todo Konsta will write this.
 */
int boozer_init_offload(boozer_offload_data* offload_data,
                        real** offload_array) {
    return 0;
}

/**
 * @brief Initialize boozer data struct on target
 *
 * @param boozerdata pointer to data struct on target
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array
 *
 * @todo Konsta will write this.
 */
void boozer_init(boozer_data* boozerdata, boozer_offload_data* offload_data,
                     real* offload_array) {

}

/**
 * @brief Boozer magic here.
 *
 * @todo This is just a dummy.
 */
a5err boozer_eval(boozer_data* boozerdata) {
    return 0;
}
