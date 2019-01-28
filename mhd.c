/**
 * @file mhd.c
 * @brief Module for evaluating MHD parameters.
 */
#include "ascot5.h"
#include "error.h"
#include "mhd.h"

/**
 * @brief Load MHD data and prepare parameters for offload.
 *
 * This function fills the MHD offload struct with parameters and allocates
 * and fills the offload array. Sets offload array length in the offload struct.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 *
 * @return zero if initialization succeeded.
 *
 * @todo Konsta will write this.
 */
int mhd_init_offload(mhd_offload_data* offload_data,
                     real** offload_array) {
    return 0;
}

/**
 * @brief Initialize MHD data struct on target
 *
 * @param mhddata pointer to data struct on target
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array
 *
 * @todo Konsta will write this.
 */
void mhd_init(mhd_data* MHDdata, mhd_offload_data* offload_data,
              real* offload_array) {

}

/**
 * @brief MHD magic here.
 *
 * @todo This is just a dummy.
 */
a5err mhd_eval(mhd_data* mhddata) {
    return 0;
}
