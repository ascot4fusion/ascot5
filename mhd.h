/**
 * @file mhd.h
 * @brief Header file for mhd.c
 *
 * Contains a list declaring all mhd_types, and declaration of
 * mhd_offload_data and mhd_data structs.
 */
#ifndef MHD_H
#define MHD_H

#include "ascot5.h"
#include "error.h"
#include "B_field.h"
#include "boozer.h"
#include "mhd/mhd_stat.h"
#include "mhd/mhd_nonstat.h"

/** @brief includemode parameter to include all modes (default) */
#define MHD_INCLUDE_ALL -1

/**
 * @brief MHD input types
 *
 * MHD types are used in the MHD interface (mhd.c) to direct function calls to
 * correct MHD instances. Each MHD instance must have a corresponding type.
 */
typedef enum mhd_type {
    mhd_type_stat,   /**< MHD where mode amplitude does not depend on time */
    mhd_type_nonstat /**< MHD where mode amplitude depends on time         */
} mhd_type;

/**
 * @brief MHD offload data
 *
 * This struct holds data necessary for offloading. The struct is initialized in
 * mhd_init_offload().
 *
 * The intended usage is that only single offload data is used at the time, and
 * the type of the data is declared with the "type" field.
 */
typedef struct {
    mhd_type type;                   /**< MHD type wrapped by this struct     */
    mhd_stat_offload_data stat;      /**< Stat field or NULL if not active    */
    mhd_nonstat_offload_data nonstat;/**< Nonstat field or NULL if not active */
    int offload_array_length;        /**< Allocated offload array length      */
} mhd_offload_data;

/**
 * @brief MHD simulation data
 *
 * This struct holds data necessary for simulation. The struct is initialized
 * from the mhd_offload_data in mhd_init().
 *
 * The intended usage is that only single mhd_data is used at the time, and
 * the type of the data is declared with the "type" field.
 */
typedef struct {
    mhd_type type;            /**< MHD type wrapped by this struct     */
    mhd_stat_data stat;       /**< Stat field or NULL if not active    */
    mhd_nonstat_data nonstat; /**< Nonstat field or NULL if not active */
} mhd_data;

int mhd_init_offload(mhd_offload_data* offload_data,
                     real** offload_array);
void mhd_free_offload(mhd_offload_data* offload_data,
                      real** offload_array);

#pragma omp declare target
int mhd_init(mhd_data* data, mhd_offload_data* offload_data,
             real* offload_array);
#pragma omp declare simd uniform(boozerdata, mhddata, Bdata, includemode)
a5err mhd_eval(real mhd_dmhd[10], real r, real phi, real z, real t,
               int includemode, boozer_data* boozerdata, mhd_data* mhddata,
               B_field_data* Bdata);
#pragma omp declare simd uniform(boozerdata, mhddata, Bdata, pertonly,\
                                 includemode)
a5err mhd_perturbations(real pert_field[7], real r, real phi, real z,
                        real t, int pertonly, int includemode,
                        boozer_data* boozerdata, mhd_data* mhddata,
                        B_field_data* Bdata);

#pragma omp end declare target
#endif
