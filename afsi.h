/**
 * @file afsi.c
 * @brief ASCOT Fusion Source Integrator AFSI
 */
#ifndef AFSI_H
#define AFSI_H

#include <math.h>
#include "ascot5.h"
#include "consts.h"
#include "random.h"
#include "boschhale.h"
#include "diag/dist_5D.h"

/**
 * @brief Structure for passing in 2D thermal temperature and density
 */
typedef struct {
    int n_r;          /**< number of r bins       */
    real min_r;       /**< value of lowest r bin  */
    real max_r;       /**< value of highest r bin */

    int n_phi;        /**< number of r bins       */
    real min_phi;     /**< value of lowest r bin  */
    real max_phi;     /**< value of highest r bin */

    int n_z;          /**< number of z bins       */
    real min_z;       /**< value of lowest z bin  */
    real max_z;       /**< value of highest z bin */

    real* temperature; /**< pointer to start of histogram array */
    real* density;     /**< pointer to start of histogram array */
} afsi_thermal_data;


/**
 * @brief Wrapper around input data structures
 */
typedef struct {
    int type;                        /**< Fusion reaction type */
    dist_5D_data* dist_5D;           /**< Distribution data    */
    afsi_thermal_data* dist_thermal; /**< Thermal data         */
} afsi_data;


/**
 * @brief Struct for storing 5D output data
 */
typedef struct {
    int n_r;         /**< number of r bins       */
    real min_r;      /**< value of lowest r bin  */
    real max_r;      /**< value of highest r bin */

    int n_phi;       /**< number of r bins       */
    real min_phi;    /**< value of lowest r bin  */
    real max_phi;    /**< value of highest r bin */

    int n_z;         /**< number of z bins       */
    real min_z;      /**< value of lowest z bin  */
    real max_z;      /**< value of highest z bin */

    int n_pitch;     /**< number of v_parallel bins            */
    real min_pitch;  /**< value of lowest v_parallel bin       */
    real max_pitch;  /**< value of highest v_parallel bin      */

    int n_energy;    /**< number of v_perpendicular bins       */
    real min_energy; /**< value of lowest v_perpendicular bin  */
    real max_energy; /**< value of highest v_perpendicular bin */

    real* histogram; /**< pointer to start of histogram array  */
} afsi_dist_5D;


/**
 * @brief Struct for storing 6D output data
 */
typedef struct {
    int n_r;       /**< number of r bins       */
    real min_r;    /**< value of lowest r bin  */
    real max_r;    /**< value of highest r bin */

    int n_phi;     /**< number of r bins       */
    real min_phi;  /**< value of lowest r bin  */
    real max_phi;  /**< value of highest r bin */

    int n_z;       /**< number of z bins       */
    real min_z;    /**< value of lowest z bin  */
    real max_z;    /**< value of highest z bin */

    int n_vr;      /**< number of v_parallel bins            */
    real min_vr;   /**< value of lowest v_parallel bin       */
    real max_vr;   /**< value of highest v_parallel bin      */

    int n_vphi;    /**< number of v_perpendicular bins       */
    real min_vphi; /**< value of lowest v_perpendicular bin  */
    real max_vphi; /**< value of highest v_perpendicular bin */

    int n_vz;      /**< number of v_perpendicular bins       */
    real min_vz;   /**< value of lowest v_perpendicular bin  */
    real max_vz;   /**< value of highest v_perpendicular bin */

    real* histogram;  /**< pointer to start of histogram array */
} afsi_dist_6D;

void afsi_mc(int react, int n, afsi_data* dist1, afsi_data* dist2,
             afsi_dist_5D* fusion_dist);
void afsi_create_reaction_products(int react, int n, real* v1x, real* v1y,
                                   real* v1z, real* v2x, real* v2y, real* v2z);
void afsi_test_dist(dist_5D_data* dist1);
void afsi_test_thermal();

#endif
