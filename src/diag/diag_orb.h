/**
 * @file diag_orb.h
 * @brief Header file for diag_orb.c.
 *
 * This file also contains definitions for orbit diagnostics data structures.
 */
#ifndef DIAG_ORB_H
#define DIAG_ORB_H

#include <stdio.h>
#include "../particle.h"
#include "../options.h"

#define DIAG_ORB_POINCARE 0      /**< Poincare mode flag                 */
#define DIAG_ORB_INTERVAL 1      /**< Interval mode flag                 */
#define DIAG_ORB_MAXPOINCARES 30 /**< Maximum number of Poincare planes  */

#define DIAG_ORB_FOFIELDS     16 /**< Number of coordinates in FO output     */
#define DIAG_ORB_GCFIELDS     16 /**< Number of coordinates in GC output     */
#define DIAG_ORB_MLFIELDS     11 /**< Number of coordinates in ML output     */
#define DIAG_ORB_HYBRIDFIELDS 19 /**< Number of coordinates in hybrid output */

#define DIAG_ORB_FO 1 /**< Data stored in FO mode */
#define DIAG_ORB_GC 2 /**< Data stored in GC mode */
#define DIAG_ORB_ML 3 /**< Data stored in ML mode */

DECLARE_TARGET_SIMD_UNIFORM(ang0)
real diag_orb_check_plane_crossing(real fang, real iang, real ang0);
DECLARE_TARGET_SIMD_UNIFORM(r0)
real diag_orb_check_radial_crossing(real fr, real ir, real r0);

/**
 * @brief Orbit diagnostics data struct.
 *
 * The pointers are asssigned to the offload array but only those pointers
 * are used which corresponds to the simulation mode (e.g. GC simulation).
 *
 * The marker orbit data are stored in maxpoints length chunks in the
 * offload array. The chuncks are in no particular order. Once chunk is
 * filled and the marker is still recording, new points replace the old
 * ones from the start.
 */
typedef struct{
    real* id;     /**< Marker ID                                            */
    real* mileage;/**< Time marker has been simulated for [s]               */
    real* r;      /**< Marker R coordinate [m]                              */
    real* phi;    /**< Marker phi coordinate [rad]                          */
    real* z;      /**< Marker z coordiante [m]                              */
    real* p_r;    /**< Particle momentum R component [kg m/s]               */
    real* p_phi;  /**< Particle momentum phi component [kg m/s]             */
    real* p_z;    /**< Particle momentum z component [kg m/s]               */
    real* ppar;   /**< Guiding center parallel momentum [kg m/s]            */
    real* mu;     /**< Guiding center magnetic moment [J/T]                 */
    real* zeta;   /**< Guiding center gyroangle [rad]                       */
    real* weight; /**< Marker weight [1]                                    */
    real* charge; /**< Marker charge [C]                                    */
    real* rho;    /**< Normalized poloidal flux at marker position [1]      */
    real* theta;  /**< Marker poloidal angle [rad]                          */
    real* B_r;    /**< Magnetic field R component at marker position [T]    */
    real* B_phi;  /**< Magnetic field phi component at marker position [T]  */
    real* B_z;    /**< Magnetic field z component at marker position [T]    */
    real* simmode;/**< In what simulation mode data point was recorded      */
    real* pncrid; /**< Id for the poincare plot a point corresponds to      */
    real* pncrdi; /**< Direction in which Poincare plane was crossed        */

    integer* mrk_pnt;     /**< Index of the last recorded point             */
    real* mrk_recorded;   /**< Last time (in seconds) a marker was updated  */

}diag_orb_data;

void diag_orb_init(diag_orb_data* data, sim_parameters* params, size_t nmarkers);

void diag_orb_free(diag_orb_data* data, sim_parameters* params);

void diag_orb_update_fo(diag_orb_data* data, sim_parameters* params,
                        particle_simd_fo* p_f, particle_simd_fo* p_i);

void diag_orb_update_gc(diag_orb_data* data, sim_parameters* params,
                        particle_simd_gc* p_f, particle_simd_gc* p_i);

void diag_orb_update_ml(diag_orb_data* data, sim_parameters* params,
                        particle_simd_ml* p_f, particle_simd_ml* p_i);

#endif
