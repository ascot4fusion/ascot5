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

#define DIAG_ORB_POINCARE 0      /**< Poincare mode flag                 */
#define DIAG_ORB_INTERVAL 1      /**< Interval mode flag                 */
#define DIAG_ORB_MAXPOINCARES 30 /**< Maximum number of Poincare planes  */

#define DIAG_ORB_FOFIELDS     15 /**< Number of coordinates in FO output     */
#define DIAG_ORB_GCFIELDS     15 /**< Number of coordinates in GC output     */
#define DIAG_ORB_MLFIELDS     10 /**< Number of coordinates in ML output     */
#define DIAG_ORB_HYBRIDFIELDS 18 /**< Number of coordinates in hybrid output */

/**
 * @brief Orbit diagnostics offload data struct.
 */
typedef struct{
    int record_mode;    /**< Defines what fields are initialized           */
    int mode;           /**< Defines condition for recording markers       */
    int Npnt;           /**< Maximum number of points to keep recorded     */
    int Nmrk;           /**< Number of markers to record                   */
    int Nfld;           /**< Number of fields the record contains          */
    real writeInterval; /**< Interval at which markers are recorded        */
    int ntoroidalplots; /**< Number of toroidal Poincare planes            */
    int npoloidalplots; /**< Number of toroidal Poincare planes            */
    real toroidalangles[DIAG_ORB_MAXPOINCARES]; /**< Toroidal plane angles */
    real poloidalangles[DIAG_ORB_MAXPOINCARES]; /**< Poloidal plane angles */

}diag_orb_offload_data;

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
    real* time;   /**< Marker time [s]                                      */
    real* r;      /**< Marker R coordinate [m]                              */
    real* phi;    /**< Marker phi coordinate [rad]                          */
    real* z;      /**< Marker z coordiante [m]                              */
    real* rdot;   /**< Particle dR/dt [m/s]                                 */
    real* phidot; /**< Particle dphi/dt [rad]                               */
    real* zdot;   /**< Particle dz/dt [m/s]                                 */
    real* vpar;   /**< Guiding center parallel velocity [m/s]               */
    real* mu;     /**< Guiding center magnetic moment [J/T]                 */
    real* theta;  /**< Guiding center gyroangle [rad]                       */
    real* weight; /**< Marker weight [1]                                    */
    real* charge; /**< Marker charge [C]                                    */
    real* rho;    /**< Normalized poloidal flux at marker position [1]      */
    real* pol;    /**< Marker poloidal angle [rad]                          */
    real* B_r;    /**< Magnetic field R component at marker position [T]    */
    real* B_phi;  /**< Magnetic field phi component at marker position [T]  */
    real* B_z;    /**< Magnetic field z component at marker position [T]    */
    real* pncrid; /**< Id for the poincare plot a point corresponds to      */

    integer* mrk_pnt;     /**< Index of the last recorded point             */
    real* mrk_recorded;   /**< Last time (in seconds) a marker was updated  */

    int mode;           /**< Defines condition for recording markers        */
    int Npnt;           /**< Maximum number of points to keep recorded      */
    int Nmrk;           /**< Number of markers to record                    */
    real writeInterval; /**< Interval at which markers are recorded         */
    int ntoroidalplots; /**< Number of toroidal Poincare planes             */
    int npoloidalplots; /**< Number of toroidal Poincare planes             */
    real toroidalangles[DIAG_ORB_MAXPOINCARES]; /**< Toroidal plane angles  */
    real poloidalangles[DIAG_ORB_MAXPOINCARES]; /**< Poloidal plane angles  */

}diag_orb_data;

#pragma omp declare target
void diag_orb_init(diag_orb_data* data, diag_orb_offload_data* offload_data,
                   real* offload_array);

void diag_orb_free(diag_orb_data* data);

void diag_orb_update_fo(diag_orb_data* data,
                        particle_simd_fo* p_f, particle_simd_fo* p_i);

void diag_orb_update_gc(diag_orb_data* data,
                        particle_simd_gc* p_f, particle_simd_gc* p_i);

void diag_orb_update_ml(diag_orb_data* data,
                        particle_simd_ml* p_f, particle_simd_ml* p_i);
#pragma omp end declare target

#endif
