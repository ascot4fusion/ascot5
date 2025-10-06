/**
 * @file diag_orb.h
 * @brief Header file for diag_orb.c.
 *
 * This file also contains definitions for orbit diagnostics data structures.
 */
#ifndef DIAG_ORB_H
#define DIAG_ORB_H

#include <stdio.h>
#include "particle.h"
#include "options.h"
#include "diag.h"

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


void diag_orb_init(diag_orb_data* data, sim_parameters* params, size_t nmarkers);

void diag_orb_free(diag_orb_data* data, sim_parameters* params);

void diag_orb_update_fo(diag_orb_data* data, sim_parameters* params,
                        particle_simd_fo* p_f, particle_simd_fo* p_i);

void diag_orb_update_gc(diag_orb_data* data, sim_parameters* params,
                        particle_simd_gc* p_f, particle_simd_gc* p_i);

void diag_orb_update_ml(diag_orb_data* data, sim_parameters* params,
                        particle_simd_ml* p_f, particle_simd_ml* p_i);

#endif
