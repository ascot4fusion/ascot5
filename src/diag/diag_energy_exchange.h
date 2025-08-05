/**
 * @file diag_energy_exchange.h
 * @brief Header file for the energy exchange diagnostic module.
 */

#ifndef DIAG_ENERGY_EXCHANGE_H
#define DIAG_ENERGY_EXCHANGE_H

#include "../ascot5.h"
#include "../mhd.h"
#include "../boozer.h"
#include "../spline/interp.h"
#include "../particle.h"
#include "../B_field.h"
#include "../plasma.h"


#ifdef __cplusplus
extern "C" {
#endif
/**
 * @brief energy exchange diagnostic data structure.
 */
typedef struct {
    real* S1;      /** Source term for mode evolution */
    real* S2;      /** Source term for mode evolution */
    real* dEnergy; /** Change in nergy. */
    int* is_mode_evol; /** Whether to evolve the mode */
    int n_modes;   /** Number of modes. */
    int nprt;      /** Number of particles. */

    // MHD mode parameters and Boozer coordinates.
    mhd_stat_data* mhd_data;
    boozer_data* boozerdata;
    B_field_data* B_data;
    plasma_data* plasma_data; /** Plasma data for displacement */
    real thrmass; /** Mass of the thermal particle */

    // Internal variables.
    int enabled;

} diag_energy_exchange_data;

void diag_energy_exchange_init(diag_energy_exchange_data* data, 
                               int nprt, boozer_data* boozer_data,
                               B_field_data* B_data,
                               mhd_stat_data* mhd_data, 
                               plasma_data* plasma_data);
void diag_energy_exchange_update_nprt(diag_energy_exchange_data* data, 
                                        int nprt);
void diag_energy_exchange_free(diag_energy_exchange_data* data);
void diag_energy_exchange_offload(diag_energy_exchange_data* data);
void diag_energy_exchange_onload(diag_energy_exchange_data* data);
void diag_energy_exchange_update_fo(diag_energy_exchange_data* data, 
                                    particle_simd_fo* p_f, 
                                    particle_simd_fo* p_i);
void diag_energy_exchange_update_gc(diag_energy_exchange_data* data, 
                                     particle_simd_gc* p_f, 
                                     particle_simd_gc* p_i);

void diag_energy_exchange_compact(diag_energy_exchange_data* data, 
                                  real* S1, real* S2, int clear);

#ifdef __cplusplus
}
#endif

#endif 