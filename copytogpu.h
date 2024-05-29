/**
 * @file copytogpu.h
 * @brief Header file for copytogpu.c
 */
#ifndef COPYTOGPU_H
#define COPYTOGPU_H
void simulate_fo_fixed_copy_to_gpu(sim_data* sim, particle_simd_fo *p_ptr, particle_simd_fo *p0_ptr, B_field_data* Bdata, E_field_data* Edata, particle_loc*  p_loc, real* hin, real* rnd, int n_queue_size );

void simulate_fo_fixed_copy_from_gpu(sim_data* sim, particle_simd_fo *p_ptr, int n_queue_size);

#endif
