/**
 * @file copytogpu.h
 * @brief Header file for copytogpu.c
 */
#ifndef COPYTOGPU_H
#define COPYTOGPU_H
void simulate_fo_fixed_copy_to_gpu(sim_data* sim, particle_simd_fo *p_ptr, particle_simd_fo *p0_ptr, real* hin, real* rnd);

void simulate_fo_fixed_copy_from_gpu(sim_data* sim, particle_simd_fo *p_ptr);

#endif
