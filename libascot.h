#ifndef LIBASCOT_H_
#define LIBASCOT_H_

#include "ascot5.h"
#include "simulate.h"

void libascot_B_field_eval_B_dB(sim_offload_data* sim0, real* B_offload_array,
                                int Neval, real* R, real* phi, real* z, real* t,
                                real* BR, real* Bphi, real* Bz,
                                real* BR_dR, real* BR_dphi, real* BR_dz,
                                real* Bphi_dR, real* Bphi_dphi, real* Bphi_dz,
                                real* Bz_dR, real* Bz_dphi, real* Bz_dz);

#endif
