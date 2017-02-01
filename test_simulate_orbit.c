/**
 * @file test_simulate_orbit.c
 * @brief Tests that different orbit integrators yield same results
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ascot5.h"
#include "math.h"
#include "B_field.h"
#include "E_field.h"
#include "particle.h"
#include "step_fo_lf.h"
#include "step_fo_vpa.h"

int main(void) {

  /* Init background */
  B_field_offload_data offload_Bdata;
  real* offload_array;
  B_field_init_offload(&offload_Bdata, &offload_array);

  B_field_data Bdata;
  B_field_init(&Bdata, &offload_Bdata, offload_array);

  E_field_offload_data offload_Edata;
  E_field_init_offload(&offload_Edata, &offload_array);

  E_field_data Edata;
  E_field_init(&Edata, &offload_Edata, offload_array);

  real tstep = 1.e-9;
  real tsimend = 1.e-8;
  real tsim;

  /* Init a bunch of identical test particles */
  particle_simd_fo* p = (particle_simd_fo*) malloc(4*sizeof(particle_simd_fo));
  int i = 0;
  real Brpz[3];
  real Erpz[3];
  for(i=0;i<4;i++){
    p[i].r[0] = 8.01;
    p[i].phi[0] = 0.0;
    p[i].z[0] = 0.0;
    p[i].rdot[0] = 0.0;
    p[i].phidot[0] = 0.0;
    p[i].zdot[0] = 1.e7;
    p[i].mass[0] = 6.6447e-27;
    p[i].charge[0] = 1.6022e-19;
    p[i].weight[0] = 1;
    p[i].id[0] = i+1;
    p[i].running[0] = 1;
    B_field_eval_B(Brpz, p[i].r[0], p[i].phi[0], p[i].z[0], &Bdata);
    p[i].B_r[0] = Brpz[0];
    p[i].B_phi[0] = Brpz[1];
    p[i].B_z[0] = Brpz[2];
    E_field_eval_E(Erpz, p[i].r[0], p[i].phi[0], p[i].z[0], &Edata);
    p[i].E_r[0] = Erpz[0];
    p[i].E_phi[0] = Erpz[1];
    p[i].E_z[0] = Erpz[2];
  }

  /* Test leap-frog */
  tsim=0.0;
  printf("%d, %le, %le, %le, %le, %le, %le, %le\n",
	 (int) p[0].id[0], tsim, p[0].r[0], p[0].phi[0], p[0].z[0],
	 p[0].rdot[0], p[0].phidot[0], p[0].zdot[0]);
  while(tsim < tsimend){
    step_fo_lf(&p[0], tsim, tstep, &Bdata);
    tsim += tstep;

    printf("%d, %le, %le, %le, %le, %le, %le, %le\n",
	   (int) p[0].id[0], tsim, p[0].r[0], p[0].phi[0], p[0].z[0],
                p[0].rdot[0], p[0].phidot[0], p[0].zdot[0]);
  }

  /* Test volume-preserving algorithm */
  tsim=0.0;
  printf("%d, %le, %le, %le, %le, %le, %le, %le\n",
	 (int) p[1].id[0], tsim, p[1].r[0], p[1].phi[0], p[1].z[0],
	 p[1].rdot[0], p[1].phidot[0], p[1].zdot[0]);
  while(tsim < tsimend){
    step_fo_vpa(&p[1], tsim, tstep, &Bdata, &Edata);
    tsim += tstep;

    printf("%d, %le, %le, %le, %le, %le, %le, %le\n",
	   (int) p[1].id[0], tsim, p[1].r[0], p[1].phi[0], p[1].z[0],
                p[1].rdot[0], p[1].phidot[0], p[1].zdot[0]);
  }

  /* Test RK4 */


  /* Test Cash-Karp */
  
}

/*
void test_simulate_orbit_new_particle(particle_simd_fo porig){
  p[i].r[0] = 8.01;
  p[i].phi[0] = 0.0;
  p[i].z[0] = 0.0;
  p[i].rdot[0] = 0.0;
  p[i].phidot[0] = 0.0;
  p[i].zdot[0] = 1.e7;
  p[i].mass[0] = 6.6447e-27;
  p[i].charge[0] = 1.6022e-19;
  p[i].weight[0] = 1;
  p[i].id[0] = i+1;
  p[i].running[0] = 1;
  
}
*/
