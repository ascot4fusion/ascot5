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
#include "orbit_write.h"
#include "step_fo_lf.h"
#include "step_fo_vpa.h"
#include "step_gc_rk4.h"
#include "step_gc_cashkarp.h"

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

    real tstep[NSIMD];
    tstep[0] = 1.e-10;
    real tsimend = 1.e-6;

    /* Init a bunch of identical test particles */
    particle p;
    particle_simd_fo p_fo;
    particle_simd_gc p_gc;

    p.r = 8.01;
    p.phi = 0.0;
    p.z = 0.0;
    p.rdot = 2.e4;
    p.phidot = 0.0;
    p.zdot = 1.e4;
    p.mass = 6.6447e-27;
    p.charge = 1.6022e-19;
    p.weight = 1;
    p.id = 1;
    p.running = 1;
  
    /* Output files */
    FILE* f_particle = fopen("orbits_particle.test","w");
    FILE* f_guidingcenter = fopen("orbits_guidingcenter.test","w");

    /* Test leap-frog */
    p.id = 1;
    p.time = 0.0;
    particle_to_fo(&p, 0, &p_fo, 0, &Bdata, &Edata);
    write_fo_as_particle(f_particle, &p_fo);
    write_fo_as_guidingcenter(f_guidingcenter, &p_fo, &Bdata);
    while(p_fo.time[0] < tsimend){
	step_fo_lf(&p_fo, p_fo.time[0], tstep[0], &Bdata);
	p_fo.time[0] += tstep[0];
	write_fo_as_particle(f_particle, &p_fo);
	write_fo_as_guidingcenter(f_guidingcenter, &p_fo, &Bdata);
    }

    /* Test volume-preserving algorithm */
    p.id = 2;
    p.time = 0.0;
    particle_to_fo(&p, 0, &p_fo, 0, &Bdata, &Edata);
    write_fo_as_particle(f_particle, &p_fo);
    write_fo_as_guidingcenter(f_guidingcenter, &p_fo, &Bdata);
    while(p_fo.time[0] < tsimend){
	step_fo_vpa(&p_fo, p_fo.time[0], tstep[0], &Bdata, &Edata);
	p_fo.time[0] += tstep[0];
	write_fo_as_particle(f_particle, &p_fo);
	write_fo_as_guidingcenter(f_guidingcenter, &p_fo, &Bdata);
    }

    /* Test RK4 */
    p.id = 3;
    p.time = 0.0;
    particle_to_gc(&p, 0, &p_gc, 0, &Bdata);
    write_gc_as_guidingcenter(f_guidingcenter, &p_gc);
    while(p_gc.time[0] < tsimend){
	step_gc_rk4(&p_gc, p_fo.time[0], tstep[0], &Bdata, &Edata);
	p_gc.time[0] += tstep[0];
	write_gc_as_guidingcenter(f_guidingcenter, &p_gc);
    }

    /* Test Cash-Karp */
    real tol = 1.0;
    real tnext[NSIMD];
    p.id = 4;
    p.time = 0.0;
    particle_to_gc(&p, 0, &p_gc, 0, &Bdata);
    write_gc_as_guidingcenter(f_guidingcenter, &p_gc);
    while(p_gc.time[0] < tsimend){
	step_gc_cashkarp(&p_gc, p_fo.time, tstep, tnext, tol, &Bdata, &Edata);
	p_gc.time[0] += tstep[0];
	write_gc_as_guidingcenter(f_guidingcenter, &p_gc);
    }
    
  
    /* Done! Close files and clean */
    fclose(f_particle);
    fclose(f_guidingcenter);
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
