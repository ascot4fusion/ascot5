/**
 * @file orbit_write.c
 * @brief !!THIS MODULE IS REDUNDANT!! Functions to write particle and guiding center information. 
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "ascot5.h"
#include "consts.h"
#include "phys_orbit.h"
#include "particle.h"
#include "B_field.h"
#include "orbit_write.h"

void write_particle(FILE* out, particle* p){
    fprintf(out, "%d, %le, %le, %le, %le, %le, %le, %le\n",
	   (int) p->id, p->time, p->r, p->phi, p->z,
	   p->v_r, p->v_phi, p->v_z);
    
}


void write_guidingcenter(FILE* out, particle* p){

    
}

void write_fo_as_particle(FILE* out, particle_simd_fo* p){
    fprintf(out, "%d, %le, %le, %le, %le, %le, %le, %le\n",
	    (int) p->id[0], p->time[0], p->r[0], p->phi[0], p->z[0],
	    p->rdot[0], p->r[0]*p->phidot[0], p->zdot[0]);
}

void write_fo_as_guidingcenter(FILE* out, particle_simd_fo* p, B_field_data* Bdata){
    real gcpos[5];
    real B_dB[12];
    B_field_eval_B_dB(B_dB, p->r[0], p->phi[0], p->z[0], Bdata);

    real gamma = phys_gammaprtv(sqrt(p->rdot[0]*p->rdot[0] + pow(p->phidot[0]*p->r[0],2) + p->zdot[0]*p->zdot[0]));
    phys_prttogc(p->mass[0], p->charge[0], p->r[0], p->phi[0], p->z[0], 
		 gamma*p->mass[0]*p->rdot[0], gamma*p->mass[0]*p->phidot[0]*p->r[0], gamma*p->mass[0]*p->zdot[0], B_dB, gcpos);

    B_field_eval_B_dB(B_dB, gcpos[0], gcpos[1], gcpos[2], Bdata);
    gamma = phys_gammagcp(p->mass[0], gcpos[3], gcpos[4]);
    gcpos[3] = gcpos[3]/(p->mass[0]*gamma);

    fprintf(out, "%d, %le, %le, %le, %le, %le, %le\n",
	    (int) p->id[0], p->time[0], gcpos[0], gcpos[1], gcpos[2],
	    gcpos[3], gcpos[4]);

    
}

void write_gc_as_particle(FILE* out, particle* p){

}

void write_gc_as_guidingcenter(FILE* out, particle_simd_gc* p){
    fprintf(out, "%d, %le, %le, %le, %le, %le, %le\n",
	    (int) p->id[0], p->time[0], p->r[0], p->phi[0], p->z[0],
	    p->vpar[0], p->mu[0]);
}

void writeorbit_store_gc2guidingcenter(particle_simd_gc p, writeorbit_guidingcenter* w, int* write, int Nwrite, int slot){
    if(Nwrite == 0){
	return;
    }

    
    int i;
    for(i=0; i < NSIMD; i=i+1){
	int k = write[i];
	if(k >= 0){
	    w->id[slot+k]   = p.id[i];
	    w->time[slot+k] = p.time[i];
	    w->r[slot+k]    = p.r[i];
	    w->phi[slot+k]  = p.phi[i];
	    w->z[slot+k]    = p.z[i];
	    w->vpar[slot+k] = p.vpar[i];
	    w->mu[slot+k]   = p.mu[i];
	}
	
    }
}
