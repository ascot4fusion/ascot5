/**
 * @file orbit_write.h
 * @brief Functions to write particle and guiding center information.
 */
#ifndef DIAG_ORB_H
#define DIAG_ORB_H

#include <stdio.h>
#include "particle.h"

#define DIAG_ORB_POINCARE 0
#define DIAG_ORB_ORBIT 1
#define DIAG_ORB_WRITELAST 2
#define DIAG_ORB_MAXPOINCARES 30

typedef struct{ 
    integer id;
    real time;
    real r;
    real phi;
    real z;
    real rdot;
    real phidot;
    real zdot;
    real weight;
    real charge;
    real mass;
    real B_r;
    real B_phi;
    real B_z;

} diag_orb_dat_fo;

typedef struct{ 
    integer id;
    real time;
    real r;
    real phi;
    real z;
    real vpar;
    real mu;
    real theta;
    real weight;
    real charge;
    real mass;
    real B_r;
    real B_phi;
    real B_z;

} diag_orb_dat_gc;

typedef struct{ 
    integer id;
    real time;
    real r;
    real phi;
    real z;
    real weight;
    real B_r;
    real B_phi;
    real B_z;

} diag_orb_dat_ml;

typedef enum diag_orb_dat_type {
    diag_orb_type_fo, diag_orb_type_gc, diag_orb_type_ml
} diag_orb_dat_type;


typedef struct diag_orb_dat {
    union {
	diag_orb_dat_fo fo;
	diag_orb_dat_gc gc;
	diag_orb_dat_ml ml;
    };

    struct diag_orb_dat* prev;
    struct diag_orb_dat* next;

    int poincareId;
} diag_orb_dat;


typedef struct{

    /* Options */
    real writeInterval;
    int mode;
    int ntoroidalplots;
    real toroidalangles[DIAG_ORB_MAXPOINCARES];
    int npoloidalplots;
    real poloidalangles[DIAG_ORB_MAXPOINCARES];
    int writeNlast;

}diag_orb_offload_data;

typedef struct{


    /* Data storage */
    diag_orb_dat* writelist;
    diag_orb_dat* poincarelist;
    diag_orb_dat** lastlist;
    int size;
    diag_orb_dat_type type;

    /* Options */
    real writeInterval;
    int mode;
    int ntoroidalplots;
    real toroidalangles[DIAG_ORB_MAXPOINCARES];
    int npoloidalplots;
    real poloidalangles[DIAG_ORB_MAXPOINCARES];
    int writeNlast;

    /* Particle specific data */
    int particleId[NSIMD];
    real prevWriteTime[NSIMD];
    
}diag_orb_data;

#pragma omp declare target

void diag_orb_init_offload(diag_orb_offload_data* data);

void diag_orb_init(diag_orb_data* data, diag_orb_offload_data* offload_data);

void diag_orb_update_gc(particle_simd_gc* p_f, particle_simd_gc* p_i, diag_orb_data* data);

void diag_orb_update_fo(particle_simd_fo* p_f, particle_simd_fo* p_i, diag_orb_data* data);

void diag_orb_update_ml(particle_simd_ml* p_f, particle_simd_ml* p_i, diag_orb_data* data);

void diag_orb_clean(diag_orb_data* data);
#pragma omp end declare target

#endif
