/**
 * @file orbit_write.h
 * @brief Functions to write particle and guiding center information.
 */
#ifndef DIAG_ORB_H
#define DIAG_ORB_H

#include <stdio.h>
#include "particle.h"
#include "B_field.h"

#define DIAG_ORB_SLOTS 100000
#define DIAG_ORB_ORBIT 1
#define DIAG_ORB_POINCARE 2

typedef struct{ 
    int id[DIAG_ORB_SLOTS];
    real time[DIAG_ORB_SLOTS];
    real r[DIAG_ORB_SLOTS];
    real phi[DIAG_ORB_SLOTS];
    real z[DIAG_ORB_SLOTS];
    real vpar[DIAG_ORB_SLOTS];
    real mu[DIAG_ORB_SLOTS];

    int size;
} diag_orb_gcarray;

typedef struct diag_orb_gclist{ 
    int id;
    real time;
    real r;
    real phi;
    real z;
    real vpar;
    real mu;

    struct diag_orb_gclist * prev;
} diag_orb_gclist;

typedef struct{ 
    int id[DIAG_ORB_SLOTS];
    real time[DIAG_ORB_SLOTS];
    real r[DIAG_ORB_SLOTS];
    real phi[DIAG_ORB_SLOTS];
    real z[DIAG_ORB_SLOTS];
    real rdot[DIAG_ORB_SLOTS];
    real phidot[DIAG_ORB_SLOTS];
    real zdot[DIAG_ORB_SLOTS];

    int size;
} diag_orb_foarray;

typedef struct diag_orb_folist{ 
    int id;
    real time;
    real r;
    real phi;
    real z;
    real rdot;
    real phidot;
    real zdot;

    struct diag_orb_folist * prev;
} diag_orb_folist;

typedef struct{ 
    int id[DIAG_ORB_SLOTS];
    real distance[DIAG_ORB_SLOTS];
    real r[DIAG_ORB_SLOTS];
    real phi[DIAG_ORB_SLOTS];
    real z[DIAG_ORB_SLOTS];

    int size;
} diag_orb_mlarray;

typedef struct diag_orb_mllist{ 
    int id;
    real distance;
    real r;
    real phi;
    real z;

    struct diag_orb_mllist * prev;
} diag_orb_mllist;

typedef struct{

    /* Essentials */
    int storewhat;
    int islist;
    int markertype;// 0: fo, 1: gc, 2: ml

    /* Options */
    real writeInterval;

}diag_orb_offload_data;

typedef struct{

    /* Essentials */
    int storewhat;
    int islist;
    int markertype;// 0: fo, 1: gc, 2: ml

    /* Data storage */
    diag_orb_gclist* gclist;
    diag_orb_gcarray* gcarray;
    diag_orb_folist* folist;
    diag_orb_foarray* foarray;
    diag_orb_mllist* mllist;
    diag_orb_mlarray* mlarray;

    /* Options */
    real writeInterval;

    /* Particle specific data */
    int particleId[NSIMD];
    real lastWriteTime[NSIMD];
    
}diag_orb_data;

#pragma omp declare target
void diag_orb_writegc(FILE* out, particle_simd_gc* p);

void diag_orb_writegcarray(FILE* out, diag_orb_gcarray* w);

void diag_orb_writegclist(FILE* out, diag_orb_gclist** w);

void diag_orb_storegcarray(particle_simd_gc* p, diag_orb_gcarray* w, int* write);

void diag_orb_storegclist(particle_simd_gc* p, diag_orb_gclist** winout, int* write);

diag_orb_gcarray* diag_orb_initgcarray();

diag_orb_gclist* diag_orb_initgclist();

void diag_orb_poincaregc(particle_simd_gc* p_f, particle_simd_gc* p_i, diag_orb_data* data, int* write);

void diag_orb_orbitgc(particle_simd_gc* p_f, particle_simd_gc* p_i, diag_orb_data* data, int* write);

void diag_orb_updategc(particle_simd_gc* p_f, particle_simd_gc* p_i, diag_orb_data* data);

void diag_orb_writefo(FILE* out, particle_simd_fo* p);

void diag_orb_writefoarray(FILE* out, diag_orb_foarray* w);

void diag_orb_writefolist(FILE* out, diag_orb_folist** w);

void diag_orb_storefoarray(particle_simd_fo* p, diag_orb_foarray* w, int* write);

void diag_orb_storefolist(particle_simd_fo* p, diag_orb_folist** winout, int* write);

diag_orb_foarray* diag_orb_initfoarray();

diag_orb_folist* diag_orb_initfolist();

void diag_orb_poincarefo(particle_simd_fo* p_f, particle_simd_fo* p_i, diag_orb_data* data, int* write);

void diag_orb_orbitfo(particle_simd_fo* p_f, particle_simd_fo* p_i, diag_orb_data* data, int* write);

void diag_orb_updatefo(particle_simd_fo* p_f, particle_simd_fo* p_i, diag_orb_data* data);

void diag_orb_writeml(FILE* out, particle_simd_ml* p);

void diag_orb_writemlarray(FILE* out, diag_orb_mlarray* w);

void diag_orb_writemllist(FILE* out, diag_orb_mllist** w);

void diag_orb_storemlarray(particle_simd_ml* p, diag_orb_mlarray* w, int* write);

void diag_orb_storemllist(particle_simd_ml* p, diag_orb_mllist** winout, int* write);

diag_orb_mlarray* diag_orb_initmlarray();

diag_orb_mllist* diag_orb_initmllist();

void diag_orb_poincareml(particle_simd_ml* p_f, particle_simd_ml* p_i, diag_orb_data* data, int* write);

void diag_orb_orbitml(particle_simd_ml* p_f, particle_simd_ml* p_i, diag_orb_data* data, int* write);

void diag_orb_updateml(particle_simd_ml* p_f, particle_simd_ml* p_i, diag_orb_data* data);

void diag_orb_write(diag_orb_data* data);

diag_orb_data* diag_orb_init();

void diag_orb_clean(diag_orb_data* data);
#pragma omp end declare target

#endif
