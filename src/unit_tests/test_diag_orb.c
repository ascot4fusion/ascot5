/*
 * test_diag_orb.c
 *
 *  Created on: Oct 12, 2020
 *      Author: sjjamsa
 */

#include "../diag/diag_orb.h"
#include "../consts.h"
#include "../simulate.h"

#include <stdio.h>

#define TOLERANCE 1.0e-12


int test_radial_crossings(){

    int nFails;

    nFails = 0;


    real correct,fr,ir,r0,got;

#define RADIAL_CROSSING(FR,IR,R0,CORRECT,GOT) \
    got = diag_orb_check_radial_crossing( FR, IR, R0);\
    if ( fabs( CORRECT - GOT)  > TOLERANCE ) {\
        nFails += 1;\
        printf(" Radial: ini=%g r0=%g final=%g FAIL: got=%g, expected=%g\n    ---> FAIL!\n",IR,R0,FR,GOT,CORRECT);\
    }\
    else {\
    printf(" Radial: ini=%g r0=%g final=%g got=%g, expected=%g ---> OK \n",IR,R0,FR,GOT,CORRECT);\
    }

    printf("Testing radial crossings.\n");

    ir = 0.1; fr = 0.3; r0 = 0.2; correct = 0.5; RADIAL_CROSSING(fr,ir,r0,correct,got)
    ir = 0.3; fr = 0.1; r0 = 0.2; correct = 0.5; RADIAL_CROSSING(fr,ir,r0,correct,got)
    ir = 0.4; fr = 0.3; r0 = 0.2; correct = 0.0; RADIAL_CROSSING(fr,ir,r0,correct,got)
    ir = 0.9; fr = 1.1; r0 = 1.0; correct = 0.5; RADIAL_CROSSING(fr,ir,r0,correct,got)
    ir = 0.0; fr = 0.2; r0 = 0.0; correct = 0.0; RADIAL_CROSSING(fr,ir,r0,correct,got)
    ir = 0.0; fr = 0.2; r0 = 0.2; correct = 0.0; RADIAL_CROSSING(fr,ir,r0,correct,got)
    ir = 0.0; fr = 0.200000000000001; r0 = 0.2; correct = 1.0; RADIAL_CROSSING(fr,ir,r0,correct,got)
    ir = 0.0; fr = 0.0; r0 = 0.2; correct = 0.0; RADIAL_CROSSING(fr,ir,r0,correct,got)
    ir = 0.4; fr = 0.4; r0 = 0.4; correct = 0.0; RADIAL_CROSSING(fr,ir,r0,correct,got)

    printf("Testing finished with radial crossings.\n");


    return nFails;



}

int test_plane_crossings(){

    // Remember, we work in radians here!

    int nFails = 0;
    real correct,fa,ia,a0,got;

#define PLANE_CROSSING(FA,IA,A0,CORRECT,GOT) \
    got = diag_orb_check_plane_crossing( FA, IA, A0);\
    if ( fabs( CORRECT - GOT)  > TOLERANCE ) {\
        nFails += 1;\
        printf(" Plane: ini=%g r0=%g final=%g FAIL: got=%g, expected=%g\n    ---> FAIL!\n",IA,A0,FA,GOT,CORRECT);\
    }\
    else {\
    printf(" Plane: ini=%g r0=%g final=%g got=%g, expected=%g ---> OK \n",IA,A0,FA,GOT,CORRECT);\
    }

    printf("Testing plane crossings.\n");


    ia = 0.1; fa = 0.3; a0 = 0.2; correct = 0.5; PLANE_CROSSING(fa,ia,a0,correct,got)
    ia = 0.1; fa = 0.1; a0 = 0.2; correct = 0.0; PLANE_CROSSING(fa,ia,a0,correct,got)
    ia = CONST_PI-0.1; fa = CONST_PI+0.1; a0 = CONST_PI; correct = 0.5; PLANE_CROSSING(fa,ia,a0,correct,got)
    ia = CONST_2PI-0.1; fa = CONST_2PI+0.1; a0 = CONST_2PI; correct = 0.5; PLANE_CROSSING(fa,ia,a0,correct,got)
    ia = -0.1; fa =0.1; a0 = 0.0; correct = 0.5; PLANE_CROSSING(fa,ia,a0,correct,got)
    ia = 0.1+3*CONST_2PI; fa = 0.3+3*CONST_2PI; a0 = 0.2; correct = 0.5; PLANE_CROSSING(fa,ia,a0,correct,got)
    ia = 0.1+3*CONST_2PI; fa = 0.1+3*CONST_2PI; a0 = 0.2; correct = 0.0; PLANE_CROSSING(fa,ia,a0,correct,got)
    ia = 0.1-3*CONST_2PI; fa = 0.3-3*CONST_2PI; a0 = 0.2; correct = 0.5; PLANE_CROSSING(fa,ia,a0,correct,got)
    ia = 0.1-3*CONST_2PI; fa = 0.1-3*CONST_2PI; a0 = 0.2; correct = 0.0; PLANE_CROSSING(fa,ia,a0,correct,got)



    printf("Testing finished with plane crossings.\n");

    return nFails;



}


int test_update_gc(){

    int nFail = 0;

    int iSIMD;

//diag_orb_update_gc(diag_orb_data* data, particle_simd_gc* p_f,
//                        particle_simd_gc* p_i)

#define NFLD  DIAG_ORB_GCFIELDS
#define NMRK  NSIMD
#define NPNT  1
#define NPOLPLOTS 0
#define NTORPLOTS 0
#define NRADPLOTS 2
#define RADIALSTEP 0.025

#define nData   NFLD * NMRK * NPNT

    diag_orb_data data;
    diag_orb_offload_data offload_data;
    real offload_array[nData];
    real correct, got;
    particle_simd_gc p_f, p_i;

    printf("Testing update_gc in poincare mode for radial steps.\n");


    offload_data.mode              = DIAG_ORB_POINCARE;   /**< Defines condition for recording markers       */
    offload_data.Nfld              = NFLD;/**< Number of fields the record contains          */
    offload_data.Nmrk              = NMRK;/**< Number of markers to record                   */
    offload_data.Npnt              = NPNT;/**< Maximum number of points to keep recorded     */
    offload_data.record_mode       = simulate_mode_gc;    /**< Defines what fields are initialized           */
    offload_data.writeInterval     = 0.01; /**< Interval at which markers are recorded        */
    offload_data.ntoroidalplots    = NTORPLOTS; /**< Number of toroidal Poincare planes            */
    offload_data.npoloidalplots    = NPOLPLOTS;/**< Number of toroidal Poincare planes            */
    offload_data.nradialplots      = NRADPLOTS;   /**< Number of radial Poincare planes              */
    offload_data.toroidalangles[0] = 0.0; /**< Toroidal plane angles */
    offload_data.toroidalangles[1] = 0.1;
    offload_data.poloidalangles[0] = 0.0; /**< Poloidal plane angles */
    offload_data.poloidalangles[1] = 0.1;
    offload_data.radialdistances[0]= 0.11;    /**< Radial plane angles*/
    offload_data.radialdistances[1]= 0.21;

    // Initialize the data
    diag_orb_init(&data, &offload_data, offload_array);


    for ( iSIMD=0; iSIMD<NSIMD; iSIMD++ ){
        p_i.rho[iSIMD]    =  RADIALSTEP*iSIMD;    p_f.rho[iSIMD]     = RADIALSTEP*(iSIMD+1);
        //printf("rho %g -- %g \n",RADIALSTEP*iSIMD,RADIALSTEP*(iSIMD+1) );
        p_i.time[iSIMD]   =  RADIALSTEP*iSIMD;          p_f.time[iSIMD]    = RADIALSTEP*(iSIMD+1);
        p_i.phi[iSIMD]    = -0.05;         p_f.phi[iSIMD]     = 0.05;
        p_i.theta[iSIMD]  = -0.05;         p_f.theta[iSIMD]   = 0.05;
        p_i.mileage[iSIMD]=  RADIALSTEP*iSIMD;          p_f.mileage[iSIMD] = RADIALSTEP*(iSIMD+1);
        p_i.id[iSIMD]     =  iSIMD+20;     p_f.id[iSIMD]      =  iSIMD+20;
        p_i.index[iSIMD]  = iSIMD;         p_f.index[iSIMD]  = iSIMD;
    }


    diag_orb_update_gc(&data, &p_f, &p_i);


    real rhoi,rhof;
    integer imrk,ipoint,idx;

    for ( iSIMD=0; iSIMD<NSIMD; iSIMD++ ){
        imrk   = p_f.index[iSIMD];
        ipoint = data.mrk_pnt[imrk];
        idx    = imrk * data.Npnt + ipoint;
        rhoi = p_i.rho[iSIMD];
        rhof = p_f.rho[iSIMD];
        if(       rhoi < offload_data.radialdistances[0] &&  rhof > offload_data.radialdistances[0]    ) {
            correct = offload_data.radialdistances[0];
        }else if (rhoi < offload_data.radialdistances[1] &&  rhof > offload_data.radialdistances[1]   ) {
            correct = offload_data.radialdistances[1];
        }else{
            correct = 0.0;
        }

        got = data.rho[idx];
        if ( fabs( got-correct ) > TOLERANCE ){
            nFail++;
            printf(" Rho %6g->%6g expected=%g but got=%g.\n FAIL\n",rhoi,rhof,correct,got);
            printf("   time=%g\n",data.mileage[idx]);
            printf(" Index %ld\n",idx);
        }else{
            printf(" Rho %-6g->%-6g expected=%g and got=%g. OK\n",rhoi,rhof,correct,got);
        }




        if (nFail > 0) {
                /* We don't need to print errors again and again...*/
                //break;
        }
    }


    diag_orb_free(&data);


    printf("   The test does not cover the whole modality. Just a single quantity in one modality.\n");

    printf("Finished with update_gc.\n");

    return nFail;
}



int main(int argc, char **argv) {

    printf("\n\nTesting diag_orb.c\n\n");


    int fails=0;

    fails += test_radial_crossings();

    fails += test_plane_crossings();

    fails += test_update_gc();

    printf("\nFinished testing diag_orb.c\n\n");


    return fails;

}
