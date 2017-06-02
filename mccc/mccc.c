#include <stdio.h>
#include "../ascot5.h"
#include "../.h"
#include "mccc.h"
#include "mccc_wiener.h"
#include "mccc_push.h"
#include "mccc_coefs.h"

void mccc_init(){

}

void mccc_update_fo(particle_simd_fo* fo){

}

void mccc_update_gc(particle_simd_fo* fo){

}

void mccc_step_fo_fixed(particle_simd_fo* p, real* dtin){

}

void mccc_step_gc_fixed(particle_simd_gc* p, real* dtin){

}

void mccc_step_gc_adaptive(particle_simd_gc* p, real* dtin, real* dtout){

}

/** Prints error description if err != 0.
 *
 * input:
 *
 * int err -- error flag
 */
void mccc_printerror(int err){

    if(err == 0){
	return;
    }
    else if(err == MCCC_WIENER_EXCEEDEDCAPACITY){
	printf("Error: Number of slots in Wiener array exceeded.\n");
    }
    else if(err == MCCC_WIENER_NOASSOCIATEDPROCESS){
	printf("Error: No associated process found.\n");
    }
    else if(err == MCCC_PUSH_ISNAN){
	printf("Error: Collision operator yields NaN or Inf.\n");
    }
    else{
	printf("Error: Unknown error\n");
    }    


}
