#include <stdlib.h>
#include <math.h>
#include "../ascot5.h"
#include "math.h"
#include "../consts.h"

void mccc_special_G(real x, real* G, int exact){

    if(x == 0){
	G[0] = 0;
	G[1] = 0;
	return;
    }
	
    if(exact) {
	real expm2x = exp(-x*x);
	real erfx = erf(x);
	    
	G[0] = ( erfx - ( 2 * x / sqrt(CONST_PI) ) * expm2x )/ (x*x);
	G[1] = erfx - 0.5 * G[0];
    }
    else{
	//TODO implement
	G[0] = 0;
	G[1] = 0;
    }
	
}

void mccc_special_GdG(real x, real* GdG, int exact){

    if(x == 0){
	GdG[0] = 0;
	GdG[1] = 0;
	GdG[2] = 4/(3*sqrt(CONST_PI));
	GdG[3] = GdG[2];
	return;
    }
	
    if(exact) {
	real expm2x = exp(-x*x);
	real erfx = erf(x);
	    
	GdG[0] = ( erfx - ( 2 * x / sqrt(CONST_PI) ) * expm2x )/ (x*x);
	GdG[1] = erfx - 0.5 * GdG[0];

	GdG[2] = 4*expm2x/sqrt(CONST_PI) - 2*GdG[0]/x;
	GdG[3] = GdG[0]/x;
    }
    else{
	//TODO implement
	GdG[0] = 0;
	GdG[1] = 0;
	GdG[2] = 0;
	GdG[3] = 0;
    }
	
}

void mccc_special_fo(real x, real* fdf, int exact){
#ifdef MCCC_RELATIVISTIC

#else
    if(x == 0){
	fdf[0] = 0;
	fdf[1] = 0;
	fdf[2] = 4/(3*sqrt(CONST_PI));
	return;
    }
	
    if(exact) {
	real expm2x = exp(-x*x);
	real erfx = erf(x);
	    
	fdf[0] = ( erfx - ( 2 * x / sqrt(CONST_PI) ) * expm2x )/ (x*x);
	fdf[1] = erfx - 0.5 * fdf[0];

	fdf[2] = 4*expm2x/sqrt(CONST_PI) - 2*fdf[0]/x;
    }
    else{
	//TODO implement
	fdf[0] = 0;
	fdf[1] = 0;
        fdf[2] = 0;
    }
#endif
	
}


void mccc_special_mu(real u, real th, real* mu, int exact){
    //TODO implement
    mu[0] = 0;
    mu[1] = 0;
    mu[2] = 0;
    if(exact) {

    }
    else{
	
    }
}

void mccc_special_mudmu(real u, real th, real* mudmu, int exact){
    //TODO implement
    mudmu[0] = 0;
    mudmu[1] = 0;
    mudmu[2] = 0;
    mudmu[3] = 0;
    mudmu[4] = 0;
    mudmu[5] = 0;
    if(exact) {
	
    }
    else{
	
    }
}
