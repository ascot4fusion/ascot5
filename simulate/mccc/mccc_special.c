/**
 * @author Konsta Sarkimaki konsta.sarkimaki@aalto.fi
 * @file mccc_special.c
 * @brief Special mathematical functions needed in mccc package
 * @todo Relativistic coefficients and table lookup
 */
#include <stdlib.h>
#include <math.h>
#include "../../ascot5.h"
#include "../../math.h"
#include "../../consts.h"
#include "mccc_special.h"
#include "mccc_coefs.h"

/**
 * @brief Special functions for non-relativistic coefficients
 */
void mccc_special_fo(real x, real* fdf, int exact, real* coldata){
    
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

	if( x > G_STEP*(G_NSLOT-1) ) {
	    fdf[0] = 1/(2*x*x);
	    fdf[1] = 1.0 - 1/(2*x*x);
	    fdf[2] = -1/(x*x*x);
	}
	else {
	    int i = (int)floor(x/G_STEP);
	    fdf[0] = coldata[i] + (x-(i-1)*G_STEP)*( coldata[i+1] - coldata[i] )/ G_STEP;

	    i = i+G_NSLOT;
	    fdf[1] = coldata[i] + (x-(i-1)*G_STEP)*( coldata[i+1] - coldata[i] )/ G_STEP;

	    i = i+G_NSLOT;
	    fdf[2] = coldata[i] + (x-(i-1)*G_STEP)*( coldata[i+1] - coldata[i] )/ G_STEP;
	}
    }
	
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
