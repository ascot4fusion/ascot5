/**
* This module contains tools to evaluate Fokker-Planck coefficients
* from plasma parameters.
*
* @author Konsta Sarkimaki konsta.sarkimaki@aalto.fi
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../math.h"
#include "../consts.h"
#include "../ascot5.h"
#include "mccc_special.h"

int MCCC_COEFS_EXACT = 1;

void mccc_coefs_init(){

}

void mccc_coefs_fo(real ma, real qa, real va, real* mb, real* qb, real* nb, real* Tb, real* clogab, int nspec, 
		   real* F, real* Dpara, real* Dperp, real* K, real* nu){

    int i;
    real Q, dDpara;
    real x, vth, cab;

#ifdef MCCC_RELATIVISTIC

#else
    real fdf[3];

    for(i = 0; i < nspec; i=i+1){
	cab = nb[i]*qa*qa*qb[i]*qb[i]*clogab[i]/(4*CONST_PI*CONST_E0*CONST_E0);
	vth = sqrt(2*Tb[i]/mb[i]);
	x = va/vth;
	mccc_special_fo(x, fdf, MCCC_COEFS_EXACT);

	Q = -cab*fdf[0]/(ma*mb[i]*vth*vth);
	dDpara = (cab/(2*ma*ma*va)) * (fdf[2]/vth - fdf[0]/va);
	
	if(va == 0){
	    F[i] = 0;
	    Dpara[i] = (cab/(2*ma*ma))*4/(3*sqrt(CONST_PI)*vth);
	    Dperp[i] = (cab/(2*ma*ma))*4/(3*sqrt(CONST_PI)*vth);
	}
	else{
	    F[i] = (1+mb[i]/ma)*Q;
	    Dpara[i] = cab*fdf[0]/(2*ma*ma*va);
	    Dperp[i] = cab*fdf[1]/(2*ma*ma*va);

	}

	K[i] = Q + dDpara + 2*Dpara[i]/va;
	nu[i] = 2*Dperp[i]/(va*va);
    }

#endif
	
}

void mccc_coefs_gcfixed(real ma, real qa, real va, real xi, real* mb, real* qb, real* nb, real* Tb, real B, real* clogab, int nspec, 
			real* Dpara, real* DX, real* K, real* nu){

    int i;
    real Q, dDpara, Dperp;
    real x, vth, cab, gyrofreq;

#ifdef MCCC_RELATIVISTIC

#else
    real fdf[3];

    gyrofreq = qa*B/ma;

    for(i = 0; i < nspec; i=i+1){
	cab = nb[i]*qa*qa*qb[i]*qb[i]*clogab[i]/(4*CONST_PI*CONST_E0*CONST_E0);
	vth = sqrt(2*Tb[i]/mb[i]);
	x = va/vth;
	mccc_special_fo(x, fdf, MCCC_COEFS_EXACT);

	Q = -cab*fdf[0]/(ma*mb[i]*vth*vth);
	dDpara = (cab/(2*ma*ma*va)) * (fdf[2]/vth - fdf[0]/va);
	
	
	if(va == 0){
	    Dpara[i] = (cab/(2*ma*ma))*4/(3*sqrt(CONST_PI)*vth);
	    Dperp = (cab/(2*ma*ma))*4/(3*sqrt(CONST_PI)*vth);
	}
	else{
	    Dpara[i] = cab*fdf[0]/(2*ma*ma*va);
	    Dperp = cab*fdf[1]/(2*ma*ma*va);

	}

	K[i] = Q + dDpara + 2*Dpara[i]/va;
	nu[i] = 2*Dperp/(va*va);
	DX[i] = ( (Dpara[i] - Dperp)*(1-xi*xi)/2 + Dperp )/(gyrofreq*gyrofreq);
    }

#endif
	
}
   
void mccc_coefs_gcadaptive(real ma, real qa, real va, real xi, real* mb, real* qb, real* nb, real* Tb, real B, real* clogab, int nspec, 
			   real* Dpara, real* DX, real* K, real* nu, real* dQ, real* dDpara){
    int i;
    real Q, Dperp;
    real x, vth, cab, gyrofreq;

#ifdef MCCC_RELATIVISTIC

#else
    real fdf[3];

    gyrofreq = qa*B/ma;

    for(i = 0; i < nspec; i=i+1){
	cab = nb[i]*qa*qa*qb[i]*qb[i]*clogab[i]/(4*CONST_PI*CONST_E0*CONST_E0);
	vth = sqrt(2*Tb[i]/mb[i]);
	x = va/vth;
	mccc_special_fo(x, fdf, MCCC_COEFS_EXACT);

	Q = -cab*fdf[0]/(ma*mb[i]*vth*vth);
	dQ[i] = -cab*fdf[2]/(ma*mb[i]*vth*vth);
	dDpara[i] = (cab/(2*ma*ma*va)) * (fdf[2]/vth - fdf[0]/va);
	
	
	if(va == 0){
	    Dpara[i] = (cab/(2*ma*ma))*4/(3*sqrt(CONST_PI)*vth);
	    Dperp = (cab/(2*ma*ma))*4/(3*sqrt(CONST_PI)*vth);
	}
	else{
	    Dpara[i] = cab*fdf[0]/(2*ma*ma*va);
	    Dperp = cab*fdf[1]/(2*ma*ma*va);

	}

	K[i] = Q + dDpara[i] + 2*Dpara[i]/va;
	nu[i] = 2*Dperp/(va*va);
	DX[i] = ( (Dpara[i] - Dperp)*(1-xi*xi)/2 + Dperp )/(gyrofreq*gyrofreq);
    }

#endif
    
}



/**	Evalutes the Coulomb logarithm
*/
void mccc_coefs_clog(real ma, real qa, real va, real* mb, real* qb, real* nb, real* Tb, real* clogab, int nspec){

    int i;

    real u = (va/CONST_C)/sqrt(1 - (va*va)/CONST_C2 );
    real vbar[MAX_SPECIES];
    real s = 0;
    for(i = 0; i < nspec; i=i+1){
	s = s + (nb[i] * qb[i] * qb[i])/(Tb[i]);
	vbar[i] = va*va + 2*Tb[i]/mb[i];
    }
    real debyeLength = sqrt(CONST_E0/s);

    real mr, bcl, bqm;
    for(i=0; i < nspec; i=i+1){
	mr = ma*mb[i]/(ma+mb[i]);
	bcl = fabs(qa*qb[i]/(4 * CONST_PI * CONST_E0 * mr * vbar[i]));
	bqm = fabs(CONST_HBAR/(2*mr*sqrt(vbar[i])));
       
	if(bcl > bqm){
	    clogab[i] = log(debyeLength/bcl);
	}
	else{
	    clogab[i] = log(debyeLength/bqm);
	}
    }
}


