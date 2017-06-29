/**
 * @author Konsta Sarkimaki konsta.sarkimaki@aalto.fi
 * @file mccc_coefs.c
 * @brief Tools to evaluate collision coefficients
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../../math.h"
#include "../../consts.h"
#include "../../ascot5.h"
#include "mccc_special.h"

int MCCC_COEFS_EXACT = 1;

/**
 * @brief Initializes lookup tables for more(?) efficient evaluation
 *
 * @todo Not implemented
 */
void mccc_coefs_init(){

}

/**
 * @brief Evaluates coefficients in particle picture
 * 
 * @param ma test particle mass [kg]
 * @param qa test particle charge [C]
 * @param va test particle velocity [m/s]
 * @param mb pointer to array storing background species mass [kg]
 * @param qb pointer to array storing background species charge [C]
 * @param nb pointer to array storing background species density [1/m^3]
 * @param Tb pointer to array storing background species temperature [J]
 * @param clogab pointer to array storing the evaluated Coulomb logarithms
 * @param nspec number of background species
 * @param F pointer to array storing the evaluated Fokker-Planck coefficient [m/s^2]
 * @param Dpara pointer to array the evaluated storing parallel diffusion coefficient [m^2/s^3] 
 * @param Dperp pointer to array the evaluated storing perpendicular diffusion coefficient [m²/s^3]
 * @param K pointer to array storing the evaluated friction coefficient [m/s^2]
 * @param nu pointer to array storing the evaluated pitch collision frequency [1/s]
 *
 * @todo Implement relativistic coefficients
 */
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

/**
 * @brief Evaluates coefficients in guiding center picture for fixed scheme
 * 
 * @param ma guiding center mass [kg]
 * @param qa guiding center charge [C]
 * @param va guiding center velocity [m/s]
 * @param va guiding center pitch
 * @param mb pointer to array storing background species mass [kg]
 * @param qb pointer to array storing background species charge [C]
 * @param nb pointer to array storing background species density [1/m^3]
 * @param Tb pointer to array storing background species temperature [J]
 * @param B magnetic field magnitude [T]
 * @param clogab pointer to array storing the evaluated Coulomb logarithms
 * @param nspec number of background species
 * @param Dpara pointer to array the evaluated storing parallel diffusion coefficient [m^2/s^3]
 * @param DX pointer to array storing the classical diffusion coefficient [m^2/s]
 * @param K pointer to array storing the evaluated friction coefficient [m/s^2]
 * @param nu pointer to array storing the evaluated pitch collision frequency [1/s]
 *
 * @todo Implement relativistic coefficients
 */
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

   
/**
 * @brief Evaluates coefficients in guiding center picture for adaptive scheme
 * 
 * @param ma guiding center mass [kg]
 * @param qa guiding center charge [C]
 * @param va guiding center velocity [m/s]
 * @param va guiding center pitch
 * @param mb pointer to array storing background species mass [kg]
 * @param qb pointer to array storing background species charge [C]
 * @param nb pointer to array storing background species density [1/m^3]
 * @param Tb pointer to array storing background species temperature [J]
 * @param B magnetic field magnitude [T]
 * @param clogab pointer to array storing the evaluated Coulomb logarithms
 * @param nspec number of background species
 * @param Dpara pointer to array the evaluated storing parallel diffusion coefficient [m^2/s^3]
 * @param DX pointer to array storing the classical diffusion coefficient [m^2/s]
 * @param K pointer to array storing the evaluated friction coefficient [m/s^2]
 * @param nu pointer to array storing the evaluated pitch collision frequency [1/s]
 * @param dQ pointer to array storing the evaluated derivative dQ/dv [1/s]
 * @param dDpara pointer to array storing the evaluated derivative dDpara/dv [m/s^2]
 *
 * @todo Implement relativistic coefficients
 */
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



/**
 * @brief Evaluate Coulomb logarithm
 *
 * @param ma guiding center mass [kg]
 * @param qa guiding center charge [C]
 * @param va guiding center velocity [m/s]
 * @param va guiding center pitch
 * @param mb pointer to array storing background species mass [kg]
 * @param qb pointer to array storing background species charge [C]
 * @param nb pointer to array storing background species density [1/m^3]
 * @param Tb pointer to array storing background species temperature [J]
 * @param clogab pointer to array storing the evaluated Coulomb logarithms
 * @param nspec number of background species
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


