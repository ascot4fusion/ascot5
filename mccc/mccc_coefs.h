/**
 * @file mccc_coefs.h
 * @brief Header file for mccc_coefs.c
*/
#ifndef MCCC_COEFS_H
#define MCCC_COEFS_H

#include "../ascot5.h"

void mccc_coefs_fo(real ma, real qa, real va, real* mb, real* qb, real* nb, real* Tb, real* clogab, int nspec, 
		   real* F, real* Dpara, real* Dperp, real* K, real* nu);

void mccc_coefs_gcfixed(real ma, real qa, real va, real xi, real* mb, real* qb, real* nb, real* Tb, real B, real* clogab, int nspec, 
			real* Dpara, real* Dx, real* K, real* nu);

void mccc_coefs_gcadaptive(real ma, real qa, real va, real xi, real* mb, real* qb, real* nb, real* Tb, real B, real* clogab, int nspec, 
			   real* Dpara, real* DX, real* K, real* nu, real* dQ, real* dDpara);

void mccc_coefs_clog(real ma, real qa, real va, real* mb, real* qb, real* nb, real* Tb, real* clogab, int nspec);

#endif
