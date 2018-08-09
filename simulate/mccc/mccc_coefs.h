/**
 * @file mccc_coefs.h
 * @brief Header file for mccc_coefs.c
*/
#ifndef MCCC_COEFS_H
#define MCCC_COEFS_H

#define G_STEP 2.e-2
#define G_NSLOT 150

#include "../../ascot5.h"
#include "../../error.h"

#pragma omp declare target
void mccc_coefs_init(real* coldata);
    
#pragma omp declare simd
a5err mccc_coefs_fo(real ma, real qa, real va, real* mb, real* qb, real* nb, real* Tb, real* clogab, int nspec, real* coldata, 
		    real* F, real* Dpara, real* Dperp, real* K, real* nu);

#pragma omp declare simd
a5err mccc_coefs_gcfixed(real ma, real qa, real va, real xi, real* mb, real* qb, real* nb, real* Tb, real B, real* clogab, int nspec, real* coldata, 
			 real* Dpara, real* Dx, real* K, real* nu);

#pragma omp declare simd
a5err mccc_coefs_gcadaptive(real ma, real qa, real va, real xi, real* mb, real* qb, real* nb, real* Tb, real B, real* clogab, int nspec, real* coldata, 
			    real* Dpara, real* DX, real* K, real* nu, real* dQ, real* dDpara);

#pragma omp declare simd
a5err mccc_coefs_clog(real ma, real qa, real va, real* mb, real* qb, real* nb, real* Tb, real* clogab, int nspec);

#pragma omp end declare target

#endif
