/** MCCC operator for guiding center coordinates.
* These coordinates are location X, momentum p (scalar), and
* pitch = p_para/p. The collisions can be evaluated separately
* for each coordinate. Both fixed and adaptive scheme are supported,
* and Milstein method is used for both.
*
*	@author Konsta Sarkimaki konsta.sarkimaki@aalto.fi
*/
#include <stdio.h>
#include <math.h>
#include "../ascot5.h"
#include "../math.h"
#include "mccc_push.h"

/** Computes the value for particle momentum after collisions with
* background species using Euler-Maruyama method with a fixed time step.
* 
* @param mccc_special dat  -- mccc_special struct that is obtained with the mccc_init() call
* @param  ma        -- test particle mass [kg]
* @param  qa        -- test particle charge [C]
* @param  clogab(:) -- list of coulomb logarithms for species a colliding with b [1] 
* @param  mb(:)     -- list of background species masses [kg]
* @param  qb(:)     -- list of background species charges [C]
* @param  nb(:)     -- list of background densities [1/m^3]
* @param  thb(:)    -- list of normalized background temperatures thb=T_b/(m_b*c^2) [1] 
* @param  dt        -- time step length [s]
* @param  rnd(3)    -- array with three elements of standard normal random numbers ~ N(0,1)
* @param  uin(3)    -- normalized test particle momentum u=p/mc before collisions [1] 
* @param  uout(3)   -- normalized test particle momentum u=p/mc after collisions [1]
* @param err       -- error flag, negative indicates something went wrong
*/
void mccc_push_foEM(real F, real Dpara, real Dperp, real dt, real* rnd, real* vin, real* vout, int* err){
    
    real dW[3],vhat[3];
    *err = 0;

    // Wiener process for this step
    dW[0]=sqrt(dt)*rnd[0];
    dW[1]=sqrt(dt)*rnd[1];
    dW[2]=sqrt(dt)*rnd[2];
    
    math_unit(vin,vhat);

    // Use Euler-Maruyama method to get vout
    real t1 = math_dot(vhat,dW);
    real k1 = F*dt;
    real k2 = sqrt(2*Dpara)*t1;
    real k3 = sqrt(2*Dperp);
    vout[0] = vin[0]+k1*vhat[0]+k2*vhat[0]+k3*(dW[0]-t1*vhat[0]);
    vout[1] = vin[1]+k1*vhat[1]+k2*vhat[1]+k3*(dW[1]-t1*vhat[1]);
    vout[2] = vin[2]+k1*vhat[2]+k2*vhat[2]+k3*(dW[2]-t1*vhat[2]);
	
    if(isnan(vout[0]) || isnan(vout[0]) || isnan(vout[0])){ *err = MCCC_PUSH_ISNAN;}
    if(isinf(vout[0]) || isinf(vout[0]) || isinf(vout[0])){ *err = MCCC_PUSH_ISNAN;}

}

void mccc_push_gcEM(real K, real nu, real Dpara, real DX, real* B, real dt, real* rnd, 
		    real vin, real* vout, real xiin, real* xiout, real* Xin, real* Xout, real cutoff, int* err){

    real dW[5], bhat[3];
    *err = 0;

    // Wiener process for this step
    dW[0]=sqrt(dt)*rnd[0]; // For X_1
    dW[1]=sqrt(dt)*rnd[1]; // For X_2
    dW[2]=sqrt(dt)*rnd[2]; // For X_3
    dW[3]=sqrt(dt)*rnd[3]; // For v
    dW[4]=sqrt(dt)*rnd[4]; // For xi

    math_unit(B,bhat);

    // Use Euler-Maruyama method to get Xout, vout, and xiout

    real k1 = sqrt(DX);
    real k2 = math_dot(bhat,dW);

    Xout[0] = Xin[0] + k1*(dW[0]-k2*bhat[0]);
    Xout[1] = Xin[1] + k1*(dW[1]-k2*bhat[1]);
    Xout[2] = Xin[2] + k1*(dW[2]-k2*bhat[2]);

    *vout = vin + K*dt + sqrt(2*Dpara)*dW[3];

    *xiout = xiin - xiin*nu*dt + sqrt((1-xiin*xiin)*nu)*dW[4];

    // Enforce boundary conditions
    if(*vout < cutoff){
	*vout = 2*cutoff-(*vout);
    }

    if(fabs(*xiout) > 1){
	*xiout = ((*xiout > 0) - (*xiout < 0))*(2-fabs(*xiout));
    }
    

    if(isnan(Xout[0]) || isnan(Xout[0]) || isnan(Xout[0]) || isnan(*vout) || isnan(*xiout)){ *err = MCCC_PUSH_ISNAN;}
    if(isinf(Xout[0]) || isinf(Xout[0]) || isinf(Xout[0]) || isnan(*vout) || isnan(*xiout)){ *err = MCCC_PUSH_ISNAN;}

    if(*vout <= 0 || fabs(*xiout) > 1){ *err = MCCC_PUSH_NOTPHYSICAL;}
}

void mccc_push_gcMI(real K, real nu, real Dpara, real DX, real* B, real dt, real* dW, real dQ, real dDpara, real vin, real* vout, 
		    real xiin, real* xiout, real* Xin, real* Xout, real cutoff, real tol, real* kappa_k, real* kappa_d, int* err){

    int rejected = 0;
    real bhat[3];
    *err = 0;

    math_unit(B,bhat);

    // Use Euler-Maruyama method to get Xout, vout, and xiout

    real k1 = sqrt(DX);
    real k2 = math_dot(bhat,dW);

    Xout[0] = Xin[0] + k1*(dW[0]-k2*bhat[0]);
    Xout[1] = Xin[1] + k1*(dW[1]-k2*bhat[1]);
    Xout[2] = Xin[2] + k1*(dW[2]-k2*bhat[2]);

    *vout = vin + K*dt + sqrt(2*Dpara)*dW[3] + 0.5*dDpara*(dW[3]*dW[3]-dt);
    *xiout = xiin - xiin*nu*dt + sqrt((1-xiin*xiin)*nu)*dW[4] - 0.5*xiin*nu*(dW[4]*dW[4]-dt);

    // Enforce boundary conditions
    if(*vout < cutoff){
	*vout = 2*cutoff-(*vout);
    }

    if(fabs(*xiout) > 1){
	*xiout = ((*xiout > 0) - (*xiout < 0))*(2-fabs(*xiout));
    }
    
    // Error estimates for drift and diffusion limits
    real erru = tol*(fabs(K)*dt + sqrt(2*Dpara));

    k1 = (1/(2*erru))*fabs(K*dQ);
    k2 = (1/(2*tol))*fabs(xiin*nu*nu);

    if(k1 > k2){
	*kappa_k = k1*dt*dt;
    }
    else{
	*kappa_k = k2*dt*dt;
    }

    kappa_d[0] = (1/(6*erru))*fabs(dW[3]*dW[3]*dW[3]*dDpara*dDpara/sqrt(Dpara));
    kappa_d[1] = sqrt(1-xiin*xiin)*nu*sqrt(nu)*fabs(dW[4] + sqrt(dt/3))*dt/(2*tol);

    if(isnan(Xout[0]) || isnan(Xout[0]) || isnan(Xout[0]) || isnan(*vout) || isnan(*xiout)){ *err = MCCC_PUSH_ISNAN;}
    if(isinf(Xout[0]) || isinf(Xout[0]) || isinf(Xout[0]) || isnan(*vout) || isnan(*xiout)){ *err = MCCC_PUSH_ISNAN;}

    if(*vout <= 0 || fabs(*xiout) > 1){ *err = MCCC_PUSH_NOTPHYSICAL;}

}

/** Computes the value for guiding center phase space position after collisions with
* background species using Milstein method with a fixed or an adaptive time step.
* Returns a suggestion for the next time step. A negative suggestion indicates 
* the current time step was rejected, and absolute value should be used for the 
* next try. If no tolerance is given, Milstein method with fixed time step will
* be used instead. Wiener processes are generated automatically and only cleaning
* is required after each accepted time step.
*
* @param mcccgc_coefstruct coefs -- struct containing the coefficients from mccca_evalCoefs() call
* @param w     -- array that stores the Wiener processes (Ndim=5)
* @param B             -- magnetic field at guiding center position before collisions [T]
* @param cutoff           -- value (>= 0) used to mirror u as "u = 2*cutoff - u" if u < cutoff otherwise [1]
* @param uin              -- normalized test particle momentum u=p/mc before collisions [1]
* @param xiin             -- guiding center pitch before collisions [1]
* @param dtin             -- candidate time step length [s]
* @param tol           -- relative tolerance for the guiding center momentum and absolute tolerance for the pitch. 
*                            Fixed time step is used if negative
* @param uout  -- normalized guiding center momentum u=p/mc after collisions [1]
* @param xiout -- guiding center pitch after collisions [1]
* @param dx -- change in guiding center spatial position due to collisions [m]
* @param dtout -- suggestion for the next time step. dtout = dtin in case fixed step is used [s]
* @param err  -- error flag, negative value indicates error
*/	
/**							
void mccc_push_gcM(coefs,mccc_wienarr* w, real* B, real cutoff, 
					real uin, real* uout, real xiin, real* xiout,
					real* dx, real dtin, real* dtout, real tol, int* err){
						
    real*8 :: bhat(3)
    real*8 :: F, gu, gxi
    real*8 :: kappa_k, kappa_d(2), erru, gx
    real*8 :: dWopt(2),dti,alpha
    int i, windex, tindex, ki, kmax;
    logical rejected
    
	int Ndim = 5;

    // Generate or retrieve a Wiener process for this step (t = time+dtin)
    real time w->time[0];
    wiener_generate(w, time+dtin, &windex, &err);
    if{err < 0} return;
	
	// change in the Wiener process during dt
	real dW[Ndim], dW2[Ndim];
	for{i = 0; i < Ndim; 1}{
		dW[i] = w->wienarr[i + windex*Ndim] - w->wienarr[i];
		dW2[i] = dW[i]*dW[i] - dtin;
	}
    windex = 0;

	bhat = math_unit(B);
    gx = sqrt(2*abs(coefs%gx/dot_product(B,B)));
    F = coefs%kappa+coefs%dDpar+2*coefs%Dpar/uin
    gu = sqrt(2*coefs%Dpar)
    gxi = sqrt((1-xiin**2)*coefs%nu)
    
    if{tol <= 0} {
		real t1 = dot_product(dW,bhat);
		dx[0] = gx * dW[0] - t1[0]*bhat[0];
		dx[1] = gx * dW[1] - t1[1]*bhat[1];
		dx[2] = gx * dW[2] - t1[2]*bhat[2];
		uout = uin + F*dtin + gu*dW[3] + coefs%dDpar*dW2[3]/2;
		xiout = xiin - xiin*coefs%nu*dtin + gxi*dW[4] - xiin*coefs%nu*dW2[4]/2;

       // Enforce boundary conditions
       if{uout < cutoff}{
          uout = 2*cutoff-uout;
       }

       if{abs(xiout) > 1.0}{
          xiout = sign(2-abs(xiout),xiout);
       }

       time = time + dtin;
       dtout = dtin;
       wiener_generate(w, time, windex, err);
       if{err < 0} return;
       if{isnan(uout) || isnan(xiout) || all(isnan(dx)) || isnan(dtout)} err = -1;
       if{(uout < 0) || (abs(xiout) > 1)} err = -1;
       return
    }

    // Error estimates for drift and diffusion limits
    real erru = tol*(abs(F)*dtin + gu*sqrt(dtin));
    real kappa_k = max((1.0/(2*erru))*abs(coefs%kappa*coefs%dkappa),&
         (1.0/(2*tol))*abs(xiin*coefs%nu**2))*dtin**2
    real kappa_d = max((1.0/(6*erru))*abs(dW[3]*dW[3]*dW[3])*abs(coefs%dDpar)*abs(coefs%dDpar/sqrt(coefs%Dpar)),&
         (1.0/(1.15**2*tol))*dtin*(abs(dW[4]+sqrt(dtin/3)))*abs(coefs%nu*gxi))

    // If the time step is accepted, use Milstein method to calculate new momentum
    rejected = .true.;
    if{(kappa_k < 1) && all(kappa_d < 1)} {
       real t1 = dot_product(dW,bhat);
		dx[0] = gx * dW[0] - t1[0]*bhat[0];
		dx[1] = gx * dW[1] - t1[1]*bhat[1];
		dx[2] = gx * dW[2] - t1[2]*bhat[2];
		uout = uin + F*dtin + gu*dW[3] + coefs%dDpar*dW2[3]/2;
		xiout = xiin - xiin*coefs%nu*dtin + gxi*dW[4] - xiin*coefs%nu*dW2[4]/2;

       // Enforce boundary conditions
       if{uout < cutoff}{
          uout = 2*cutoff-uout;
       }

       if{abs(xiout) > 1.0}{
          xiout = sign(2-abs(xiout),xiout);
       }

       time = time + dtin;
       dtout = dtin;
       wiener_generate(w, time, windex, err);
       if{err < 0} return;
       if{isnan(uout) || isnan(xiout) || all(isnan(dx)) || isnan(dtout)} err = -1;
       if{(uout < 0) || (abs(xiout) > 1)} err = -1;
       
       rejected = .false.
    }
     
    // Different time step estimates are used depending which error estimate dominates
    // This scheme automatically takes care of time step reduction (increase) when time step is rejected (accepted)
	real dWopt[2];
	real k = pow(kappa_d,-1.0/3);
    dWopt[0] = 0.9*abs(dW[3])*k;
	dWopt[1] = 0.9*abs(dW[4])*k;
    real alpha = MAX(abs(dW[3]), abs(dW[4]))/sqrt(dtin);
    if{kappa_k .gt. maxval(kappa_d)} {
       dti = min(1.5*dtin,0.8*dtin/sqrt(kappa_k))
       for{ki=1,ki < 3,1}{
          mccc_wiener_generate(w, time+ki*dti/3, tindex, err)
          if{err < 0} return
		  dW[3] = abs(w->wienarr[3 + windex*Ndim] - w->wienarr[3 + tindex*Ndim]);
		  if(dW[3] > dWopt[0]){exit;}
          dW[4] = abs(w->wienarr[4 + windex*Ndim] - w->wienarr[4 + tindex*Ndim]);
		  if(dW[4] > dWopt[1]){exit;}
       }
       dtout = MAX(MIN(ki-1,3),1)*(dti/3);
	}
    else{
       int kmax = 6;
       if (rejected) {
          kmax = 2;
	   }
       else if (alpha > 2) {
          kmax = 4;
       }

       for{ki=1,ki < kmax,1}{
          mccc_wiener_generate(w, time+ki*dtin/3, tindex, err);
		  if(err < 0){return;}
		  dW[3] = abs(w->wienarr[3 + windex*Ndim] - w->wienarr[3 + tindex*Ndim]);
		  if(dW[3] > dWopt[0]){exit;}
          dW[4] = abs(w->wienarr[4 + windex*Ndim] - w->wienarr[4 + tindex*Ndim]);
		  if(dW[4] > dWopt[1]){exit;}
       }
       dtout = MAX(MIN(ki-1,kmax),1)*(dtin/3);
    }

    // Rejected step is indicated by a negative time step suggestion
    if(rejected){
       dtout = -dtout;
    }

    if(isnan(uout) || isnan(xiout) || all(isnan(dx)) || isnan(dtout)){err = -1;}
	if(abs(dtout) < MCCCGC_EXTREMELYSMALLDTVAL){ err = -1;}
    if((uout < 0) || (abs(xiout) > 1)){ err = -1;}
						
}
*/
