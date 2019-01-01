/**
 * @file mccc_push.c
 * @brief Coulomb collision operators
 */
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "../../ascot5.h"
#include "../../math.h"
#include "../../consts.h"
#include "../../error.h"
#include "mccc_push.h"

/**
 * @brief Collision operator for particle with fixed time-step
 *
 * Computes the value for particle velocity after collisions with
 * background species using Euler-Maruyama method.
 *
 * @param F friction coefficient
 * @param Dpara parallel diffusion coefficient
 * @param Dperp perpendicular diffusion coefficient
 * @param dt time step
 * @param rnd three element array of normally distributed random numbers
 * @param vin pointer to particle velocity vector before collisions
 * @param vout pointer to particle velocity vector after collisions
 */
a5err mccc_push_foEM(real F, real Dpara, real Dperp, real dt, real* rnd, real* vin, real* vout){
    a5err err = 0;
    real dW[3],vhat[3];

    /* Wiener process for this step */
    real sdt = sqrt(dt);
    dW[0]=sdt*rnd[0];
    dW[1]=sdt*rnd[1];
    dW[2]=sdt*rnd[2];

    math_unit(vin,vhat);

    /* Use Euler-Maruyama method to get vout */
    real t1 = math_dot(vhat,dW);
    real k1 = F*dt;
    real k2 = sqrt(2*Dpara)*t1;
    real k3 = sqrt(2*Dperp);
    vout[0] = vin[0] + k1*vhat[0] + k2*vhat[0] + k3*(dW[0] - t1*vhat[0]);
    vout[1] = vin[1] + k1*vhat[1] + k2*vhat[1] + k3*(dW[1] - t1*vhat[1]);
    vout[2] = vin[2] + k1*vhat[2] + k2*vhat[2] + k3*(dW[2] - t1*vhat[2]);

    if( (vout[0]*vout[0] + vout[1]*vout[1] + vout[2]*vout[2]) >= CONST_C2 ) {err = error_raise(ERR_INTEGRATION, __LINE__, EF_MCCC_PUSH);}

    return err;
}

/**
 * @brief Collision operator for guiding center with fixed time-step
 *
 * Computes the value for guiding center position after collisions with
 * background species using Euler-Maruyama method.
 *
 * @param K friction coefficient
 * @param nu pitch collision frequency coefficient
 * @param Dpara parallel diffusion coefficient
 * @param DX classical diffusion coefficient
 * @param B pointer to local magnetic field vector
 * @param dt time step
 * @param rnd three element array of normally distributed random numbers
 * @param vin pointer to guiding center velocity vector before collisions
 * @param vout pointer to guiding center velocity vector after collisions
 * @param xiin pointer to guiding center pitch before collisions
 * @param xiout pointer to guiding center pitch after collisions
 * @param Xin pointer to Cartesian guiding center position before collisions
 * @param Xout pointer to Cartesian guiding center position after collisions
 * @param cutoff value below which particle velocity is mirrored
 */
a5err mccc_push_gcEM(real K, real nu, real Dpara, real DX, real* B, real dt, real* rnd, 
                     real vin, real* vout, real xiin, real* xiout, real* Xin, real* Xout, real cutoff){
    a5err err = 0;
    real dW[5], bhat[3];

    /* Wiener process for this step */
    real sdt = sqrt(dt);
    dW[0]=sdt*rnd[0]; // For X_1
    dW[1]=sdt*rnd[1]; // For X_2
    dW[2]=sdt*rnd[2]; // For X_3
    dW[3]=sdt*rnd[3]; // For v
    dW[4]=sdt*rnd[4]; // For xi

    math_unit(B,bhat);

    /* Use Euler-Maruyama method to get Xout, vout, and xiout */

    real k1 = sqrt(2*DX);
    real k2 = math_dot(bhat,dW);

    Xout[0] = Xin[0] + k1*(dW[0]-k2*bhat[0]);
    Xout[1] = Xin[1] + k1*(dW[1]-k2*bhat[1]);
    Xout[2] = Xin[2] + k1*(dW[2]-k2*bhat[2]);

    vout[0] = vin + K*dt + sqrt(2*Dpara)*dW[3];

    xiout[0] = xiin - xiin*nu*dt + sqrt((1 - xiin*xiin)*nu)*dW[4];

    /* Enforce boundary conditions */
    if(vout[0] < cutoff){
        vout[0] = 2*cutoff-(vout[0]);
    }

    if(fabs(xiout[0]) > 1){
        xiout[0] = ((xiout[0] > 0) - (xiout[0] < 0))*(2-fabs(xiout[0]));
    }

    if(!err && ( vout[0] < 0 || vout[0] >= CONST_C) ) {err = error_raise(ERR_INTEGRATION, __LINE__, EF_MCCC_PUSH);}
    if(!err && fabs(xiout[0]) >= 1)              {err = error_raise(ERR_INTEGRATION, __LINE__, EF_MCCC_PUSH);}

    return err;
}

/**
 * @brief Collision operator for guiding center with adaptive time-step
 *
 * Computes the value for guiding center position after collisions with
 * background species using Milstein method. Returns error estimates.
 *
 * @param K friction coefficient
 * @param nu pitch collision frequency coefficient
 * @param Dpara parallel diffusion coefficient
 * @param DX classical diffusion coefficient
 * @param B pointer to local magnetic field vector
 * @param dt time step
 * @param dW pointer to 5D Wiener process increments
 * @param dQ derivate of coefficient Q with respect to velocity
 * @param dDpara derivate of coefficient Dpara with respect to velocity
 * @param vin pointer to guiding center velocity vector before collisions
 * @param vout pointer to guiding center velocity vector after collisions
 * @param xiin pointer to guiding center pitch before collisions
 * @param xiout pointer to guiding center pitch after collisions
 * @param Xin pointer to Cartesian guiding center position before collisions
 * @param Xout pointer to Cartesian guiding center position after collisions
 * @param cutoff value below which particle velocity is mirrored
 * @param tol error tolerance
 * @param kappa_k pointer to drift part of the error estimate
 * @param kappa_d pointer to diffusion part of the error estimate
 */
a5err mccc_push_gcMI(real K, real nu, real Dpara, real DX, real* B, real dt, real* dW, real dQ, real dDpara, real vin, real* vout,
                     real xiin, real* xiout, real* Xin, real* Xout, real cutoff, real tol, real* kappa_k, real* kappa_d0, real* kappa_d1){
    a5err err = 0;
    real bhat[3];

    math_unit(B,bhat);

    /* Use Euler-Maruyama method to get Xout, vout, and xiout */

    real k1 = sqrt(2*DX);
    real k2 = math_dot(bhat,dW);

    Xout[0] = Xin[0] + k1*(dW[0]-k2*bhat[0]);
    Xout[1] = Xin[1] + k1*(dW[1]-k2*bhat[1]);
    Xout[2] = Xin[2] + k1*(dW[2]-k2*bhat[2]);

    vout[0]  = vin + K*dt + sqrt(2*Dpara)*dW[3] + 0.5*dDpara*(dW[3]*dW[3] - dt);
    xiout[0] = xiin - xiin*nu*dt + sqrt((1 - xiin*xiin)*nu)*dW[4] - 0.5*xiin*nu*(dW[4]*dW[4] - dt);

    /* Enforce boundary conditions */
    if(vout[0] < cutoff){
        vout[0] = 2*cutoff - vout[0];
    }

    if(fabs(xiout[0]) > 1){
        xiout[0] = ((xiout[0] > 0) - (xiout[0] < 0))*(2 - fabs(xiout[0]));
    }

    /* Error estimates for drift and diffusion limits */
    real erru = tol*(vin + fabs(K)*dt + sqrt(2*Dpara*dt)) + DBL_EPSILON;

    k1 = (1/(2*erru))*fabs(K*dQ);
    k2 = (1/(2*tol))*fabs(xiin*nu*nu);

    if(k1 > k2){
        kappa_k[0] = k1*dt*dt;
    }
    else{
        kappa_k[0] = k2*dt*dt;
    }

    kappa_d0[0] = (1/(6*erru))*fabs(dW[3]*dW[3]*dW[3]*dDpara*dDpara/sqrt(Dpara));
    kappa_d1[0] = sqrt(1-xiin*xiin)*nu*sqrt(nu)*fabs(dW[4] + sqrt(dt/3))*dt/(2*tol);

    if(!err && vout[0] < 0)                      {err = error_raise(ERR_INTEGRATION, __LINE__, EF_MCCC_PUSH);}
    if(!err && fabs(xiout[0]) >= 1)              {err = error_raise(ERR_INTEGRATION, __LINE__, EF_MCCC_PUSH);}

    return err;
}
