#!/usr/bin/env python
# -*- coding: utf-8 -*-
import math
import numpy as np
import matplotlib.pyplot  as plt
from analyticBKGpsifun import *

#analyticBfield
#Code copied and translated to python2 from https://solps-mdsplus.aug.ipp.mpg.de/wsvn/ASCOT/branches/simppa/matlab/analyticEquilibrium/
#Solution is based on this article:  http://dx.doi.org/10.1063/1.3328818

def analyticGS(epsilon,kappa,delta,Xpointx,Xpointy,D=None,R0=None):
    alpha = math.asin(delta) #alpha as defined in the article

    outerEqPoint = [1+epsilon, 0]
    innerEqPoint = [1-epsilon, 0]
    upperEqPoint = [1-epsilon*delta, kappa*epsilon]
    upperHighpoint = [1-epsilon*delta, kappa*epsilon]
    lowerXPoint = [xsep, ysep]
    outerSlope = 0 #outer equatorial point slope
    innerSlope = 0 #inner equatorial point slope
    curvatureOuterEq = -(1+alpha)**2/(epsilon*kappa**2) #curvature at the outboard midplane
    curvatureInnerEq = (1-alpha)**2/(epsilon*kappa**2) #curvature at the inboard midplane
    curvatureHighpointEq = -kappa/(epsilon*(math.cos(alpha))**2) #curvature at the top @onko alpha radiaaneissa?
    beta = -0.155 #beta is the parameter determining the beta regime of interest (Note: it's called A in the article)
    #beta=1  #: force-free equilibrium
    #beta=0 #: vacuum toroidal field ("low-beta" equilibrium)
    aGS_allInputs(epsilon,kappa,delta,Xpointx,Xpointy,outerEqPoint,innerEqPoint,upperEqPoint,upperHighpoint,lowerXPoint,outerSlope,innerSlope,curvatureOuterEq,curvatureInnerEq,curvatureHighpointEq,beta,D,R0)

def aGS_allInputs(epsilon,kappa,delta,Xpointx,Xpointy,outerEqPoint,innerEqPoint,upperEqPoint,upperHighpoint,lowerXPoint,outerSlope,innerSlope,curvatureOuterEq,curvatureInnerEq,curvatureHighpointEq,beta,D=None,R0=None):

    ###########################################################################
    #                                                                         #
    #   Construct the matrix A of the boundary conditions for the funtions    #
    #   which are solutions to the homogeneous equation                       #
    #                                                                         #
    ###########################################################################

    A = np.zeros((12,12))
    for ipsi in range(12):
        A[0][ipsi] = psi(outerEqPoint[0], outerEqPoint[1], ipsi)
        A[1][ipsi] = psi(innerEqPoint[0], innerEqPoint[1], ipsi)
        A[2][ipsi] = psi(upperEqPoint[0], upperEqPoint[1], ipsi)
        A[3][ipsi] = psi(lowerXPoint[0], lowerXPoint[1], ipsi)
        A[4][ipsi] = psix(outerEqPoint[0], outerEqPoint[1], ipsi)*outerSlope+psiy(outerEqPoint[0], outerEqPoint[1], ipsi)
        A[5][ipsi] = psix(innerEqPoint[0], innerEqPoint[1], ipsi)*innerSlope+psiy(innerEqPoint[0], innerEqPoint[1], ipsi)
        A[6][ipsi] = psix(upperHighpoint[0], upperHighpoint[1], ipsi)
        A[7][ipsi] = psix(lowerXPoint[0], lowerXPoint[1], ipsi)
        A[8][ipsi] = psiy(lowerXPoint[0], lowerXPoint[1], ipsi)
        A[9][ipsi] = psiyy(outerEqPoint[0], outerEqPoint[1], ipsi)+curvatureOuterEq*psix(outerEqPoint[0], outerEqPoint[1], ipsi)
        A[10][ipsi]= psiyy(innerEqPoint[0], innerEqPoint[1], ipsi)+curvatureInnerEq*psix(innerEqPoint[0], innerEqPoint[1], ipsi)
        A[11][ipsi]= psixx(upperEqPoint[0], upperEqPoint[1], ipsi)+curvatureHighpointEq*psiy(upperEqPoint[0], upperEqPoint[1], ipsi)

    ###########################################################################
    #                                                                         #
    #   Construct the matrix B of the boundary conditions for the particular  #
    #   solutions to the equation                                             #
    #                                                                         #
    ###########################################################################

    B = np.zeros((12,1))
    f = np.array([beta, 1-beta])
    for ipart in range(2):
        B[0] = B[0]+f[ipart]*psipart(outerEqPoint[0], outerEqPoint[1], ipart)
        B[1] = B[1]+f[ipart]*psipart(innerEqPoint[0], innerEqPoint[1], ipart)
        B[2] = B[2]+f[ipart]*psipart(upperEqPoint[0], upperEqPoint[1], ipart)
        B[3] = B[3]+f[ipart]*psipart(lowerXPoint[0], lowerXPoint, ipart)
        B[4] = B[4]+f[ipart]*(outerSlope*psipartx(outerEqPoint[0], outerEqPoint[1], ipart)+psiparty(outerEqPoint[0], outerEqPoint[1], ipart))
        B[5] = B[5]+f[ipart]*(outerSlope*psipartx(innerEqPoint[0], innerEqPoint[1], ipart)+psiparty(innerEqPoint[0], innerEqPoint[1], ipart))
        B[6] = B[6]+f[ipart]*psipartx(upperHighpoint[0], upperHighpoint[1], ipart)
        B[7] = B[7]+f[ipart]*psipartx(lowerXPoint[0], lowerXPoint[1], ipart)
        B[8] = B[8]+f[ipart]*psiparty(lowerXPoint[0], lowerXPoint[1], ipart)
        B[9] = B[9]+f[ipart]*(curvatureOuterEq*psipartx(outerEqPoint[0], outerEqPoint[1], ipart)+psipartyy(outerEqPoint[0], outerEqPoint[1], ipart))
        B[10] = B[10]+f[ipart]*(curvatureInnerEq*psipartx(innerEqPoint[0], innerEqPoint[1], ipart)+psipartyy(innerEqPoint[0], innerEqPoint[1], ipart))
        B[11] = B[11]+f[ipart]*(curvatureHighpointEq*psiparty(upperEqPoint[0], upperEqPoint[1], ipart)+psipartxx(upperEqPoint[0], upperEqPoint[1], ipart))

    B = np.negative(B)

    ###########################################################################
    #                                                                         #
    #   Solve the linear system for the coefficients C of the general         #
    #   solution to the equation                                              #
    #                                                                         #
    ###########################################################################


    C = np.linalg.lstsq(A,B)[0].flatten()

    #Write C-data in a text file:
    np.savetxt('analyticBKG.txt',C)

    ###########################################################################
    #                                                                         #
    #   Plot the solution if parameters r0 and D were given                   #
    #                                                                         #
    ###########################################################################

    if R0 is not None and D is not None:
        R= np.linspace((1-epsilon-0.1),1+epsilon+0.5,200) 
        Z= np.linspace(-1.2*kappa*epsilon,1.2*kappa*epsilon,200) 
        X,Y = np.meshgrid(R,Z) 
        X[X < 0.2] = 0
        Z = psi0(X,Y,C[0],C[1],C[2],C[3],C[4],C[5],C[6],C[7],C[8],C[9],C[10],C[11],beta) 
        plt.contour(X*R0, Y*R0, Z, D) 
        plt.plot([0, 0],[-1.5*kappa*epsilon, 1.5*kappa*epsilon], '--')
        plt.xlim(np.array([(1-epsilon-0.5), 1+epsilon+0.5])*R0)
        plt.ylim(np.array([-kappa*epsilon-0.2, kappa*epsilon+0.2])*R0)
        plt.show()

if __name__ == "__main__":
    R0=6.2
    method='itertok'

    epsilon=kappa=delta=xsep=ysep=D=0    
    if method == 'itertok':
        epsilon = 2/6.2 #Inverse aspect ratio
        kappa = 1.7 #Elongation
        delta  = 0.33 #Triangularity
        xsep = 1-1.1*delta*epsilon #x-location of the separatrix
        ysep = -1.1*kappa*epsilon #y-location of the separatrix
        D=[0,-0.005,-0.01,-0.015,-0.02,-0.025,-0.03,-0.035] #Values of the contours for plotting
    elif method == 'nstx':
        epsilon = 0.67/0.86 #Inverse aspect ratio
        kappa = 2 #Elongation
        delta  = 0.35 #Triangularity
        xsep = 1-1.1*delta*epsilon #x-location of the separatrix
        ysep = -1.1*kappa*epsilon #y-location of the separatrix
        D=[0,-0.015,-0.04,-0.075,-0.115,-0.16,-0.2,-0.23,-0.245] #Values of the contours for plotting
    elif method == 'aug':
        R0=1.65 
        epsilon = 0.5/R0 #Inverse aspect ratio
        kappa = 1.6 #Elongation
        delta  = 0.2 #Triangularity
        ysep = -0.9041 /R0 #x-location of the separatrix (x-point)
        xsep = 1.4528 /R0 #y-location of the separatrix
        D= -np.linspace(-0.005,0.04,15) #Values of the contours for plotting

    analyticGS(epsilon,kappa,delta,xsep,ysep,D,R0)
