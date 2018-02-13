# -*- coding: utf-8 -*-
import matplotlib as plt
import numpy as np
import scipy.interpolate as sci

FONTTITLE = 12
FONTLABEL = 10
FONTTICK  = 8

def plotsurfandcontour(R, z, psi, title, ax=None, fig=None):
    if ax == None:
        ax = plt.figure()  
        fig = ax
    
    x,y = np.meshgrid(z,R)   
    
    surf = ax.pcolor(R,z,psi)
    cont = ax.contour(R,z,psi)
    
    ax.set_title(title, fontsize=FONTTITLE)
    ax.set_xlabel('R (m)', fontsize=FONTLABEL)
    ax.set_ylabel('z (m)', fontsize=FONTLABEL)
    ax.axis('image')
    plt.rc('font', size=FONTTICK)
    
    cbar = fig.colorbar(surf,ax=ax)
        
    return ax
    
def plotsurf(R, z, psi, title, ax=None, fig=None):
    if ax == None:
        ax = plt.figure()
        fig = ax
    
    x,y = np.meshgrid(z,R)   
    
    surf = ax.pcolor(R,z,psi)
    
    ax.set_title(title, fontsize=FONTTITLE)
    ax.set_xlabel('R (m)', fontsize=FONTLABEL)
    ax.set_ylabel('z (m)', fontsize=FONTLABEL)
    ax.axis('image')
    plt.rc('font', size=FONTTICK)
    
    cbar = fig.colorbar(surf,ax=ax)
        
    return ax
    
def plotripplemap(R, z, Bphi, ax=None, fig=None):
    if ax == None:
        ax = plt.figure()
        fig = ax
        
    Bphi = np.abs(Bphi)
    ripple = (np.amax(Bphi,axis=1) - np.amin(Bphi,axis=1)) / (np.amax(Bphi,axis=1) + np.amin(Bphi,axis=1))
    print(( np.amax(ripple),np.amin(ripple) ) )
    cont = ax.contour(R,z,np.log10(ripple))
    
    ax.set_title('Ripple', fontsize=FONTTITLE)
    ax.set_xlabel('R (m)', fontsize=FONTLABEL)
    ax.set_ylabel('z (m)', fontsize=FONTLABEL)
    ax.axis('image')
    plt.rc('font', size=FONTTICK)
    
    
def plotripplephi(Rvec, zvec, phivec, R, z, B, ax=None, fig=None):
    if ax == None:
        ax = plt.figure()
        fig = ax
        
    nphi = len(phivec)
    Bi = np.zeros(phivec.shape)
    for i in range(nphi):
        fB = sci.interp2d(Rvec, zvec, np.squeeze(B[:,i,:]) )
        Bi[i] = fB(R,z)
    
    ax.plot(phivec, Bi)
    ax.set_title('Ripple', fontsize=FONTTITLE)
    ax.set_xlabel('phi (deg)', fontsize=FONTLABEL)
    ax.set_ylabel('B (T)', fontsize=FONTLABEL)
    #ax.axis('tight')
    plt.rc('font', size=FONTTICK)
    