# -*- coding: utf-8 -*-
import matplotlib as plt
import numpy as np

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