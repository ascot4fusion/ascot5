# -*- coding: utf-8 -*-
import matplotlib as plt
import numpy as np

FONTTITLE = 12
FONTLABEL = 10
FONTTICK  = 8

def plottemperature(rho, T, labels, ax=None, logy=False):
    if ax == None:
        ax = plt.figure()
        
    handles = ax.plot(rho,T)
    ax.set_title('Temperature', fontsize=FONTTITLE)
    ax.set_xlabel('rho', fontsize=FONTLABEL)
    ax.set_ylabel('T (eV)', fontsize=FONTLABEL)
    plt.rc('font', size=FONTTICK)
    
    ax.legend(handles, labels, loc='upper right')
    
    if not logy:
        ax.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    else:
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        
    return ax
    
    
def plotdensity(rho, n, labels, ax=None, logy=False):
    if ax == None:
        ax = plt.figure()
       
    #if np.shape(rho)[0] != np.shape(n)[0]:
     #   n = n.transpose()
        
    handles = ax.plot(rho, n)
    ax.set_title('Density', fontsize=FONTTITLE)
    ax.set_xlabel('rho', fontsize=FONTLABEL)
    ax.set_ylabel('n (m^-3)', fontsize=FONTLABEL)
    plt.rc('font', size=FONTTICK)
    
    ax.legend(handles, labels, loc='upper right')
    
    if not logy:
        ax.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    else:
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        
    return ax