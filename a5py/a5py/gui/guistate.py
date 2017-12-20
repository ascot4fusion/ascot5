import matplotlib as plt
import numpy as np

FONTTITLE = 12
FONTLABEL = 10
FONTTICK  = 8

def plotEkin(E, ax=None, logy=False, weights=None, nbin=50):
    if ax == None:
        ax = plt.figure()
        
    ax.hist(E, nbin, normed=True, stacked=True, log=logy, weights=weights, rwidth=1)
    ax.set_title('Energy', fontsize=FONTTITLE)
    ax.set_xlabel('Ekin (eV)', fontsize=FONTLABEL)
    plt.rc('font', size=FONTTICK)
    
    if not logy:
        ax.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    else:
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        
    return ax

def plotpitch(xi, ax=None, logy=False, weights=None, nbin=50):
    if ax == None:
        ax = plt.figure()

    ax.hist(xi, nbin, normed=True, stacked=True, log=logy, weights=weights, rwidth=1)
    ax.set_title('Pitch', fontsize=FONTTITLE)
    ax.set_xlabel('pitch = v_par/v', fontsize=FONTLABEL)
    ax.set_ylabel('distribution (a.u.)', fontsize=FONTLABEL)
    ax.set_xlim(-1,1)
    plt.rc('font', size=FONTTICK)
    
    if not logy:
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    
    return ax

def plottime(time, ax=None, logy=False, weights=None, nbin=50):
    if ax == None:
        ax = plt.figure()
        
    ax.hist(time, nbin, normed=True, stacked=True, log=logy, weights=weights, rwidth=1)
    ax.set_title('Time', fontsize=FONTTITLE)
    ax.set_xlabel('time (s)', fontsize=FONTLABEL)
    ax.set_ylabel('distribution (a.u.)', fontsize=FONTLABEL)
    plt.rc('font', size=FONTTICK)
    
    if not logy:
        ax.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    else:
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        
    return ax

def plotcputime(cputime, ax=None, logy=False, weights=None, nbin=50):
    if ax == None:
        ax = plt.figure()
        
    ax.hist(cputime, nbin, normed=True, stacked=True, log=logy, weights=weights, rwidth=1)
    ax.set_title('CPU-time', fontsize=FONTTITLE)
    ax.set_xlabel('time (s)', fontsize=FONTLABEL)
    ax.set_ylabel('distribution (a.u.)', fontsize=FONTLABEL)
    plt.rc('font', size=FONTTICK)
    
    if not logy:
        ax.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    else:
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    
    return ax

def plotcharge(charge, ax=None, logy=False, weights=None, nbin=50):
    if ax == None:
        ax = plt.figure()
        
    ax.hist(charge, nbin, normed=True, stacked=True, log=logy, weights=weights, rwidth=1)
    ax.set_title('Charge', fontsize=FONTTITLE)
    ax.set_xlabel('charge (e)', fontsize=FONTLABEL)
    ax.set_ylabel('distribution (a.u.)', fontsize=FONTLABEL)
    plt.rc('font', size=FONTTICK)
    
    if not logy:
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    return ax

def plotphi(phi, ax=None, logy=False, weights=None, nbin=50):
    if ax == None:
        ax = plt.figure()

    for i in range(len(phi)):
        phi[i] = np.mod(phi[i],360)
        
    ax.hist(phi, nbin, normed=True, stacked=True, log=logy, weights=weights, rwidth=1)
    ax.set_title('Toroidal angle', fontsize=FONTTITLE)
    ax.set_xlabel('phi (deg)', fontsize=FONTLABEL)
    ax.set_ylabel('distribution (a.u.)', fontsize=FONTLABEL)
    ax.set_xlim(0,360)
    ax.set_xticks(np.arange(0,361,60))
    plt.rc('font', size=FONTTICK)
    
    if not logy:
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    
    return ax

def plotendcond(endcond, endcondlabel, ax=None, logy=False, weights=None):
    if ax == None:
        ax = plt.figure()

    bins = 0
    ticks = []
    for i in range(len(endcond)):
        bins = bins + 1
        ticks.append(endcond[i][0][0])
        
    ax.hist(endcond, bins, normed=False, stacked=True, log=logy, weights=weights, rwidth=1)
    ax.set_title('End condition', fontsize=FONTTITLE)
    ax.set_ylabel('Number of particles/markers', fontsize=FONTLABEL)
    plt.rc('font', size=FONTTICK)
    
    
    if not logy:
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    ax.set_xticks(ticks)
    ax.set_xticklabels(endcondlabel)
    return ax

def plotrho(rho, ax=None, logy=False, weights=None, nbin=50):
    if ax == None:
        ax = plt.figure()
        
    ax.hist(rho, nbin, normed=True, stacked=True, log=logy, weights=weights, rwidth=1)
    ax.set_title('Radial position', fontsize=FONTTITLE)
    ax.set_xlabel('rho', fontsize=FONTLABEL)
    ax.set_ylabel('distribution (a.u.)', fontsize=FONTLABEL)
    plt.rc('font', size=FONTTICK)
    
    if not logy:
        ax.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    else:
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

    return ax

def plotRz(R,z, ax=None):
    if ax == None:
        ax = plt.figure()
    for i in range(len(R)):
        ax.plot(R[i],z[i],'.')
        
    ax.set_title('Poloidal position', fontsize=FONTTITLE)
    ax.set_xlabel('R (m)', fontsize=FONTLABEL)
    ax.set_ylabel('z (m)', fontsize=FONTLABEL)
    ax.axis('equal')
    plt.rc('font', size=FONTTICK)
    
    return ax
