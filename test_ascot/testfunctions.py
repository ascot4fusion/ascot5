import numpy as np
from scipy import constants
from scipy.signal import argrelextrema

##################
#Unit conversions#
##################

#Atomic mass unit to kilograms
def amu2kg(val):
    amu = constants.physical_constants['atomic mass unit-kilogram relationship'][0]
    return np.multiply(amu, val)

#Electron volts to Joules
def eV2J(val):
    eV = constants.physical_constants['electron volt-joule relationship'][0]
    return np.multiply(eV, val)

#Elementary charge to Coulombs 
def e2C(val):
    return np.multiply(constants.e,val)

#Degrees to radians: use np.deg2rad

########################
#Coordinate converisons# 
########################

#Get cartesian x-coordinate from polar coordinates
def xFromPol(r, phi):
    return np.multiply(r,np.cos(phi))

#Get cartesian y-coordinate from polar coordinates
def yFromPol(r, phi):
    return np.multiply(r,np.sin(phi))

#Get cartesian components of vr = dr/dt 
#phi needs to be in deg
def vrToCart(vr, phi):
    i = np.multiply(np.cos(phi),vr)
    j = np.multiply(np.sin(phi),vr)
    return (i,j)

#Get cartesian components of vphi = dphi/dt
#phi needs to be in deg
def vphiToCart(vphi, phi):
    i = -1*np.multiply(np.sin(np.deg2rad(phi)),vphi)
    j = np.multiply(np.cos(np.deg2rad(phi)),vphi)
    return (i,j)

#############################
#Velocity/other vector stuff#
#############################

#Velocity in x and y directions 
def v(vr_pol, vphi_pol, phi):
    vr_cart = vrToCart(vr_pol, phi)
    vphi_cart = vphiToCart(vphi_pol, phi)
    i = vr_cart[0]+vphi_cart[0] #velocity in x direction
    j = vr_cart[1]+vphi_cart[1] #velocity in y direction
    return (i, j)

#Velocity on a cartesian 2D plane (or Rz-plane)
def v2D(v):
    return  np.sqrt(np.power(v[0],2)+np.power(v[1],2))

#Velocity on a cartesian 3D field
#v = Tuple of x and y velocity (returned by function v)
def v3D(v,z):
    vxy = v2D(v)
    return v2D((vxy,z))

#Average velocity
def v_avg(start,end,time):
    return np.divide(end-start,time)

#Uni vector of 3D vector
#v = vector in (x,y,z) tuple format 
def unit3(v):
    v_magn = v3D((v[0],v[1]),v[2])
    return np.transpose(v/v_magn)

#Unit vector of  2D vector
#v = (x,y)
def unit2(v):
    magn = v2D(v)
    return np.transpose(v/magn)

#Vector projection of v1 in the direction of v2
#Returns 2x3 tuple in the form ((x1,y1,z1),(x2,y2,z2))
def vectProjection(v1, v2):
    u2 = unit3(v2) #unit vector of v2
    scalar = np.dot(v1,u2) #scalar projection
    component = np.multiply(scalar,u2)
    rejection = np.subtract(v1,component)
    return (component, rejection)

####################
#Analytic solutions#
####################

#Larmor radius
def Larmor_a(m,v_xy,charge,B):
    return (m*v_xy)/(charge*B)

#Larmor radius when it isn't constant
def Larmor2_a(mass,v_perp,charge,mu):
    return (charge*v_perp)/(2*mass*mu)

#Period of one rotation 
#p_xy = momentum on xy_plane
def T_a(radius, p_xy, mass):
    return (2*np.pi*radius)/(p_xy/mass)

#Velocity of GC
#v_par = velocity parallel to the magnetic field
#b = magnetic field vector
#e = electric field vector
def vGC_a(v_par,b,e):
    return (np.multiply(np.array([0,0,v_par[0]]),b)+np.cross(e,b))[0]

def totalEnergy(mass,vpara,vperp):
    return 0.5*mass*(np.power(vpara,2)+np.power(vperp,2))

def totalEnergyMu(B,mu,mass,vpar):
    return B*mu+0.5*mass*np.power(vpar,2)

def canonicalMomentum(mass,R,v_para,charge,psi):
    return mass*R*v_para+charge*psi

#####Collisions#####

#Maxwell distribution (for thermal particles)
#E = energies in energy distribution
def MaxwellDistr(E):
    return np.sqrt(E)*np.exp(-E/1.0e4)

#Slowing down distribution (for fast particles) (@Wesson s.250)
#S = source rate
#E = energies in energy distribution (in J not eV)
#n_i = number of ions
#n_e = number of electrons
#T = temperature (in K not eV)
def slowingDownDistr(S,E,n_i,n_e,T,lnLambda):
    m_i = constants.m_p * n_i #mass of ions (given that there is only one proton per ion and no neutrons)
    m_e = constants.m_e * n_e #mass of electrons
    m_b = m_i + m_e           #mass of beam
    epsilon_c = np.power((3*np.sqrt(np.pi))/4,2/3) * np.power(m_i/m_e,1.0/3) * (m_b/m_i) * T
    A_D = (n*np.power(constants.e,4)*lnLambda) / (2*np.pi*np.power(constants.epsilon_0,2) * np.power(m_b,2))
    tau_se = (3*np.power(2*np.pi,0.5)*np.power(T,1.5)) / (np.power(m_e,0.5)*m_b*A_D) #slowing down time for electrons
    return (tau_se*S) / (2*E*(1+np.power(epsilon_c/E,3/2)))

#################
#ASCOT solutions#
#################

#x = x or y coordinates
#def Larmor(x):
#    xmax = np.amax(x)
#    xmin = np.amin(x)
#    xmid = xmin+(xmax-xmin)/2.0
#    return np.sqrt((xmax-xmid)/2.0)

#x = x coordinate, vy = velocity in y direction
#y coordinate and x velocity can be used too as a pair
def Larmor2(x, vy):
    vmax = argrelextrema(vy, np.greater)[0][0]
    vmin = argrelextrema(vy, np.less)[0][0]
    xmax = x[vmax]
    xmin = x[vmin]
    if xmin > xmax:
        temp = xmax
        xmax = xmin
        xmin = temp
    xmid = xmin+(xmax-xmin)/2.0
    return np.sqrt((xmax-xmid)**2.0) 

#works only with fixed timestep
def Larmor(x, y, T, time):
    t = int((T/2.0)/(time[0]-time[1]))  #number of timesteps it takes to make a half orbit
    return np.sqrt((x[t]-x[0])**2+(y[t]-y[0])**2)/2
    
#Period of rotation calculated from velocity or location
#v = velocity in either x or y direction (or location in the direction GC doesn't drift) 
def T(v, time):
    vmax = argrelextrema(v, np.greater)
    if vmax[0].size < 2:
        print 'Simulation time is too short, cannot calculate period of rotation'
    return np.multiply(vmax[0][0]-vmax[0][1],time[1]-time[0])

#Calculates GC drift velocity with the help of larmor radius
def vGClarmor(x, y, vr, vphi, phi, charge, radius, time):
    vxy = v(vr, vphi, phi) #v[0]:v_x, v[1]:v_y
    #velocity unit vectors
    vNorm = np.sqrt(np.diag(np.dot(np.transpose(vxy),vxy)))
    vUnit = np.divide(np.transpose(vxy), np.transpose(np.expand_dims(vNorm, axis = 0)))
    #unit vectors pointing to GC
    gcUnit = np.transpose(np.flipud(vUnit))
    gcHelper = np.array([[1],[-1]])
    if charge < 0: #Direction of unit vector depends on the charge of the particle
        gcHelper = np.array([[-1],[1]])
    gcUnit = np.multiply(gcUnit,gcHelper)
    #GCs location in each timestep
    X = np.array([x,y])+np.multiply(radius,gcUnit)
    #Average velocity of GC
    vGCavg = np.polyfit(time[::-1],X[0],1)
    return vGCavg[0]

#Calculates GC drift from the difference between ini and endstate GC transformation coordinates
def vGC(ini, end, time):
    return velocityAVG(ini,end,time[0])
