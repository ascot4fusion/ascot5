import numpy as np

ITM_EV           = 1.602176565e-19              # electron volt [eV]
ITM_AMU          = 1.660538921e-27          # atomic mass unit [kg]
ITM_ME           = 9.10938291e-31           # electron mass [kg]
ITM_MP           = 1.672621777e-27          # proton mass [kg]
ITM_MD           = 3.34358348e-27           # deuteron mass [kg]
ITM_MT           = 5.00735630e-27           # triton mass [kg]
ITM_MA           = 6.64465675e-27           # alpha mass [kg]
ITM_MN           = 1.674927351e-27          # neutron mass [kg]
ITM_MASS_He_3    = 3.0160293191 * ITM_AMU   # isotope mass [kg]

consts_e         = ITM_EV
consts_amu       = ITM_AMU
consts_mElectron = ITM_ME
consts_mProton   = ITM_MP
consts_mDeuteron = ITM_MD
consts_mTriton   = ITM_MT
consts_mAlpha    = ITM_MA
consts_mNeutron  = ITM_MN
consts_mHe3      = ITM_MASS_He_3

consts_Aelectron = 0
consts_AProton   = 1
consts_ANeutron  = 1
consts_ADeuteron = 2
consts_ATriton   = 3
consts_AHe3      = 3
consts_AAlpha    = 4

consts_Zelectron = -1
consts_Zneutron  = 0
consts_ZProton   = 1
consts_ZDeuteron = 1
consts_ZTriton   = 1
consts_ZHe3      = 2
consts_ZAlpha    = 2                                                            

def read_particles(fname):
    '''
    Read ASCOT4 input.particles file
    '''
    import numpy
    
    # Read Header part
    #-----------------
    with open(fname,'r') as f:
        headerLength = 0
        
        data={'filename':fname,'fields':{}}
        
        line = f.readline()
        headerLength += 1
        if line[:-1] != ' PARTICLE INITIAL DATA FOR ASCOT' and line[:-1] != ' PARTICLE OUTPUT DATA FROM ASCOT':
            print( '"'+line[:-1]+'"')
            raise NameError('Unrecognized first line in "'+fname+'".')
        
        line = f.readline()
        headerLength += 1
        if line.split()[0] != '4':
            print( '"'+line+'"')
            raise NameError('Bad file version in "'+fname+'".')
    
        line = f.readline()
        headerLength += 1
        
        line = f.readline()
        headerLength += 1
        nComments = int(line.split()[0])
        #print nComments
        
        comments=[]
        for i in range(0,nComments):
            comments.append(f.readline()[:-1])
            headerLength += 1
        data['comments']=comments
        
        line = f.readline()
        headerLength += 1
        
        line = f.readline()
        headerLength += 1
        nParticles = int(line.split()[0])
        #print nParticles
        
        line = f.readline()
        headerLength += 1
        
        line = f.readline()
        headerLength += 1
        nFields = int(line.split()[0])
        #print nFields
        
        fields = []
        for i in range(0,nFields):
            fields.append(f.readline()[:10].strip())
            headerLength += 1
        data['fieldNames']=fields
        
        line = f.readline()
        headerLength += 1
           
        #print headerLength
    
    # Read the actual data table
    # --------------------------
    
    print('Reading', nFields, 'fields for', nParticles, 'particles.')
    columns = numpy.loadtxt(fname,skiprows=headerLength)
    #nParticles = columns.shape[0]

    print('Read ', nParticles, 'particles')
    
    for i in range(0,nFields):
        data['fields'][data['fieldNames'][i]]=-999.0 * numpy.ones(nParticles)    

        
    for i in range(0,nParticles):
        for j in range(0,nFields):
            data['fields'][data['fieldNames'][j]][i] = columns[i,j]
        
    return data
            
def magneticMoment(Ekin,pitch,mass,B):
    v = Epitch2vparaperp(mass,Ekin,pitch)
    gammar = physics_eval.lorentzfactor(sqrt(v(1)^2+ v(2)^2))
    pperp = gammar*v(2)*mass
    val = pperp^2/(2*mass*B)
    return val
    
def Epitch2vparaperp(mass,Ekin,pitch):
    # (E,pitch) as (v_para,v_perp)
    magnv = physics_eval.Ekin2velocity(mass,Ekin);
    v_para = magnv*pitch;
    v_perp = magnv*sqrt(1-pitch^2);
    return (v_para, v_perp)

def Ekin2velocity(mass,Ekin):
    # Kinetic energy as velocity: v = sqrt( 1-1/( 1+E_kin/mc^2 )^2 )
    vperc = sqrt( 1-1/( 1+Ekin/(mass*physics_const.c^2) )^2 );
    val = physics_const.c*vperc;
    return val

# Return the mass of element with mass number A, charge number Z and
# a given charge. The mass is returned in [amu]
def guessMass(A,Z,charge):
    mass=-1.0
    
    if (Z == consts_Zelectron):
        if (charge < 0.0):
            mass = consts_mElectron
            return mass / consts_amu
    elif (Z == consts_Zneutron):
        if (A == consts_Aneutron):
            mass = consts_mNeutron
            return mass / consts_amu
    elif (Z == 1): # Hydrogen
        if (A == consts_AProton):
            mass = massPlusElectrons( consts_mProton, Z, charge )
            return mass / consts_amu
        elif (A == consts_ADeuteron):
            mass = massPlusElectrons( consts_mDeuteron, Z, charge )
            return mass / consts_amu
        elif (A == consts_ATriton):
            mass = massPlusElectrons( consts_mTriton, Z, charge )
            return mass / consts_amu
    elif (Z == 2): # Helium
        if (A == consts_AHe3):
            mass = massPlusElectrons( consts_mHe3, Z, charge )
            return mass / consts_amu
        elif (A == consts_AAlpha):
            mass = massPlusElectrons( consts_mAlpha, Z, charge )
            return mass / consts_amu
    # Hmm. we better just estimate with other  elements...ignoring the binding energy
    print("Warning! Unknown marker isotope, guesstimating mass.")
    mass = massPlusElectrons( Z*consts_mProton + (A-Z)*consts_mNeutron, Z, charge )
    return mass / consts_amu

def massPlusElectrons(massIn, Z, charge):
    if( int(round(charge/consts_e)) == Z):
        return massIn
    else:
        return massIn + consts_mElectron * (round(-charge/consts_e)-Z)
