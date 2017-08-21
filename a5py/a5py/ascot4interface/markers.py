import numpy as np

import sys
sys.path.append('../ui')
import ui_markers as markers

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
            print '"'+line[:-1]+'"'
            raise NameError('Unrecognized first line in "'+fname+'".')
        
        line = f.readline()
        headerLength += 1
        if line.split()[0] != '4':
            print '"'+line+'"'
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
    
    print 'Reading', nFields, 'fields for', nParticles, 'particles.'
    columns = numpy.loadtxt(fname,skiprows=headerLength)
    nParticles = columns.shape[0]

    print 'Read ', nParticles, 'particles'
    
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

def write_particles(fn, data):
    if 'vphi' in data['fieldNames']:
        # We have particles
        data = data["fields"]
        markers.write_hdf5_particles(fn, data["id"], data["mass"], data["charge"], data["Rprt"], data["phiprt"], data["zprt"], data["vR"], data["vphi"], data["vz"], data["weight"], data["weight"]*0)
    elif 'energy' in data['fieldNames']:
        # We have guiding centers
        markers.write_hdf5_guidingcenters(fn, data["id"], data["mass"], data["charge"], data["R"], data["phi"], data["z"], data["energy"], data["pitch"], data["theta"]*0, data["weight"], data["weight"]*0)
