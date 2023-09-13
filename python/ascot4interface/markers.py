import numpy as np
from a5py.postprocessing.physicslib import guessMass
import a5py.ascot5io.mrk_gc
import a5py.ascot5io.mrk_prt
import time

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
        nComments = int(float(line.split()[0]))

        comments=[]
        for i in range(0,nComments):
            comments.append(f.readline()[:-1])
            headerLength += 1
        data['comments']=comments

        line = f.readline()
        headerLength += 1

        line = f.readline()
        headerLength += 1
        nParticles = int(float(line.split()[0]))

        line = f.readline()
        headerLength += 1

        line = f.readline()
        headerLength += 1
        nFields = int(float(line.split()[0]))

        fields = []
        for i in range(0,nFields):
            fields.append(f.readline()[:10].strip())
            headerLength += 1
        data['fieldNames']=fields

        line = f.readline()
        headerLength += 1

    # Read the actual data table
    # --------------------------

    print('Reading', nFields, 'fields for', nParticles, 'particles.')
    columns = numpy.loadtxt(fname,skiprows=headerLength)
    
    # -1 means unkown number of lines.
    if nParticles == -1:
        nParticles = columns.shape[0]

    print('Read ', nParticles, 'particles')

    for i in range(0,nFields):
        data['fields'][data['fieldNames'][i]]=-999.0 * numpy.ones(nParticles)

    for i in range(0,nParticles):
        for j in range(0,nFields):
            data['fields'][data['fieldNames'][j]][i] = columns[i,j]

    if 'charge' not in data['fieldNames'] and 'Znum' in data['fieldNames']:
        print("Converting Znum to charge.")
        data["fields"]['charge'] = data["fields"]['Znum'].astype('float')
    if 'mass' not in data['fieldNames'] and 'Anum' in data['fieldNames']:
        print("Converting Anum to mass.")
        data["fields"]['mass'] = np.array(
            list(map(guessMass,
                     data["fields"]['Anum'],
                     data["fields"]['Znum'],
                     data["fields"]['charge'])))
    if 'id' not in data['fieldNames']:
        print("Generating unique ids.")
        data["fields"]['id'] = np.array(
            range(1,data["fields"]['charge'].size + 1))
    if (min(data["fields"]["id"]) <= 0):
        zero_ind = np.where(data["fields"]["id"] == 0)[0];
        data["fields"]["id"][zero_ind] = max(data["fields"]["id"] ) + 1
        print("Converting id 0 to new unique id: "
              + str(int(max(data["fields"]["id"]))))

    return data



def write_particles(filename,dataObject):
    ''' 
        Write ASCOT4 input.particles 
        @param filename: The file to write
        @param dataObject: The particle data to write, e.g. a5py.ascot5io.mrk_gc.mrk_gc 
    '''
        
    if isinstance(dataObject, a5py.ascot5io.mrk_gc.mrk_gc):
        print('Writing guiding centres')
        
        d=dataObject.read()
        # ascot5 name, ascot4 name, conversion factor
        fieldTuples=[('anum','Anum',    1.0),
                     ('charge','charge',1.0),
                     ('mass','mass',    1.0),
                     ('znum','Znum',    1.0),
                     ('r','R',          1.0),
                     ('phi','phi',      1.0),
                     ('z','z',          1.0),
                     ('pitch','pitch',  1.0),
                     ('energy','energy',1.0),
                     ('time','time',    1.0),
                     ('weight','weight',1.0),
                     ('ids','id',       1.0),   ]
        
        data={'fields':{}}
        for t in fieldTuples:
            data['fields'][t[1]] = t[2] * d[t[0]]
            #print(t[0], d[t[0]].shape)
        
    elif isinstance(dataObject, a5py.ascot5io.mrk_prt.mrk_prt):
        raise ValueError('Unsupported dataObject to write.')
        print('Writing particles')
        
    else:
        raise ValueError('Unknown dataObject to write.')
        
        
        
    data['comments']=[time.strftime('Written with a5py on %Y-%m-%d at %H:%M:%S%z.'),
                      'The input file and group are "{}" : "{}"'.format(dataObject._file, dataObject._path),
                      'Desc: "{}"'.format(dataObject.get_desc())]

    
    
    
    nFields   = len(data['fields'])
    fields    = list(data['fields'])
    nParticles= len(data['fields'][fields[0]])
    
    header=''' PARTICLE INITIAL DATA FOR ASCOT
 4 VERSION =====================
 
 {0}  # Number of comment lines, max length 256 char, (defined in prtRead_lineLength)
'''.format(len(data['comments']))

    # ADD THe COMMENTS
    
    for i in range(0,len(data['comments'])):
        header= header + data['comments'][i] + '\n'

    # ADD NUMBER OF PARTICLES
        
    header = header + '\n{0} # Number of particles (-1 means unknown number)\n\n'.format(nParticles)

    
    # ADD THE FIELD LIST
    
    header =  header + '{0} # Number of different fields for each particle [10 first letters are significant]\n'.format(nFields)

    fmt = []
    D   = np.empty((nParticles,0))
    
    descs={ 'Anum'   : 'mass number of particle        (integer)',
            'mass'   : 'mass of the particle           (amu)',
            'Znum'   : 'charge number of particle      (integer)',
            'charge' : 'charge of the particle         (elemental charge)',
            'energy' : 'kinetic energy of particle     (eV)',
            'rho'    : 'normalized poloidal flux       (real)',
            'phi'    : 'toroidal angle of particle     (deg)',
            'R'      : 'R of particle                  (m)',
            'z'      : 'z of particle                  (m)',
            'pitch'  : 'pitch angle cosine of particle (vpar/vtot)',
            'origin' : 'origin of the particle         ()',
            'weight' : 'weight factor of particle      (particle/second)',
            'id'     : 'unique identifier of particle  (integer)',
            'Tmax'   : 'maximum time to follow the prt (s)',
            'time'   : 'elapsed time of the particle   (s)',
            }

    
    #( D.shape, np.transpose(np.array([data['fields'][fields[0]]])).shape)
    integerFields = ['Anum','Znum','id','origin']
    for i in range(0,nFields):
        description = descs.get(fields[i],' ---')
        header = header + '{:10s}- {}\n'.format(fields[i],description)
        if any( fields[i] in s   for s in integerFields ):
            fmt.append('%5d')
        else:
            fmt.append('%20.10e')
        D = np.hstack((D,np.transpose(np.array([data['fields'][fields[i]]]))))
    
    footer = '#EOF'
    
    
    #print '"' + header + '"' 
    #print fmt
    #print D.shape
    #print footer
    print ('writing', nFields, 'fields for', nParticles,'particles to "'+filename+'".')
    np.savetxt(filename,D,fmt=fmt,delimiter=' ',comments='',header=header,footer=footer)


