import numpy as np
from a5py.postprocessing.physicslib import guessMass

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
