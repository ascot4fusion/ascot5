import numpy as np
import matplotlib.pyplot as plt
import colormaps as cmaps

def plot(asc):
    colors = [ cmaps.viridis(x) for x in np.linspace(0.0, 1.0, 10) ]
    colors = (np.array(colors))[np.random.permutation(10).astype('i8')]

    for val in asc["orbits"]:
        if "tor" in val and asc["orbits"][val]["N"] > 0:
            plt.figure()
            ids = asc["orbits"][val]["id"]
            uniqueId = asc["orbits"][val]["uniqueId"]
            for i in (np.linspace(0,asc["orbits"][val]["N"]-1)).astype('i8'):
                color = colors[np.mod(i,10)]
                plt.plot(asc["orbits"][val]["rho"][ids==uniqueId[i]],np.mod(asc["orbits"][val]["phi"][ids==uniqueId[i]],360),marker='.',linestyle='none',color=color)
                
            plt.show(block=False)


        elif "pol" in val and asc["orbits"][val]["N"] > 0:
            plt.figure()
            ids = asc["orbits"][val]["id"]
            uniqueId = asc["orbits"][val]["uniqueId"]
            for i in (np.linspace(0,asc["orbits"][val]["N"]-1)).astype('i8'):
                color = colors[np.mod(i,10)]
                plt.plot(asc["orbits"][val]["R"][ids==uniqueId[i]],asc["orbits"][val]["z"][ids==uniqueId[i]],marker='.',linestyle='none',color=color)
                
            plt.show(block=False)

