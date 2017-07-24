import numpy as np
import matplotlib.pyplot as plt

def plot(asc):
    for i in asc["orbits"]:
        if "tor" in i and asc["orbits"][i]["N"] > 0:
            plt.figure()
            plt.plot(asc["orbits"][i]["rho"],asc["orbits"][i]["phi"])
            plt.show(block=False)

        elif "pol" in i and asc["orbits"][i]["N"] > 0:
            plt.figure()
            plt.plot(asc["orbits"][i]["R"],asc["orbits"][i]["z"])
            plt.show(block=False)

