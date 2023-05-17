'''
Created on Apr 14, 2020

@author: sjjamsa
'''

import copy
import numpy as np
import h5py

from a5py.ascot5io.ascot5data import AscotData
from a5py.ascot5io.ascot5file import read_data
from a5py.marker.plot import plot_histogram

from a5py.misc import openfigureifnoaxes
from a5py.physlib.species import species as getspecies

def read_hdf5(fn, qid, prefix):
    """
    Read particle marker input from HDF5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            QID of the data to be read.

    Returns:
        Dictionary containing input data.
    """

    path = "marker/"+prefix+"_" + qid

    out = {}
    with h5py.File(fn,"r") as f:
        for key in f[path]:
            # Make all read-in datasets 1-d arrays,
            # regardless of their original dimensionality
            d=f[path][key][:]
            out[key] = np.reshape(d,newshape=(d.size,))

    out["ids"] = out["id"]
    del out["id"]
    return out


def generatemrk(nmrk, mrktype, species=None):
    """
    Generate dummy marker input of given type and species.
    """
    mrk = {
        "n"      : nmrk,
        "ids"    : ( 1 + np.arange(nmrk) ),
        "r"      : np.zeros((nmrk,)),
        "z"      : np.zeros((nmrk,)),
        "phi"    : np.zeros((nmrk,)),
        "weight" : np.ones((nmrk,)),
        "time"   : np.zeros((nmrk,)),
    }
    if species is not None:
        species = getspecies(species)

    if mrktype == "particle":
        mrk["vphi"]   = np.zeros((nmrk,))
        mrk["vz"]     = np.zeros((nmrk,))
        mrk["vphi"]   = np.zeros((nmrk,))
        mrk["mass"]   = ( species["mass"]   * np.ones((nmrk,)) ).to_value("amu")
        mrk["charge"] = species["charge"] * np.ones((nmrk,), dtype=np.int16)
        mrk["anum"]   = species["anum"]   * np.ones((nmrk,), dtype=np.int16)
        mrk["znum"]   = species["znum"]   * np.ones((nmrk,), dtype=np.int16)
    if mrktype == "gc":
        mrk["pitch"]  = np.zeros((nmrk,))
        mrk["energy"] = np.zeros((nmrk,))
        mrk["zeta"]   = np.zeros((nmrk,))
        mrk["mass"]   = ( species["mass"]   * np.ones((nmrk,)) ).to_value("amu")
        mrk["charge"] = species["charge"] * np.ones((nmrk,), dtype=np.int16)
        mrk["anum"]   = species["anum"]   * np.ones((nmrk,), dtype=np.int16)
        mrk["znum"]   = species["znum"]   * np.ones((nmrk,), dtype=np.int16)
    if mrktype == "ml":
        mrk["pitch"]  = np.zeros((nmrk,))

    return mrk


class mrk(AscotData):
    '''
    A class acting as a superclass for all marker types.
    '''

    def prune(self,keepNmarkers, data=None, probabilities=None):
        '''
        Keep a subset of markers.

        Args:
            keepNmarkers : int <br>
                How many markers to keep.

            data=None         : dict <br>
                The data from .read() -method.
                If None, will be automatically read.

            probabilities=None
                What is the probability to include each marker
                If None, use the weights, (sum normalized to 1.0).

        Returns:
            A dictionary as if from .read()
        '''

        if data is None:
            data = self.read()
        else:
            data = copy.deepcopy(data)

        n = data['n'][0]
        totalWeight = np.sum(data["weight"])

        if probabilities is None:
            probabilities = data["weight"] / totalWeight

        fortune = np.random.choice(np.arange(n),size=(keepNmarkers,),replace=False,p=probabilities)

        for k in data.keys():
            if k=='n':
                data[k][0]=keepNmarkers
                continue
            data[k] = data[k][fortune]

        newWeight = np.sum(data["weight"])

        data['weight'] *= totalWeight / newWeight

        return data


    @openfigureifnoaxes(projection=None)
    def plot_hist_rhophi(self, ascotpy, rbins=10, pbins=10, weighted=False,
                         axes=None):
        """
        Plot marker rho-phi histogram
        """
        weights = None
        if weighted:
            weights = self.read()["weight"]

        rho = self.eval_rho(ascotpy)
        phi = self.eval_phi(ascotpy)
        plot_histogram(rho, xbins=rbins, y=phi, ybins=pbins,
                       weights=weights,
                       logscale=False, xlabel="Normalized poloidal flux",
                       ylabel="Toroidal angle [deg]",
                       axes=axes)


    @openfigureifnoaxes(projection=None)
    def plot_hist_energypitch(self, ascotpy, ebins=10, pbins=10, weighted=False,
                              axes=None):
        """
        Plot marker energy-pitch histogram

        Energy is on logscale.
        """
        weights = None
        if weighted:
            weights = self.read()["weight"]

        pitch  = self.eval_pitch(ascotpy)
        energy = self.eval_energy(ascotpy)
        if energy is None:
            # Field lines don't have energy
            plot_histogram(x=pitch, xbins=pbins, weights=weights,
                           logscale=False, xlabel=r"Pitch [$p_\parallel/p$]",
                           axes=axes)
        else:
            energy = np.log10(energy.to("eV"))
            plot_histogram(energy, xbins=ebins, y=pitch, ybins=pbins,
                           weights=weights,
                           logscale=False, xlabel="Energy [eV]",
                           ylabel=r"Pitch [$p_\parallel/p$]",
                           axes=axes)


    def eval_rho(self, ascotpy):
        """
        Evaluate rho coordinate.
        """
        with self as h5:
            r   = read_data(h5, "r")
            phi = read_data(h5, "phi")
            z   = read_data(h5, "z")
        rho = ascotpy.evaluate(r, phi, z, 0, "rho")
        return rho


    def eval_phi(self, ascotpy):
        """
        Evaluate toroidal angle.
        """
        with self as h5:
            phi = read_data(h5, "r")

        return phi


    def eval_energy(self, ascotpy):
        """
        Evaluate energy.

        Implement in subclass.
        """
        return None


    def eval_pitch(self, ascotpy):
        """
        Evaluate pitch.

        Implemented in subclass.
        """
        return None
