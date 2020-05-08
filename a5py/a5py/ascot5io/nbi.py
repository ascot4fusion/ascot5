"""
NBI injector HDF5 IO

File: nbi.py
"""
import numpy as np
import h5py

from . ascot5file import add_group
from a5py.ascot5io.ascot5data import AscotData

def write_hdf5(fn, nbi, desc=None):
    """
    Write NBI input in HDF5 file.

    Args:
        fn : str
            Full path to the HDF5 file.
        nbi : array
            Array of dictionaries describing each injector with fields:
                id : int
                    Numerical identifier for the injector
                nbeamlet : int
                    Number of beamlets in this injector
                beamletx : float
                    X coordinates of beamlets [m]
                beamlety : float
                    Y coordinates of beamlets [m]
                beamletz : float
                    Z coordinates of beamlets [m]
                beamletdx : float
                    X components of the unit direction vector of beamlets
                beamletdy : float
                    Y components of the unit direction vector of beamlets
                beamletdz : float
                    Z components of the unit direction vector of beamlets
                div_h : float
                    Horizontal divergence [radians]
                div_v : float
                    Vertical divergence [radians]
                div_halo_frac : float
                    Fraction of particles with halo divergence
                div_halo_h : float
                    Horizontal divergence in halo [radians]
                div_halo_v : float
                    Vertical divergence in halo [radians]
                anum : int
                    Mass number of injected species
                znum : int
                    Nuclear charge number of injected species
                mass : float
                    Mass of the injected species [kg]
                energy : float
                    Full injection energy [J]
                efrac : array_like (3)
                    Particle fractions for full, 1/2 and 1/3 energies
                power : float
                    Injected power [W]
        desc : str, optional <br>
            Input description.

    Returns:
        Name of the new input that was written.
    """

    parent = "nbi"
    group  = "nbi"
    gname  = ""

    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)
        gname = g.name.split("/")[-1]

        g.create_dataset("ninj",  (1,), data=len(nbi),    dtype="i4")

        for i in range(len(nbi)):
            ginj = g.create_group("inj"+str(i+1))

            nbeamlet = nbi[i]["nbeamlet"]

            ginj.create_dataset("id",            (1,),
                                data=nbi[i]["id"],            dtype="i4")
            ginj.create_dataset("nbeamlet",      (1,),
                                data=nbi[i]["nbeamlet"],      dtype="i4")
            ginj.create_dataset("beamletx",      (nbeamlet,),
                                data=nbi[i]["beamletx"],      dtype="f8")
            ginj.create_dataset("beamlety",      (nbeamlet,),
                                data=nbi[i]["beamlety"],      dtype="f8")
            ginj.create_dataset("beamletz",      (nbeamlet,),
                                data=nbi[i]["beamletz"],      dtype="f8")
            ginj.create_dataset("beamletdx",     (nbeamlet,),
                                data=nbi[i]["beamletdx"],     dtype="f8")
            ginj.create_dataset("beamletdy",     (nbeamlet,),
                                data=nbi[i]["beamletdy"],     dtype="f8")
            ginj.create_dataset("beamletdz",     (nbeamlet,),
                                data=nbi[i]["beamletdz"],     dtype="f8")
            ginj.create_dataset("div_h",         (1,),
                                data=nbi[i]["div_h"],         dtype="f8")
            ginj.create_dataset("div_v",         (1,),
                                data=nbi[i]["div_v"],         dtype="f8")
            ginj.create_dataset("div_halo_frac", (1,),
                                data=nbi[i]["div_halo_frac"], dtype="f8")

            ginj.create_dataset("div_halo_h", (1,), data=nbi[i]["div_halo_h"],
                                dtype="f8")
            ginj.create_dataset("div_halo_v", (1,), data=nbi[i]["div_halo_v"],
                                dtype="f8")
            ginj.create_dataset("anum",       (1,), data=nbi[i]["anum"],
                                dtype="i4")
            ginj.create_dataset("znum",       (1,), data=nbi[i]["znum"],
                                dtype="i4")
            ginj.create_dataset("mass",       (1,), data=nbi[i]["mass"],
                                dtype="f8")
            ginj.create_dataset("energy",     (1,), data=nbi[i]["energy"],
                                dtype="f8")
            ginj.create_dataset("efrac",      (3,), data=nbi[i]["efrac"],
                                dtype="f8")
            ginj.create_dataset("power",      (1,), data=nbi[i]["power"],
                                dtype="f8")

    return gname


def read_hdf5(fn, qid):
    """
    Read NBI input from HDF5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            QID of the data to be read.

    Returns:
        Dictionary containing input data.
    """

    path = "nbi/nbi_" + qid

    out = []
    with h5py.File(fn,"r") as f:
        ninj = f[path]["ninj"][0]
        for i in range(0,ninj):
            out.append({})
            for key in f[path+"/inj"+str(i+1)]:
                out[i][key] = f[path+"/inj"+str(i+1)][key][:]

    return out


def write_hdf5_dummy(fn):
    """
    Write a dummy injector.
    """

    nbi = {
        "id" : 1,
        "nbeamlet" : 1,
        "beamletx" : np.array([1.0]),
        "beamlety" : np.array([0.0]),
        "beamletz" : np.array([0.0]),
        "beamletdx" : np.array([1.0]),
        "beamletdy" : np.array([0.0]),
        "beamletdz" : np.array([0.0]),
        "div_h" : 0.0,
        "div_v" : 0.0,
        "div_halo_frac" : 0.0,
        "div_halo_h" : 0.0,
        "div_halo_v" : 0.0,
        "anum" : 1.0,
        "znum" : 1.0,
        "mass" : 1.0,
        "energy" : 1.0,
        "efrac" : [1,0,0],
        "power" : 1
    }
    return write_hdf5(fn, [nbi], desc="DUMMY")


class nbi(AscotData):
    """
    Object representing nbi data.
    """

    def read(self):
        return read_hdf5(self._file, self.get_qid())

    def write(self, fn, data=None):
        if data is None:
            data = self.read()

        return write_hdf5(fn, **data)
