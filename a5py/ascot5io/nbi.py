"""Input data representing NBI injectors.

NBI injectors are used by BBNBI5 to generate NBI-ion source and calculate
shinethrough.
"""
import numpy as np
import h5py

import a5py.nbi.plot as plot

from .coreio.fileapi import add_group
from .coreio.treedata import DataGroup

class NBI(DataGroup):
    """Object representing nbi data.
    """

    def read(self):
        """Read data from HDF5 file.

        Returns
        -------
        data : dict
            Data read from HDF5 stored in the same format as is passed to
            :meth:`write_hdf5`.
        """
        fn   = self._root._ascot.file_getpath()
        path = self._path

        out = []
        with h5py.File(fn,"r") as f:
            ninj = f[path]["ninj"][0]
            for i in range(0,ninj):
                out.append({})
                for key in f[path+"/inj"+str(i+1)]:
                    out[i][key] = f[path+"/inj"+str(i+1)][key][:]

        return out

    def plot_grid_3D(self, ax=None):
        beams = self.read()

        for beam in beams:
            ax = plot.plot_scatter_3D(beam["beamletx"], beam["beamlety"], beam["beamletz"],
                                      equal = True, axes=ax,
                                      color="red",linewidth=0.75)
        return ax

    def plot_beamlet_3D(self, ax=None):
        beams = self.read()

        for beam in beams:
            ax = plot.plot_arrow_3D(beam["beamletx"], beam["beamlety"], beam["beamletz"],
                                    beam["beamletdx"], beam["beamletdy"], beam["beamletdz"],
                                    axes=ax, arrow_length_ratio=0,
                                    color="green", linewidth=0.1, length=10)

        return ax

    @staticmethod
    def write_hdf5(fn, nbi, desc=None):
        """Write input data to the HDF5 file.

        Parameters
        ----------
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
        desc : str, optional
            Input description.

        Returns
        -------
        name : str
            Name, i.e. "<type>_<qid>", of the new input that was written.
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

                ginj.create_dataset("div_halo_h", (1,),
                                    data=nbi[i]["div_halo_h"], dtype="f8")
                ginj.create_dataset("div_halo_v", (1,),
                                    data=nbi[i]["div_halo_v"], dtype="f8")
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

    @staticmethod
    def write_hdf5_dummy(fn):
        """Write dummy data that has correct format and is valid, but can be
        non-sensical.

        This method is intended for testing purposes or to provide data whose
        presence is needed but which is not actually used in simulation.

        The dummy output is a very large rectangular wall.

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.

        Returns
        -------
        name : str
            Name, i.e. "<type>_<qid>", of the new input that was written.
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
        return NBI.write_hdf5(fn, [nbi], desc="DUMMY")
