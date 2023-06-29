"""
COM distribution module.
"""
import numpy as np
import h5py

import a5py.dist as distmod
from a5py.physlib.alias import getalias as alias

from .coreio.treedata import DataContainer

def write_hdf5(fn, run, data):
    """
    Write dist5D data in HDF5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
    """

    gname = "results/" + run + "/dist5d"

    with h5py.File(fn, "a") as f:
        g = f.create_group(gname)

        for i in range(0, len(data["abscissae"])):
            name = data["abscissae"][i]
            g.create_dataset("abscissa_nbin_0"+str(i+1), (1,),
                             data=data["n" + name], dtype="i4")
            abscissa = g.create_dataset("abscissa_vec_0"+str(i+1),
                                        (data["n" + name]+1,),
                                        data=data[name + "_edges"], dtype="f8")

            abscissa.attrs["name_0"+str(i)] = np.string_(name)
            abscissa.attrs["unit_0"+str(i)] = np.string_(data[name + "_unit"])

        g.create_dataset("abscissa_ndim", (1,), data=7, dtype="i4")
        g.create_dataset("ordinate_ndim", (1,), data=1, dtype="i4")

        ordinate = g.create_dataset(
            "ordinate", data=np.expand_dims(data["histogram"], axis=0),
            dtype="f8")
        ordinate.attrs["name_00"] = np.string_("density")
        ordinate.attrs["unit_00"] = np.string_(data["ordinate_unit"])


def read_hdf5(fn, qid):
    """
    Read 5D distribution from a HDF5 file to a dictionary.

    Args:
        fn : str <br>
            HDF5 file filename.
        qid : str <br>
            QID of the run whose distribution is read.
    Returns:
        Distribution dictionary.
    """

    with h5py.File(fn,"r") as f:

        path = "/results/run_"+qid+"/distcom/"
        dist = f[path]
        out = {}

        # A Short helper function to calculate grid points from grid edges.
        def edges2grid(edges):
            return np.linspace(0.5*(edges[0]+edges[1]),
                               0.5*(edges[-2]+edges[-1]), num=edges.size-1)

        abscissae = [0] * int(dist["abscissa_ndim"][:])
        for i in range(0, len(abscissae)):
            abscissa     = dist["abscissa_vec_0"+str(i+1)]
            name         = abscissa.attrs["name_0"+str(i)].decode("utf-8")
            abscissae[i] = name

            out[name + "_edges"] = abscissa[:]
            out[name + "_unit"]  = abscissa.attrs["unit_0"+str(i)].decode("utf-8")
            out[name]            = edges2grid(out[name + "_edges"])
            out["n" + name]      = out[name].size

        out["abscissae"] = abscissae
        out["histogram"] = dist["ordinate"][0,:,:,:,:,:,:,:]
        out["histogram_unit"] = dist["ordinate"].attrs["unit_00"].decode("utf-8")

    return out


class Dist_COM(DataContainer):
    """
    Object representing 5D distribution data.
    """

    def read(self):
        """
        Read distribution data from HDF5 file to a dictionary.

        Returns:
            Distribution dictionary.
        """
        with self as f:

            #path = "/results/run_"+qid+"/distcom/"
            #dist = f[path]
            dist = f
            out = {}

            # A Short helper function to calculate grid points from grid edges.
            def edges2grid(edges):
                return np.linspace(0.5*(edges[0]+edges[1]),
                                   0.5*(edges[-2]+edges[-1]), num=edges.size-1)

            abscissae = [0] * int(dist["abscissa_ndim"][:])
            for i in range(0, len(abscissae)):
                abscissa     = dist["abscissa_vec_0"+str(i+1)]
                name         = abscissa.attrs["name_0"+str(i)].decode("utf-8")
                abscissae[i] = name

                out[name + "_edges"] = abscissa[:]
                out[name + "_unit"]  = abscissa.attrs["unit_0"+str(i)].decode("utf-8")
                out[name]            = edges2grid(out[name + "_edges"])
                out["n" + name]      = out[name].size

                out["abscissae"] = abscissae
                out["histogram"] = dist["ordinate"][0,:,:,:]
                out["histogram_unit"] = dist["ordinate"].attrs["unit_00"].decode("utf-8")

        return out


    def write(self, fn, run, data=None):
        """
        Write dist5D data to HDF5 file.
        """
        if data is None:
            data = self.read()

        write_hdf5(fn, run, data)
