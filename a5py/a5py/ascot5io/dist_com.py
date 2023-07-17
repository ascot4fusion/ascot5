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

        
    def plot_muptor(self, E, axes=None):
        """
        Plot constant of motion distribution on (mu, Ptor) plane for a given E.
        """
        import matplotlib.pyplot as plt
        data = self.read()
        hist = data["histogram"]
        mu_edges, Ptor_edges = data["mu_edges"], data["ptor_edges"]
        E_vector, E_edges = data["ekin"], data["ekin_edges"]
        if E is None:
            plot = np.sum(hist, axis=1)
        elif type(E) == int:
            i = E
            plot = hist[:,i,:]
        else:
            i = np.argmax(E < E_edges)-1
            plot1 = hist[:,i,:]
            plot2 = hist[:,i+1,:]

            k = E/(E_vector[i+1]-E_vector[i])
            plot = plot1*k + (1-k)*plot2
        
        h = axes.pcolormesh(Ptor_edges, mu_edges, plot, shading="flat")
        axes.set_ylabel("mu (J/T)")
        axes.set_xlabel("Ptor (Js)")
        plt.colorbar(h, ax=axes)
        
    def plot_muEkin(self, Ptor, axes=None):
        """
        Plot constant of motion distribution on (mu, Ptor) plane for a given E.
        """
        import matplotlib.pyplot as plt
        data = self.read()
        hist = data["histogram"]
        mu_edges, E_edges = data["mu_edges"], data["ekin_edges"]
        Ptor_vector, Ptor_edges = data["ptor"], data["ptor_edges"]
        if Ptor is None:
            plot = np.sum(hist, axis=2)
        elif type(Ptor) == int:
            i = Ptor
            plot = hist[:,:,i]
        else:
            i = np.argmax(Ptor < Ptor_edges)-1
            plot1 = hist[:,:,i]
            plot2 = hist[:,:,i+1]

            k = Ptor/(Ptor_vector[i+1]-Ptor_vector[i])
            plot = plot1*k + (1-k)*plot2
        
        h = axes.pcolormesh(E_edges, mu_edges, plot, shading="flat")
        axes.set_ylabel("mu (J/T)")
        axes.set_xlabel("E (J)")
        plt.colorbar(h, ax=axes)

        
    def plot_EkinPtor(self, mu, axes=None):
        """
        Plot constant of motion distribution on (mu, Ptor) plane for a given E.
        """
        import matplotlib.pyplot as plt
        data = self.read()
        hist = data["histogram"]
        E_edges, Ptor_edges = data["ekin_edges"], data["ptor_edges"]
        mu_vector, mu_edges = data["mu"], data["mu_edges"]
        if mu is None:
            plot = np.sum(hist, axis=0)
        elif type(mu) == int:
            i = mu
            plot = hist[i,:,:]
        else:
            i = np.argmax(mu < mu_edges)-1
            plot1 = hist[i,:,:]
            plot2 = hist[i+1,:,:]

            k = mu/(mu_vector[i+1]-mu_vector[i])
            plot = plot1*k + (1-k)*plot2
        
        h = axes.pcolormesh(Ptor_edges, E_edges, plot, shading="flat")
        axes.set_ylabel("E (J)")
        axes.set_xlabel("Ptor (Js)")
        plt.colorbar(h, ax=axes)
