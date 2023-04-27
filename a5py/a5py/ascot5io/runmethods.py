"""
Wrapper class for RunNode that adds methods to plot and access results.

The motivation for this class is following. Consider that user wants to plot
a figure that shows i) particle orbits in (R,z), ii) their final (R,z) position,
and iii) the wall tiles. It is not clear whether this plot should be implemented
as a method for State or Orbit class as data from both are needed. Furthermore,
neither has access to the wall data. The only logical place for this method
therefore is in the RunNode that has access to all relevant data.

File: run.py
"""

import numpy as np

import a5py.plotting as a5plt

class LackingDataError(Exception):
    """
    Raised when the HDF5 file does not have necessary groups for the operation.
    """

    def __init__(self, missing):
        """
        Error informing what groups is missing.

        Args:
            missing : str <br>
                Name of the group that is missing.
        """
        msg = "Output does not contain " + missing + " which is required.\n"
        super().__init__(msg)


class RunMethods():

    def init_runmethods(self):
        self.has_inistate = False
        self.has_endstate = False
        self.has_orbit    = False

        if hasattr(self, "inistate"):
            self.has_inistate = True
        if hasattr(self, "endstate"):
            self.has_endstate = True
        if hasattr(self, "orbit"):
            self.has_orbit = True


    def getstate(self):
        pass


    def getorbit(self):
        if not has_orbit:
            pass


    def getorbit_poincareplanes(self):
        """
        Return a list of Poincaré planes that was collected in the simulation.

        The list consists of "pol", "tor", and "rad" depending on the type of
        the plane.
        """
        opt  = self.options.read()
        npol = opt["ORBITWRITE_POLOIDALANGLES"]
        ntor = opt["ORBITWRITE_TOROIDALANGLES"]
        nrad = opt["ORBITWRITE_RADIALDISTANCES"]
        npol = 0 if npol[0] < 0 else npol.size
        ntor = 0 if ntor[0] < 0 else ntor.size
        nrad = 0 if nrad[0] < 0 else nrad.size

        plane = [None] * (npol + ntor + nrad)
        plane[:npol]                    = ["pol"]
        plane[npol:npol+ntor]           = ["tor"]
        plane[npol+ntor:npol+ntor+nrad] = ["rad"]

        return plane


    def plotstate_scatter(
            self, x, y, z=None, c=None, color="C0", endcond=None,
            iniend=["i", "i", "i", "i"], log=[False, False, False, False],
            axesequal=False, axes=None, cax=None):
        """
        Make a scatter plot of marker state coordinates.

        Marker ini and end state contains the information of marker phase-space
        position at the beginning and at the end of the simulation. With this
        routine one can plot e.g. final (R,z) coordinates or mileage vs initial
        rho coordinate.

        The plot is either 2D+1D or 3D+1D, where the extra coordinate is color,
        depending on the number of queried coordinates.

        Args:
            x : str <br>
                Name of the quantity on x-axis.
            y : str <br>
                Name of the quantity on y-axis.
            z : str, optional <br>
                Name of the quantity on z-axis (makes plot 3D if given).
            c : str, optional <br>
                Name of the quantity shown with color scale.
            color : str, optional <br>
                Name of the color markers are colored with if c is not given.
            nc : int, optional <br>
                Number of colors used if c is given (the color scale is not
                continuous)
            cmap : str, optional <br>
                Name of the colormap where nc colors are picked if c is given.
            endcond : str, array_like, optional <br>
                Endcond of those  markers which are plotted.
            iniend : str, array_like, optional <br>
                Flag whether a corresponding [x, y, z, c] coordinate is taken
                from the ini ("i") or endstate ("e").
            log : bool, array_like, optional <br>
                Flag whether the corresponding axis [x, y, z, c] is made
                logarithmic. If that is the case, absolute value is taken before
                the quantity is passed to log10.
            axesequal : bool, optional <br>
                Flag whether to set aspect ratio for x and y (and z) axes equal.
            axes : Axes, optional <br>
                The Axes object to draw on. If None, a new figure is displayed.
            cax : Axes, optional <br>
                The Axes object for the color data (if c contains data),
                otherwise taken from axes.
        Raise:
            LackingDataError if queried state is missing from output.
        """
        if not self.has_inistate and "i" in iniend:
            raise LackingDataError("inistate")
        if not self.has_endstate and "e" in iniend:
            raise LackingDataError("endstate")

        # Decipher whether quantity is evaluated from ini or endstate
        state = [None] * 4
        for i in range(4):
            if iniend[i] == "i":
                state[i] = self.inistate
            elif iniend[i] == "e":
                state[i] = self.endstate
            else:
                raise ValueError("Only 'i' and 'e' are accepted for iniend.")

        # Function for processing given coordinate
        def process(crdname, crdindex):

            # Read data
            crd   = state[crdindex].get(crdname, endcond=endcond)
            label = crdname + " [" + str(crd.units) + "]"

            # Set logscale
            if log[crdindex]:
                crd   = np.log10(np.absolute(crd))
                label = "log10(| " + label + " |)"

            return (crd, label)


        # Process
        xc, xlabel = process(x, 0)
        yc, ylabel = process(y, 1)
        if z is not None:
            zc, zlabel = process(z, 2)
        if c is not None:
            cc, clabel = process(c, 3)
        else:
            cc = color
            clabel = None

        # Plot
        if z is None:
            a5plt.scatter2d(x=xc, y=yc, c=cc, xlabel=xlabel, ylabel=ylabel,
                            clabel=clabel, axesequal=axesequal, axes=axes,
                            cax=cax)
        else:
            a5plt.scatter3d(x=xc, y=yc, z=zc, c=cc, xlabel=xlabel,
                            ylabel=ylabel, zlabel=zlabel, clabel=clabel,
                            axesequal=axesequal, axes=axes, cax=cax)


    def plotstate_histogram():
        pass


    def plotorbit_poincare(self, plane, conlen=True, axes=None, cax=None):
        """
        Poincaré plot where color separates markers or shows connection length.

        Args:
            plane : int <br>
                Index number of the plane to be plotted in the list given by
                getorbit_poincareplanes.
            conlen : bool, optional <br>
                If true, trajectories of lost markers are colored (in blue
                shades) where the color shows the connection length at that
                position. Confined (or all if conlen=False) markers are shown
                with shades of red where color separates subsequent
                trajectories.
            axes : Axes, optional <br>
                The Axes object to draw on.
            cax : Axes, optional <br>
                The Axes object for the connection length (if conlen is True),
                otherwise taken from axes.
        Raise:
            LackingDataError if endstate or poincaré data is missing from
            output.
        """
        if not self.has_orbit:
            raise LackingDataError("orbit")
        try:
            self.orbit["pncrid"]
        except:
            raise LackingDataError("orbit(poincare)")
        if not self.has_endstate:
            raise LackingDataError("endstate")

        # Set quantities corresponding to the planetype
        planes = self.getorbit_poincareplanes()
        if len(planes) <= plane:
            raise ValueError("Unknown plane.")

        ptype = planes[plane]
        if ptype == "pol":
            x = "r"
            y = "z"
            xlabel = "R [m]"
            ylabel = "z [m]"
            xlim = None # TBD
            ylim = None # TBD
            axesequal = True
        elif ptype == "tor":
            x = "rho"
            y = "phimod"
            xlabel = "Normalized poloidal flux"
            ylabel = "Toroidal angle [deg]"
            xlim = [0, 1.1]
            ylim = [0, 360]
            axesequal = False
        elif ptype == "rad":
            x = "thetamod"
            y = "phimod"
            xlabel = "Poloidal angle [deg]"
            ylabel = "Toroidal angle [deg]"
            xlim = [0, 360]
            ylim = [0, 360]
            axesequal = True

        ids = self.orbit.get("id", pncrid=plane)
        x   = self.orbit.get(x,    pncrid=plane)
        y   = self.orbit.get(y,    pncrid=plane)

        if ptype == "pol":
            # Determine (R,z) limits by making sure the data fits nicely and
            # then rounding to nearest decimal.
            xlim = [np.floor(np.amin(x)*9) / 10, np.ceil(np.amax(x)*11) / 10]
            ylim = [np.floor(np.amin(y)*11) / 10, np.ceil(np.amax(y)*11) / 10]

        if conlen:
            # Calculate connection length using endstate (only markers that pass
            # this plane are considered)
            total = self.endstate.get("mileage", ids=np.unique(ids))

            # Find how many data points each marker has in orbit data arrays.
            # idx are indices of last data point for given id.
            idx = np.argwhere(ids[1:] - ids[:-1]).ravel()
            idx = np.append(idx, ids.size-1)
            # Getting the number of data points now is just a simple operation:
            idx[1:] -= idx[:-1]
            idx[0] += 1

            # Now we can repeat the total mileage by the number of data points
            # each marker has, so we can calculate the connection length at each
            # point as: conlen = total - mileage
            total  = np.repeat(total, idx)
            conlen = total - self.orbit.get("mileage", pncrid=plane)

            # Now set confined markers as having negative connection length
            lost1 = self.endstate.get("id", endcond="rhomax")
            lost2 = self.endstate.get("id", endcond="wall")

            idx = ~np.logical_or(np.in1d(ids, lost1), np.in1d(ids, lost2))
            conlen[idx] *= -1
            clabel = "Connection length [" + str(conlen.units) + "]"
        else:
            conlen = None
            clabel = None

        a5plt.poincare(x, y, ids, conlen=conlen, xlim=xlim, ylim=ylim,
                       xlabel=xlabel, ylabel=ylabel, clabel=clabel,
                       axesequal=axesequal, axes=axes, cax=cax)
