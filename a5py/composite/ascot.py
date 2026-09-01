"""Contains the definition of the Ascot class which is the main interface.

Ascot class acts as an interface to other tools e.g. AscotData (interface to the
data), BeamRun (interface to the beam simulation), etc. It also acts as an
interface to input data interpolation. While each input data itself has the
tools to process it's own data, the logical place for processing data that
requires multiple inputs is here.
"""
import os

from a5py.data import AscotData
import a5py.engine.solve as solve
from a5py.engine.solve import *

class Ascot():
    """Primary tool for executing and processing ASCOT5 simulations and data.

    Attributes
    ----------
    data : :class:`.AscotData`
        Contains the simulation inputs and outputs.
    """

    def __init__(self, inputfile=None, create=False):
        """Initialize Ascot instance.

        Parameters
        ----------
        inputfile : str, optional
            Name of the HDF5 file or `None` to create an empty instance.
        create : bool, optional
            Create a new HDF5 file with given name.
        """
        super().__init__()

        self.data: AscotData = None
        if inputfile is None:
            self.data = AscotData()
        else:
            self.data = AscotData((inputfile, not create))


    def preflight_check(self, data, params, **priority_inputs):
        """Check that inputs are consistent and return only those inputs that are
        needed.

        Parameters
        ----------
        data : :class:`.AscotData`
            Simulation data structure whose `active` inputs will be used
        params : :class:`.Options`
            Simulation parameters.
        priority_inputs : dict
            Inputs to use in place of the ones that are `active` in ``data``.

        Returns
        -------
        inputs : dict[str, :class:`.InputVariant`]
            Inputs that will be used in the simulation.
        """
        self.preflight_check_parameters(params)

        required = {"bfield", "efield", "marker"}

        if params.physics.enable_icrh:
            required.add("rfof")

        if (params.physics.enable_coulomb_collisions
                or params.endconditions.activate_energy_limits):
            required.add("plasma")

        if params.physics.enable_mhd:
            required.update({"boozer", "mhd"})

        inputs = {}
        for req in required:
            try:
                inputs[req] = data[req].active
            except AscotDataException:
                if req in priority_inputs:
                    inputs[req] = priority_inputs[req]

        missing_inputs = {req for req in required if not req in inputs.keys()}
        if missing_inputs:
            raise ValueError(f"Missing inputs: {missing_inputs}")

        return inputs

    def preflight_check_parameters(self, params):
        """Check that options are internally consistent."""
        return True

    def simulate(
        self, time=None, run=None, params=None, comm=None, **priority_inputs,
        ):
        """Let the markers loose.

        The simulation can either be newly initialized or the simulation can be
        continued from a previous run.

        For new simulation, options must be provided. By default the inputs that are
        set as `active` are used but these can be overridden with ``**inputs``:

        .. code-block:: python

            simulate(options=Options(...), bfield=mybfield)

        For continuing a simulation, the run must be provided and no ``params`` or
        ``**inputs``:

        .. code-block:: python

            simulate(run=run)

        In both cases, ``time`` can be supplied to control the real time for how
        long the simulation should run.


        Parameters
        ----------
        run : Run
            The run to simulate.
        """
        root = True if comm is None else comm.Get_rank() == 0

        inputs = {}
        if root:
            inputs = self.preflight_check(self.data, params, **priority_inputs)

        run, inputs, unstage = setup(inputs, params, comm)
        execute(run, time)
        finalize(run, unstage, comm)
        self.data._treemanager.enter_leaf(run)

        return run

    def resume(self, continue_sim_time=0, continue_real_time=0, run=None):
        """Resume a simulation.

        This method takes an existing run and continues simulating markers that
        had end conditions maximum simulation time or maximum real time.

        Parameters
        ----------
        continue_sim_time : float, optional
            Increment maximum simulation time by this amount.
        continue_real_time : float, optional
            Increment maximum real time by this amount.
        run : Run, *optional*
            The run which is continued.

            By default, the active run is continued.
        """
        pass

    def solve_nbi_source(self):
        """Solve beam ion source distribution."""
        pass

    def solve_fusion_source(
            self,
            reaction,
            product_hist1=None,
            product_hist2=None,
            beam=None,
            another_beam=None,
            samples_per_cell=1000,
            store_samples=False,
            ):
        """Solve fusion source distribution.

        This function can calculate either thermal fusion, beam-thermal fusion
        (a.k.a. beam-target), or beam-beam fusion. Here "thermal" refers to
        (background) species that have Maxwellian distribution whereas the beam
        distribution can be arbitrary.

        Parameters
        ----------
        reaction : {"DT"}
            The modelled fusion reaction.
        product_hist1 : dict, optional
            Histogram basis of the first reaction product or None if not solved.
        product_hist2 : dict, optional
            Distribution of the second reaction product.
        """
        if phi is None:
            phi = np.array([0, 360])
        dont_use_rpz_basis = ( any([r is None, z is None]) and
                               all([r is None, z is None]) )
        use_rpz_basis = ( not any([r is None, z is None]) and
                          all([not r is None, not z is None]) )
        if use_rpz_basis and dont_use_rpz_basis:
            raise ValueError(
                "Either give all of r, phi, and z to use (R,phi,z) basis or "
                "use (rho,theta,phi) basis instead"
                )
        dont_use_ekinpitch_basis = ( any([ekin1 is None, pitch1 is None]) and
                                     all([ekin1 is None, pitch1 is None]) )
        use_ekinpitch_basis = ( not any([ekin1 is None, pitch1 is None]) and
                                all([not ekin1 is None, not pitch1 is None]) )
        if use_ekinpitch_basis and dont_use_ekinpitch_basis:
            raise ValueError(
                "Either give both ekin1 and pitch1 to use (E,pitch) basis or "
                "use (ppar,pperp) basis instead"
                )

        dont_use_ekinpitch_basis2 = ( any([ekin2 is None, pitch2 is None]) and
                                      all([ekin2 is None, pitch2 is None]) )
        use_ekinpitch_basis2 = ( not any([ekin2 is None, pitch2 is None]) and
                                 all([not ekin2 is None, not pitch2 is None]) )
        if use_ekinpitch_basis2 and dont_use_ekinpitch_basis2:
            raise ValueError(
                "Either give both ekin2 and pitch2 to use (E,pitch) basis or "
                "use (ppar,pperp) basis instead"
                )
        if use_ekinpitch_basis != use_ekinpitch_basis2:
            if use_ekinpitch_basis:
                raise ValueError(
                    "Both product distributions have to use same basis system: "
                    "now the first product were to use (E,pitch) and the "
                    "second (ppar,pperp)"
                    )
            else:
                raise ValueError(
                    "Both product distributions have to use same basis system: "
                    "now the second product were to use (E,pitch) and the "
                    "first (ppar,pperp)"
                    )
        self._ascot.input_init(bfield=True, plasma=True)
        if dont_use_rpz_basis:
            rho = np.linspace(0, 1, 10) if rho is None else rho
            theta = np.array([0, 180, 360]) if theta is None else theta
            vol, rc, phic, zc = self._ascot.input_rhovolume(
                nrho=rho.size, ntheta=theta.size, nphi=phi.size, method="prism",
                return_coords=True,
                )
            phic = phic.ravel()
        else:
            phic, rc, zc = np.meshgrid(1.5*phi[:-1]-0.5*phi[1:],
                                       1.5*r[:-1]-0.5*r[1:],
                                       1.5*z[:-1]-0.5*z[1:])
            phic *= np.pi/180
            vol = ( rc * np.diff(r[:2]) * np.diff(z[:2]) * np.diff(phi[:2])
                       * np.pi/180 )
        if dont_use_ekinpitch_basis:
            ppar1 = 1.3e-19 * np.linspace(-1., 1., 80) if ppar1 is None else ppar1
            pperp1 = np.linspace(0, 1.3e-19, 40) if pperp1 is None else pperp1
            ppar2 = 1.3e-19 * np.linspace(-1., 1., 80) if ppar2 is None else ppar2
            pperp2 = np.linspace(0, 1.3e-19, 40) if pperp2 is None else pperp2


    def thermal(
            self,
            reaction,
            r=None,
            phi=None,
            z=None,
            rho=None,
            theta=None,
            ppar1=None,
            pperp1=None,
            ekin1=None,
            pitch1=None,
            ppar2=None,
            pperp2=None,
            ekin2=None,
            pitch2=None,
            nmc=1000,
            ):
        """Calculate thermonuclear fusion between two thermal (Maxwellian)
        species.

        Parameters
        ----------
        reaction : int or str
            Fusion reaction index or name.
        r : array_like
            Abscissa for the radial coordinate in (R,phi,z) basis.
        phi : array_like
            Abscissa for the toroidal coordinate in (R,phi,z) and
            (rho,theta,phi) basis.
        z : array_like
            Abscissa for the axial coordinate in (R,phi,z) basis.
        rho : array_like
            Abscissa for the radial coordinate in (rho,theta,phi) basis.
        theta : array_like
            Abscissa for the poloidal coordinate in (rho,theta,phi) basis.
        ppar1 : array_like
            Abscissa for the parallel momentum in (ppar,pperp) basis for the
            first reaction product.
        pperp1 : array_like
            Abscissa for the perpendicular momentum in (ppar,pperp) basis for
            the first reaction product.
        ekin1 : array_like
            Abscissa for the kinetic energy in (E,pitch) basis for the first
            reaction product.
        pitch1 : array_like
            Abscissa for the pitch in (E,pitch) basis for the first reaction
            product.
        ppar2 : array_like
            Abscissa for the parallel momentum in (ppar,pperp) basis for the
            second reaction product.
        pperp2 : array_like
            Abscissa for the perpendicular momentum in (ppar,pperp) basis for
            the second reaction product.
        ekin2 : array_like
            Abscissa for the kinetic energy in (E,pitch) basis for the second
            reaction product.
        pitch2 : array_like
            Abscissa for the pitch in (E,pitch) basis for the second reaction
            product.
        nmc : int, optional
            Number of MC samples used in each (R, phi, z) bin.

        Returns
        -------
        prod1 : array_like
            Source distribution of the first fusion product.
        prod2 : array_like
            Source distribution of the second fusion product.
        """
        

        m1, q1, m2, q2, _, qprod1, _, qprod2, _ = self.reactions(reaction)
        reactions = {v: k for k, v in AFSI_REACTIONS.items()}
        reaction = reactions[reaction]
        anum1 = np.round(m1.to("amu").v)
        anum2 = np.round(m2.to("amu").v)
        znum1 = np.round(q1.to("e").v)
        znum2 = np.round(q2.to("e").v)
        q1 = np.round(qprod1.to("e").v)
        q2 = np.round(qprod2.to("e").v)

        nspec, _, _, anums, znums = self._ascot.input_getplasmaspecies()
        ispecies1, ispecies2 = np.nan, np.nan
        for i in np.arange(nspec):
            if( anum1 == anums[i] and znum1 == znums[i] ):
                ispecies1 = i
            if( anum2 == anums[i] and znum2 == znums[i] ):
                ispecies2 = i
        if np.isnan(ispecies1) or np.isnan(ispecies2):
            self._ascot.input_free(bfield=True, plasma=True)
            raise ValueError("Reactant species not present in plasma input.")
        mult = 0.5 if ispecies1 == ispecies2 else 1.0

        afsi = self._init_afsi_data(
            react1=ispecies1, react2=ispecies2, reaction=reaction, mult=mult,
            r=rc, phi=phic, z=zc, vol=vol,
            )

        if use_ekinpitch_basis:
            if use_rpz_basis:
                prod1 = self._init_histogram(
                    r, phi*np.pi/180, z, ekin1, pitch1, charge=q1, exi=True)
                prod2 = self._init_histogram(
                    r, phi*np.pi/180, z, ekin2, pitch2, charge=q2, exi=True)
            else:
                prod1 = self._init_histogram(
                    rho, theta*np.pi/180, phi*np.pi/180, ekin1, pitch1,
                    charge=q1, exi=True, toroidal=True)
                prod2 = self._init_histogram(
                    rho, theta*np.pi/180, phi*np.pi/180, ekin2, pitch2,
                    charge=q2, exi=True, toroidal=True)
        else:
            if use_rpz_basis:
                prod1 = self._init_histogram(
                    r, phi*np.pi/180, z, ppar1, pperp1, charge=q1, exi=False)
                prod2 = self._init_histogram(
                    r, phi*np.pi/180, z, ppar2, pperp2, charge=q2, exi=False)
            else:
                prod1 = self._init_histogram(
                    rho, theta*np.pi/180, phi*np.pi/180, ppar1, pperp1,
                    charge=q1, exi=False, toroidal=True)
                prod2 = self._init_histogram(
                    rho, theta*np.pi/180, phi*np.pi/180, ppar2, pperp2,
                    charge=q2, exi=False, toroidal=True)

        _LIBASCOT.afsi_run(ctypes.byref(self._ascot._sim),
                           ctypes.byref(afsi), nmc, prod1, prod2,
                           )
        self._ascot.input_free(bfield=True, plasma=True)

        # Reload Ascot
        self._ascot.file_load(self._ascot.file_getpath())
        return prod1, prod2

    def beamthermal(
            self,
            reaction,
            beam,
            swap=False,
            nmc=1000,
            ppar1=None,
            pperp1=None,
            ppar2=None,
            pperp2=None,
            ):
        """Calculate beam-thermal fusion.

        Parameters
        ----------
        reaction : int
            Fusion reaction index
        beam : dict
            Beam distribution that acts as the first reactant.
        swap : bool, optional
            If True, beam distribution acts as the second reactant and
            the first reactant is a background species.
        nmc : int, optional
            Number of MC samples used in each (R, phi, z) bin.
        ppar1 : array_like
            Abscissa for the parallel momentum in (ppar,pperp) basis for the
            first reaction product.
        pperp1 : array_like
            Abscissa for the perpendicular momentum in (ppar,pperp) basis for
            the first reaction product.
        ppar2 : array_like
            Abscissa for the parallel momentum in (ppar,pperp) basis for the
            second reaction product.
        pperp2 : array_like
            Abscissa for the perpendicular momentum in (ppar,pperp) basis for
            the second reaction product.

        Returns
        -------
        prod1 : array_like
            Fusion product 1 distribution.
        prod2 : array_like
            Fusion product 2 distribution.
        """
        ppar1 = 1.3e-19 * np.linspace(-1., 1., 50) if ppar1 is None else ppar1
        pperp1 = np.linspace(0, 1.3e-19, 50) if pperp1 is None else pperp1
        ppar2 = 1.3e-19 * np.linspace(-1., 1., 50) if ppar2 is None else ppar2
        pperp2 = np.linspace(0, 1.3e-19, 50) if pperp2 is None else pperp2

        m1, q1, m2, q2, _, qprod1, _, qprod2, _ = self.reactions(reaction)
        reactions = {v: k for k, v in AFSI_REACTIONS.items()}
        reaction = reactions[reaction]
        anum1 = np.round(m1.to("amu").v)
        anum2 = np.round(m2.to("amu").v)
        znum1 = np.round(q1.to("e").v)
        znum2 = np.round(q2.to("e").v)
        q1 = np.round(qprod1.to("e").v)
        q2 = np.round(qprod2.to("e").v)

        self._ascot.input_init(bfield=True, plasma=True)
        nspec, _, _, anums, znums = self._ascot.input_getplasmaspecies()
        ispecies = np.nan
        for i in np.arange(nspec):
            if( swap and anum1 == anums[i] and znum1 == znums[i] ):
                ispecies = i
                react1 = ispecies
                react2 = self._init_dist_5d(beam)
            if( not swap and anum2 == anums[i] and znum2 == znums[i] ):
                ispecies = i
                react2 = ispecies
                react1 = self._init_dist_5d(beam)
        if np.isnan(ispecies):
            self._ascot.input_free(bfield=True, plasma=True)
            raise ValueError("Reactant species not present in plasma input.")

        mult = 1.0
        r, z, phi = ( beam.abscissa_edges("r"), beam.abscissa_edges("z"),
                      beam.abscissa_edges("phi").to("rad") )
        phic, rc, zc = np.meshgrid(1.5*phi[:-1]-0.5*phi[1:],
                                       1.5*r[:-1]-0.5*r[1:],
                                       1.5*z[:-1]-0.5*z[1:])
        vol = ( rc * np.diff(r[:2]) * np.diff(z[:2]) * np.diff(phi[:2]) )
        afsi = self._init_afsi_data(
            react1=react1, react2=react2, reaction=reaction, mult=mult,
            r=rc, phi=phic, z=zc, vol=vol,
            )

        prod1 = self._init_histogram(
            beam.abscissa_edges("r"),
            beam.abscissa_edges("phi").to("rad"),
            beam.abscissa_edges("z"),
            ppar1,
            pperp1,
            charge=q1,
            exi=False,
            )
        prod2 = self._init_histogram(
            beam.abscissa_edges("r"),
            beam.abscissa_edges("phi").to("rad"),
            beam.abscissa_edges("z"),
            ppar2,
            pperp2,
            charge=q2,
            exi=False,
            )

        _LIBASCOT.afsi_run(ctypes.byref(self._ascot._sim),
                            ctypes.byref(afsi), nmc, prod1, prod2,
                            )
        self._ascot.input_free(bfield=True, plasma=True)

        # Reload Ascot
        self._ascot.file_load(self._ascot.file_getpath())
        return prod1, prod2

    def beambeam(
            self,
            reaction,
            beam1,
            beam2=None,
            nmc=1000,
            ppar1=None,
            pperp1=None,
            ppar2=None,
            pperp2=None,
            ):
        """Calculate beam-beam fusion.

        Parameters
        ----------
        reaction : int
            Fusion reaction index.
        beam1 : dict
            Beam1 distribution.
        beam2 : dict, optional
            Beam2 distribution or None to calculate fusion generation with
            beam1 itself.
        nmc : int, optional
            Number of MC samples used in each (R, phi, z) bin.
        ppar1 : array_like
            Abscissa for the parallel momentum in (ppar,pperp) basis for the
            first reaction product.
        pperp1 : array_like
            Abscissa for the perpendicular momentum in (ppar,pperp) basis for
            the first reaction product.
        ppar2 : array_like
            Abscissa for the parallel momentum in (ppar,pperp) basis for the
            second reaction product.
        pperp2 : array_like
            Abscissa for the perpendicular momentum in (ppar,pperp) basis for
            the second reaction product.

        Returns
        -------
        prod1 : array_like
            Fusion product 1 distribution.
        prod2 : array_like
            Fusion product 2 distribution.
        """
        ppar1 = 1.3e-19 * np.linspace(-1., 1., 80) if ppar1 is None else ppar1
        pperp1 = np.linspace(0, 1.3e-19, 40) if pperp1 is None else pperp1
        ppar2 = 1.3e-19 * np.linspace(-1., 1., 80) if ppar2 is None else ppar2
        pperp2 = np.linspace(0, 1.3e-19, 40) if pperp2 is None else pperp2
        _, _, _, _, _, qprod1, _, qprod2, _ = self.reactions(reaction)
        reactions = {v: k for k, v in AFSI_REACTIONS.items()}
        reaction = reactions[reaction]
        q1 = np.round(qprod1.to("e").v)
        q2 = np.round(qprod2.to("e").v)

        self._ascot.input_init(bfield=True, plasma=True)

        r, z, phi = ( beam1.abscissa_edges("r"), beam1.abscissa_edges("z"),
                      beam1.abscissa_edges("phi").to("rad") )
        phic, rc, zc = np.meshgrid(1.5*phi[:-1]-0.5*phi[1:],
                                       1.5*r[:-1]-0.5*r[1:],
                                       1.5*z[:-1]-0.5*z[1:])
        vol = ( rc * np.diff(r[:2]) * np.diff(z[:2]) * np.diff(phi[:2]) )

        react1 = self._init_dist_5d(beam1)
        if beam2 is not None:
            react2 = self._init_dist_5d(beam2)
            mult = 1.0
        else:
            react2 = react1
            mult = 0.5
        afsi = self._init_afsi_data(
            react1=react1, react2=react2, reaction=reaction, mult=mult,
            r=rc, phi=phic, z=zc, vol=vol,
            )

        prod1 = self._init_histogram(
            beam1.abscissa_edges("r"),
            beam1.abscissa_edges("phi").to("rad"),
            beam1.abscissa_edges("z"),
            ppar1,
            pperp1,
            charge=q1,
            exi=False,
            )
        prod2 = self._init_histogram(
            beam1.abscissa_edges("r"),
            beam1.abscissa_edges("phi").to("rad"),
            beam1.abscissa_edges("z"),
            ppar2,
            pperp2,
            charge=q2,
            exi=False,
            )
        _LIBASCOT.afsi_run(ctypes.byref(self._ascot._sim),
                           ctypes.byref(afsi), nmc, prod1, prod2,
                           )

        self._ascot.input_free(bfield=True, plasma=True)

        # Reload Ascot
        self._ascot.file_load(self._ascot.file_getpath())
        return prod1, prod2