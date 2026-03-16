# coding: utf-8

# For Simppa at least, following recipe worked to get the below working:
# 0. compile ascotlib.so and add the path to LD_LIBRARY_PATH
# 1. in ascot5.yml, set the python version: python=3.6, and create the environment
# 2. module load imasenv
# 3. conda activate ascot5
# the below works.

import ctypes
from   a5py import Ascot
import a5py.ascot5io.imas
import a5py.ascot5io.wall
from   a5py.ascotpy import ascot2py
from   a5py.ascot5io.marker import Marker
from a5py.ascot5io import options
from a5py.templates.optiontemplates import OptionTemplates
import numpy as np
from types import SimpleNamespace
import os

np.int = int
np.float=float

# Read from IMAS

def distsource_run(equil_b3d, equil, core_prof, distr_sour,wall2d,wall3d,distributions,optionsxml):
    ids=SimpleNamespace()

    mrkr=a5py.ascot5io.imas.marker()
    eq=a5py.ascot5io.imas.B_2DS()
    p1d=a5py.ascot5io.imas.plasma_1d()
    w2d=a5py.ascot5io.imas.wall_2d()
    w3d=a5py.ascot5io.imas.wall_3d()
    bsts=a5py.ascot5io.imas.B_STS()

    print('Parsing distsource')
    # markers from distsource
    ids.distribution_sources=distr_sour
    mrkr.setIds(ids)
    mrkrdict2=mrkr.parse()

    print('Parsing plasma1d')
    #equilibrium
    ids.equilibrium=equil
    eq.setIds(ids)
    eqdict=eq.parse()
    ids.core_profiles=core_prof
    p1d.setIds(ids) 
    p1d_dict=p1d.parse(equilibrium_ids=eq)

    print('Reading 2D wall')
    # 2D wall
    ids.wall=wall2d
    w2d.setIds(ids)
    wdict2=w2d.parse()

    print('Reading 3D wall')
    # 3D wall
    ids.wall=wall3d
    w3d.setIds(ids)
    wdict3=w3d.parse()

    print('Reading Bfield (stellarator)')
    # 3D Bfield in Stellarator format
    ids.equilibrium=equil_b3d
    bsts.setIds(ids)
    bdict=bsts.parse()
    
    write_to_hdf5=False
    if write_to_hdf5:
        print("writing markers to 'from_imas.h5'.")
        if 'energy' in mrkrdict2:
            a5py.ascot5io.marker.GC.write_hdf5('from_imas.h5',**mrkrdict2)
        elif 'vz' in mrkrdict2:
            a5py.ascot5io.marker.Prt.write_hdf5('from_imas.h5',**mrkrdict2)
        else:
            raise(ValueError("What kind of marker is this supposed to be?"))


        a5py.ascot5io.wall_2D.write_hdf5('from_imas.h5',**wdict2)
        a5py.ascot5io.wall_3D.write_hdf5('from_imas.h5',**wdict3)
        a5py.ascot5io.B_STS.write_hdf5('from_imas.h5',**bdict)


    del bdict['desc']

    print('Initializing ascot5')
    filepath='./helloworld.h5'
    if os.path.isfile(filepath):
        a5 = Ascot('helloworld.h5')
    else:
        a5=Ascot("helloworld.h5", create=True)

    print('Initializing inputs')



    bfield=bdict # True
    wall=wdict3  # wdict2, True
    edict = a5.data.create_input("E_TC", dryrun=True)
    pdict   = a5.data.create_input("plasma_1D", dryrun=True)
    ndict   = a5.data.create_input("N0_1D", dryrun=True)
    bzrdict = a5.data.create_input("Boozer", dryrun=True)
    mhddict = a5.data.create_input("MHD_STAT", dryrun=True)
    asigmadict = a5.data.create_input("asigma_loc", dryrun=True)
    
    a5.simulation_initinputs(bfield=bdict, efield=edict, plasma=pdict, neutral=ndict,
                             wall=wdict3, boozer=bzrdict, mhd=mhddict, asigma=asigmadict)


    generate_markers = False

    if generate_markers:
        print('generating markers')
        nmarkers = 10
        mrk = Marker.generate("gc", n=nmarkers, species="alpha")
        mrk["energy"][:] = 3.5e6
        mrk["pitch"][:]  = 0.99 - 1.98 * np.random.rand(nmarkers,)
        mrk["r"][:]      = np.linspace(3.5, 5.4, nmarkers)
        a5.data.create_input("gc", **mrk)

        #print('Reading markers')
        #mrk = a5.data.marker.active.read()
        #print(mrk.keys())
        print('providing markers')
        a5.simulation_initmarkers(**mrk)

    else:
        a5.simulation_initmarkers( **mrkrdict2)




    print('generating options')
    #a5.data.create_input("options tutorial")
    #opt = OptionTemplates().options_tutorial()[1]

    #print('Reading options')
    #opt = a5.data.options.active.read()
    print('providing options')
    opt = options.Opt.convert_xml(optionsxml)
    a5.simulation_initoptions(**opt)



    # Reduce the simulation time a bit:
    max_simtime=1.0e-4
    print("Reducing max_simtime to {}s.".format(max_simtime))
    #M.sim.endcond_max_simtime = max_simtime
    a5._sim.endcond_lim_simtime=max_simtime


    #print("Injecting 2D wall from imas")
    #M.inject_wall_2d( wdict['r'], wdict['z']  )

    #print("Injecting 3D wall from imas")
    #a5.provide_wall_3d( wdict3['x1x2x3'], wdict3['y1y2y3'], wdict3['z1z2z3']  )

    #print("Injecting stellarator B-field from imas")
    #a5.provide_BSTS(bdict)


    print('starting simulation')
    vrun_output = a5.simulation_run()





    #print("offloading")
    #M.offload()

    #print("simulating")
    #M.run_simulation()

    #print("gathering output")
    #M.gather_output()

    #print("printing marker summary")
    #M.print_marker_summary()

    #print("writing output h5")
    #M.write_output_h5()


    # This could be done only for RANK=0

    print("saving output ids distributions")
    dist = a5py.ascot5io.imas.distributions()

    # Add a reference to distributions into the dist object.
    # (dist.ids=distributions)
    dist.setIds(distributions)

    # Fill in the data from ASCOT internal datastructures into dist.ids
    dist.fill(runobject=vrun_output, metadata={} )




    print("freeing memory")
    #M.free_c()
    a5.simulation_free(diagnostics=True)

    #print("finalizing MPI")
    #M.finalize()


    print('Finished test')


    #print('Stopping test')
    #quit()
