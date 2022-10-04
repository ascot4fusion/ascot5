# coding: utf-8

# For Simppa at least, following recipe worked to get the below working:
# 1. in ascot5.yml, set the python version: python=3.6, and create the environment
# 2. module load imasenv
# 3. conda activate ascot5
# the below works.


import a5py.ascot5io.imas
import a5py.ascot5io.wall_2D
import a5py.ascot5io.wall_3D
w2d=a5py.ascot5io.imas.wall_2d()
w3d=a5py.ascot5io.imas.wall_3d()
# Gateway
# wdict=w2d.read("g2jvarje","test","3",92436,272)

# ITER SDCC
wdict3=w3d.read("akaslos","mywall","3",11,11)

wdict=w2d.read("akaslos","test","3",92436,272)

a5py.ascot5io.wall_3D.write_hdf5('from_imas.h5',**wdict3)
a5py.ascot5io.wall_2D.write_hdf5('from_imas.h5',**wdict)


