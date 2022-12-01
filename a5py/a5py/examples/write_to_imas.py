# coding: utf-8

# get helloworld.h5 with
# python ../preprocessing/simpleruns.py 

from a5py.ascot5io.ascot5 import Ascot
import a5py.ascot5io.imas
import a5py.ascot5io.wall_2D
import a5py.ascot5io.wall_3D

filename='helloworld.h5'
a=Ascot(filename)
wdict3 = a.wall.active.read()
bsts   = a.bfield.active.read()

w3dimas=a5py.ascot5io.imas.wall_3d()
bstsimas=a5py.ascot5io.imas.B_STS()
w3dimas.write( wdict3,'akaslos','mywall', '3',101,101, {'comment':filename} )
bstsimas.write(  bsts,'akaslos','ggdtest','3',32,   3, {'comment':filename} )


