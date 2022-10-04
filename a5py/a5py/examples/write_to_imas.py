# coding: utf-8

# get helloworld.h5 with
# python ../preprocessing/simpleruns.py 

from a5py.ascot5io.ascot5 import Ascot
import a5py.ascot5io.imas
import a5py.ascot5io.wall_2D
import a5py.ascot5io.wall_3D

a=Ascot('helloworld.h5')
wdict3 = a.wall.active.read()

w3dimas=a5py.ascot5io.imas.wall_3d()
w3dimas.write(wdict3,'akaslos','mywall','3',101,101,{})


