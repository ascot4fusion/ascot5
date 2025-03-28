from distsource_fun import *
import imas
from imas.imasdef import MDSPLUS_BACKEND
from imas.imasdef import CLOSEST_SAMPLE

#mrkr.read("g2diy","test","3",92436,306)
#w2d.read("g2jvarje","test","3",92436,272)
#w3d.read("g2diy","mywall","3",201,101)
#bsts.read("g2diy","ggdtest","3",32,3)


def get_ids(idsname,user,database,shot,run,time):
    DB = imas.DBEntry(MDSPLUS_BACKEND, database, shot, run, user_name=user)
    DB.open()
    ids = DB.get_slice(idsname, time, CLOSEST_SAMPLE)
    DB.close()
    return ids

idsname='distribution_sources'
uname='g2diy'
dbname='test'
shot=92436
run=306
time=0.0
distr_sour=get_ids(idsname,uname,dbname,shot,run,time)

idsname='wall'
uname='g2jvarje'
dbname='test'
shot=92436
run=272
time=0.0
wall2d=get_ids(idsname,uname,dbname,shot,run,time)


idsname='wall'
uname='g2diy'
dbname='mywall'
shot=201
run=101
time=0.0
wall3d=get_ids(idsname,uname,dbname,shot,run,time)

idsname='equilibrium'
uname='g2diy'
dbname='ggdtest'
shot=32
run=3
time=0.0
equil_b3f=get_ids(idsname,uname,dbname,shot,run,time)

print('after_read_ids')

print(type(distr_sour))
print('number of sources after read', len(distr_sour.source))


distsource_run(distr_sour,wall2d,wall3d,equil_b3f)


