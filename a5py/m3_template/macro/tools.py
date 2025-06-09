import imaspy
import libmuscle
from imas_core import imasdef
from imas_core.imasdef import CLOSEST_SAMPLE
from imaspy.ids_defs import MDSPLUS_BACKEND

def send_ids(m3_instance, ids, port_name, t_cur, t_next):

    print('Sends: ', port_name, " : start")
    ids_bytes = ids.serialize(imasdef.DEFAULT_SERIALIZER_PROTOCOL)
    m3_message = libmuscle.Message(t_cur, t_next, ids_bytes)
    m3_instance.send(port_name, m3_message)
    print('Sends: ', port_name, " : end")

def get_ids(ids_name, user, database, shot, run, time):

    print(ids_name, ": READING")

    db_entry = imaspy.DBEntry(MDSPLUS_BACKEND, database, shot, run, user_name=user)
    db_entry.open()
    ids = db_entry.get_slice(ids_name, time, CLOSEST_SAMPLE)
    db_entry.close()

    print(ids_name, ": OK")

    return ids


def get_inputs():

    idsname = 'equilibrium'
    uname = 'g2diy'
    dbname = 'ggdtest'
    shot = 32
    run = 3
    time = 0.0
    equilibrium_3d = get_ids(idsname, uname, dbname, shot, run, time)

    idsname = 'equilibrium'
    uname = 'g2diy'
    dbname = 'test'
    shot = 92436
    run = 306
    time = 0.0
    equilibrium = get_ids(idsname, uname, dbname, shot, run, time)

    idsname = 'core_profiles'
    uname = 'g2diy'
    dbname = 'test'
    shot = 92436
    run = 306
    time = 0.0
    core_profiles = get_ids(idsname, uname, dbname, shot, run, time)

    idsname = 'distribution_sources'
    uname = 'g2diy'
    dbname = 'test'
    shot = 92436
    run = 306
    time = 0.0
    distribution_sources = get_ids(idsname, uname, dbname, shot, run, time)

    idsname = 'wall'
    uname = 'g2diy' 
    dbname = 'test' 
    shot = 92436
    run = 272
    time = 0.0
    wall2d = get_ids(idsname, uname, dbname, shot, run, time)


    idsname = 'wall'
    uname = 'g2diy'
    dbname = 'mywall'
    shot = 201
    run = 101
    time = 0.0
    wall3d = get_ids(idsname, uname, dbname, shot, run, time)

    return equilibrium_3d, equilibrium, core_profiles, distribution_sources, wall2d, wall3d


if __name__ == "__main__":
    equilibrium_3d, equilibrium, core_profiles, distribution_sources, wall2d, wall3d = get_inputs()
    print(equilibrium.metadata.name)
