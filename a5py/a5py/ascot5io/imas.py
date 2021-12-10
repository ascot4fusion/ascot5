import imas

class a5imas:


    def __init__(self):
        self.ids_name = "ids" #This should be overwritten by dervied classes

    def open(self, user, tokamak, version, shot, run, occurrence=0):

        # ./idsdump.py g2jvarje test 3 92436 0272 distributions
        # user = "g2jvarje"
        # tokamak = "test"
        # version = "3"
        # shot = 92436
        # run = 272
        # ids = "wall"



        self.ids = imas.ids(shot, run)
        self.ids.open_env(user, tokamak, version)

        #ids = ids.split('/')
        #if len(ids) == 1:
        #    occurrence = 0
        #else:
        #    occurrence = int(ids[1])
        #ids = ids[0]

        ids = self.ids_name

        idsdata = self.ids.__dict__[self.ids_name]
        if 'get' not in dir(idsdata):
            idsdata = self.ids.__dict__[self.ids_name + 'Array']
            idsdata.get(occurrence)
            #if idsdata.array:
            #    for slice in idsdata.array:
            #        print(slice)
        else:
            idsdata.get(occurrence)
            #print(idsdata)

        return self.ids

    def close(self):
        self.ids.close()


class wall_2d(a5imas):

    def __init__(self):
        super().__init__()
        self.ids_name = "wall"


    def read(self, user, tokamak, version, shot, run, occurrence=0 ):
        """
        Read an IMAS 2D-wall into a dictionary, that is a drop-in replacement for a dict read from hdf5.

        The data is read from
              wall.description_2d[0].limiter.unit[0].outline.r
        and
              wall.description_2d[0].limiter.unit[0].outline.z

        """

        timeIndex = 0
        unit      = 0

        itm = self.open( user, tokamak, version, shot, run, occurrence )

        #itm.wall.get()
        r = self.ids.wall.description_2d[timeIndex].limiter.unit[unit].outline.r
        z = self.ids.wall.description_2d[timeIndex].limiter.unit[unit].outline.z

        self.close()

        w = {
            "r"         : r,
            "z"         : z,
            "nelements" : len(r),
        }

        return w
