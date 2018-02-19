import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer               

class ascotpy:
    ascotlib = None
    ascotfn  = None
    bfield_initialized  = False
    efield_initialized  = False
    plasma_initialized  = False
    wall_initialized    = False

    def __init__(self, libpath, ascotfn):
        self.ascotlib = ctypes.CDLL(libpath)
        self.ascotfn = ascotfn.encode('UTF-8')

        # Since Python has weak typing while C has strong, we need to define #
        # the type of the arguments for all python wrapper functions.        #
        fun = self.ascotlib.ascotpy_bfield_eval_B
        fun.restype = None
        fun.argtypes = [ctypes.c_int,
                        ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                        ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]


    def ascotpy_init_bfield(self):
        self.ascotlib.ascotpy_init_bfield(self.ascotfn)
        self.bfield_initialized  = True

    def ascotpy_init_efield(self):
        self.ascotlib.ascotpy_init_efield(self.ascotfn)
        self.efield_initialized  = True

    def ascotpy_init_plasma(self):
        self.ascotlib.ascotpy_init_plasma(self.ascotfn)
        self.plasma_initialized  = True

    def ascotpy_init_wall(self):
        self.ascotlib.ascotpy_init_wall(self.ascotfn)
        self.wall_initialized  = True

    def ascotpy_eval_bfield(self, R, phi, z, vecB=None, jacB=None, psiB = None, axis=None):
        if(not self.bfield_initialized):
            return 0
        if(vecB is not None):
            
            BR   = np.zeros((len(R),1))
            Bphi = np.zeros((len(R),1))
            Bz   = np.zeros((len(R),1))

            self.ascotlib.ascotpy_bfield_eval_B(len(R),R,phi,z,BR,Bphi,Bz)
            vecB = np.squeeze(np.array([BR, Bphi, Bz]))
            
            return vecB
        

    def ascotpy_eval_efield(self, R, phi, z, vecE=None):
        if( (not self.bfield_initialized) or (not self.efield_initialized) ):
            return

    def ascotpy_eval_plasma(self, R, phi, z, species=None, temperatur=None, density=None, Anum=None, Znum=None):
        if( (not self.bfield_initialized) or (not self.plasma_initialized) ):
            return

    def ascotpy_eval_FOccol(self, R, phi, z, species=None, temperatur=None, density=None, Anum=None, Znum=None):
        return

    def ascotpy_eval_GCccol(self, R, phi, z, species=None, temperatur=None, density=None, Anum=None, Znum=None):
        return
        

ascot = ascotpy("/l/sarkimk1/repos/ascot5/ascotpy.so","ascot.h5")
ascot.ascotpy_init_bfield()

R = np.array([6.2, 7, 8],dtype=np.float64)
phi = np.array([0, 0, 0],dtype=np.float64)
z = np.array([0.2, 0.2, 0.2],dtype=np.float64)
vecB=1
vecB = ascot.ascotpy_eval_bfield(R, phi, z, vecB=vecB)
print(vecB)


