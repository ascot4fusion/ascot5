import numpy as np
import subprocess
import matplotlib.pyplot as plt
import a5py.ascot5io.ascot5 as ascot5
import a5py.ascot5io.options as options
import a5py.ascot5io.B_GS as B_GS
import a5py.ascot5io.markers as markers
from testcase import createbase

def run():

    # Test options
    Rmin   = 4
    Rmax   = 8
    nR     = 80
    zmin   = -5
    zmax   = 5
    nz     = 160
    phimin = 0
    phimax = 360
    nphi   = 360

    mrkNR  = 3
    mrkNz  = 3
    mrkphi = np.array([0, 180.25, 359.5])
    
    # Init testbase.
    fn = "test1.h5"
    Bxyz = np.array([1, 0, 0])
    Exyz = np.array([0, 0, 0])
    n = 1e20
    T = 1e3
    createbase(fn, Bxyz, Exyz, n, T)
    
    # We are simulating field lines for 0 seconds.
    # Only inistate is of interest.
    o = options.read_hdf5(fn)
    o["SIM_MODE"]                = 0*o["SIM_MODE"] + 4
    o["ENDCOND_SIMTIMELIM"]      = 0*o["ENDCOND_SIMTIMELIM"] + 1
    o["ENDCOND_MAX_SIM_TIME"]    = 0*o["ENDCOND_MAX_SIM_TIME"]
    o["ENABLE_ORBIT_FOLLOWING"]  = 0*o["ENABLE_ORBIT_FOLLOWING"] + 1
    options.write_hdf5(fn,o)

    # Clone base into four copies.
    fn1 = fn
    fn2 = "test2.h5"
    fn3 = "test3.h5"
    fn4 = "test4.h5"
    subprocess.call(["cp",fn,fn2])
    subprocess.call(["cp",fn,fn3])
    subprocess.call(["cp",fn,fn4])

    # Init magnetic field.
    R0 = 6.2
    z0 = 0
    B_phi0 = 5.3
    psi_mult = 200
    psi_coeff = np.array([8.629491085780416348e-02, 3.279306587723925803e-01, 5.268677701240817024e-01, -2.366208946912087274e-01, 
                          3.825826765593096646e-01, -3.573153147754407621e-01, -1.484166833037287025e-02, 1.506045943286430100e-01, 
                          7.428226459414810634e-01, -4.447153105104519888e-01, -1.084640395736786167e-01, 1.281599235951017685e-02, 
                          -0.155])
    Nripple = 18
    a0 = 2
    alpha0 = 2
    delta0 = 0.1

    B_GS.write_hdf5(fn1, R0, z0, B_phi0, psi_mult, psi_coeff)

    B_GS.write_hdf5_B_2D(fn2, R0, z0, B_phi0, psi_mult, psi_coeff, 
                         Rmin, Rmax, nR, zmin, zmax, nz)

    B_GS.write_hdf5(fn3, R0, z0, B_phi0, psi_mult, psi_coeff, 
                    Nripple=Nripple, a0=a0, alpha0=alpha0, delta0=delta0)

    B_GS.write_hdf5_B_3D(fn4, R0, z0, B_phi0, psi_mult, psi_coeff, 
                         Nripple, a0, alpha0, delta0, 
                         Rmin, Rmax, nR, zmin, zmax, nz, phimin, phimax, nphi)

    # Init markers.
    Rg = np.linspace(Rmin,Rmax,nR*mrkNR)
    Rg = Rg[1:-1]
    zg = np.linspace(zmin,zmax,nz*mrkNz)
    zg = zg[1:-1]
    R, z = np.meshgrid(Rg,zg)
    R = R.flatten()
    z = z.flatten()

    phi    = 0*np.ones(R.shape)
    pitch  = np.ones(R.shape)
    weight = np.ones(R.shape)
    time   = 0*np.ones(R.shape)
    
    N   = R.size
    ids = np.linspace(1,N,N)

    # For 2D fields.
    markers.write_hdf5_fieldlines(fn1, N, ids, R, phi, z, pitch, weight, time)
    markers.write_hdf5_fieldlines(fn2, N, ids, R, phi, z, pitch, weight, time)

    # For 3D fields.
    R, z = np.meshgrid(Rg,zg)
    R = np.tile(R,(1,1,3))
    z = np.tile(z,(1,1,3))
    phi = np.tile(phi,(1,3))
    phi[0:N]   = mrkphi[0];
    phi[N:2*N] = mrkphi[1];
    phi[2*N:]  = mrkphi[2];
    N = 3*N

    R = R.flatten()
    z = z.flatten()
    phi = phi.flatten()
    pitch  = np.ones(R.shape)
    weight = np.ones(R.shape)
    time   = 0*np.ones(R.shape)

    ids = np.linspace(1,N,N)

    markers.write_hdf5_fieldlines(fn3, N, ids, R, phi, z, pitch, weight, time)
    markers.write_hdf5_fieldlines(fn4, N, ids, R, phi, z, pitch, weight, time)

    # Simulate.
    subprocess.call(["./ascot5_main", "--in="+fn1[0:-3]])
    subprocess.call(["./ascot5_main", "--in="+fn2[0:-3]])
    subprocess.call(["./ascot5_main", "--in="+fn3[0:-3]])
    subprocess.call(["./ascot5_main", "--in="+fn4[0:-3]])

    # Read inistate psi, magnetic field, and gradients.
    is1 = ascot5.read_hdf5(fn1,"states")["states"]["inistate"]
    is2 = ascot5.read_hdf5(fn2,"states")["states"]["inistate"]
    is3 = ascot5.read_hdf5(fn3,"states")["states"]["inistate"]
    is4 = ascot5.read_hdf5(fn4,"states")["states"]["inistate"]

    
    R = np.unique(is1["R"])
    z = np.unique(is1["z"])
    
    B0 = np.reshape(is1["rho"],(z.size,R.size)) - np.reshape(is2["rho"],(z.size,R.size))
    
    # Plot if needed.
    plt.contourf(R,z,B0)
    plt.show()

    # Compare.


    # Clean.
    subprocess.call(["rm", fn1])
    subprocess.call(["rm", fn2])
    subprocess.call(["rm", fn3])
    subprocess.call(["rm", fn4])

if __name__ == '__main__':
    run()
