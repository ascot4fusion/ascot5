import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pyvista as pv

from mpl_toolkits.axes_grid1 import make_axes_locatable
from .plothelpers import openfigureifnoaxes

@openfigureifnoaxes(projection=None)
def still(wallmesh, points=None, data=None, log=False, cpos=None, cfoc=None,
          cang=None, axes=None, cax=None):
    """
    Take a still shot of the mesh and display it using matplotlib backend.

    The rendering is done using vtk but the vtk (interactive) window is not
    displayed. It is recommended to use the interactive plot to find desired
    camera position and produce the actual plot using this method. The plot
    is shown using imshow and acts as a regular matplotlib plot.

    Args:
        wallmesh : Polydata <br>
            Mesh representing the wall.
        points : array_like, optional <br>
            Array Npoint x 3 defining points (markers) to be shown. For
            each point [x, y, z] coordinates are given.
        cpos : array_like, optional <br>
            Camera position coordinates [x, y, z].
        cfoc : array_like, optional <br>
            Camera focal point coordinates [x, y, z].
        cang : array_like, optional <br>
            Camera angle [azimuth, elevation, roll].
        axes : Axes, optional <br>
            The Axes object to draw on.
        cax : Axes, optional <br>
            The Axes object for the color data (if c contains data), otherwise
            taken from axes.
    """
    p = pv.Plotter(off_screen=True)
    if data is None:
        p.add_mesh(wallmesh, color=[0.9,0.9,0.9])
    else:
        cmap = mpl.colormaps["Reds"].copy()
        cmap.set_bad(color=[0.9,0.9,0.9])
        p.add_mesh(wallmesh, scalars=data, cmap=cmap, log_scale=log)
        p.remove_scalar_bar()

        maxval = np.nanmax(wallmesh.cell_data[data])
        minval = np.nanmin(wallmesh.cell_data[data])

    if points is not None:
        p.theme.color = 'black'
        p.add_points(points, render_points_as_spheres=True, point_size=10)

    # Set camera
    if cpos is not None:
        p.camera.position = cpos
    if cfoc is not None:
        p.camera.focal_point = cfoc
    if cang is not None:
        p.camera.azimuth   = cang[0]
        p.camera.elevation = cang[1]
        p.camera.roll      = cang[2]

    p.show()
    axes.imshow(p.image)
    axes.set_xticks([])
    axes.set_yticks([])

    if data is not None:
        if log:
            norm = mpl.colors.LogNorm(vmin=minval, vmax=maxval)
        else:
            norm = mpl.colors.Normalize(vmin=0, vmax=maxval)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        cbar = plt.colorbar(sm, ax=axes, cax=cax)
        #cbar.set_label(r"Energy load J/m$^2$")


def interactive(wallmesh, *args, points=None, data=None, log=False, cpos=None,
                cfoc=None, cang=None):
    """
    Open vtk window to display interactive view of the wall mesh.

    Args:
        wallmesh : Polydata <br>
            Mesh representing the wall.
        *args : tuple (str, method), optional <br>
            Key (str) method pairs. When key is pressed when the plot is
            displayed, the associated method is called. The method should take
            Plotter instance as an argument.
        points : array_like, optional <br>
            Array Npoint x 3 defining points (markers) to be shown. For
            each point [x, y, z] coordinates are given.
        cpos : array_like, optional <br>
            Camera position coordinates [x, y, z].
        cfoc : array_like, optional <br>
            Camera focal point coordinates [x, y, z].
        cang : array_like, optional <br>
            Camera angle [azimuth, elevation, roll].
    """
    p = pv.Plotter()

    if data is None:
        p.add_mesh(wallmesh, color=[0.9,0.9,0.9], log_scale=log)
    else:
        cmap = mpl.colormaps["Reds"].copy()
        cmap.set_bad(color=[0.9,0.9,0.9])
        p.add_mesh(wallmesh, scalars=data, cmap=cmap, log_scale=log)

    if points is not None:
        p.theme.color = 'black'
        p.add_points(points, render_points_as_spheres=True, point_size=10)

    # Set events
    for i in range(len(args)):
        def wrapper(*wargs):
            args[i][1](p)

        p.add_key_event(args[i][0], wrapper)

    # Set camera
    if cpos is not None:
        p.camera.position = cpos
    if cfoc is not None:
        p.camera.focal_point = cfoc
    if cang is not None:
        p.camera.azimuth   = cang[0]
        p.camera.elevation = cang[1]
        p.camera.roll      = cang[2]

    p.show()


@openfigureifnoaxes(projection=None)
def loadvsarea(wetted, loads, axes=None):
    """
    """
    idx = np.argsort(-loads)
    wetted = np.cumsum(wetted[idx])
    loads  = loads[idx]

    axes.set_xscale('log')
    axes.set_yscale('linear')
    #axes.set_ylim((1e4, 2e8))
    #ax.set_xlim((1,14))
    #ax.spines['left'].set_visible(False)
    #ax.yaxis.set_ticks_position('right')
    #ax.yaxis.set_visible(False)

    #divider = make_axes_locatable(axes)
    #axes1 = divider.append_axes("left", size=1.8, pad=0, sharex=axes)
    #axes1.set_yscale('log')
    #ax1.set_xscale('linear')
    #ax1.set_xlim((1e-4, 9.99e-1))
    #ax1.spines['right'].set_visible(False)
    #ax1.yaxis.set_ticks_position('left')
    #plt.setp(ax1.get_xticklabels(), visible=True)
    axes.plot(loads, wetted)
    axes.set_xlabel(r"Wet area [m$^2$]")
    axes.set_ylabel(r"Wet area [m$^2$]")


def defaultcamera(wallmesh):
    """
    Get default camera (helper function for the 3D plots).

    Default camera is located at R = (Rmax + Rmin) / 2, phi = 0,
    z = (zmax, zmin) / 2, where the min/max values are taken from the bounding
    box of the wall mesh. The focal point is at same (R,z) but phi = 10. The
    camera angle is zero in all parameters.

    Args:
        wallmesh : Polydata <br>
            Mesh representing the wall.
    Returns:
        cpos : array_like, optional <br>
            Camera position coordinates [x, y, z].
        cfoc : array_like, optional <br>
            Camera focal point coordinates [x, y, z].
        cang : array_like, optional <br>
            Camera angle [azimuth, elevation, roll].
    """
    # Find min/max values
    Rmax = np.sqrt( np.amax(wallmesh.points[:,0]**2 + wallmesh.points[:,1]**2) )
    Rmin = np.sqrt( np.amin(wallmesh.points[:,0]**2 + wallmesh.points[:,1]**2) )
    zmax = np.amax( wallmesh.points[:,2] )
    zmin = np.amin( wallmesh.points[:,2] )

    # Camera position in (R,phi,z)
    cpos = np.array([( Rmax + Rmin ) / 2,  0, ( zmax + zmin ) / 2])
    cfoc = np.array([( Rmax + Rmin ) / 2, 10, ( zmax + zmin ) / 2])
    cang = np.array([0, 0, 0])

    # Coordinate transformation (R,phi,z) -> (x,y,z)
    cpos = np.array([cpos[0] * np.cos( np.pi * cpos[1] / 180 ),
                     cpos[0] * np.sin( np.pi * cpos[1] / 180 ),
                     cpos[2]])
    cfoc = np.array([cfoc[0] * np.cos( np.pi * cfoc[1] / 180 ),
                     cfoc[0] * np.sin( np.pi * cfoc[1] / 180 ),
                     cfoc[2]])

    return (cpos, cfoc, cang)
