from .. import MagnetisationData
from . import plot_tools
import numpy as np
import click
import scipy.spatial as ss
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
# import matplotlib
# from matplotlib import rcParams
from mpl_toolkits.axes_grid1 import make_axes_locatable


# -----------------------------------------------------------------------------


def plot_omf(omf_file, ax=None,
             savefig=None, cmap='hls', cbar=True, dpi=150,
             cbar_offsets=[-0.15, -0.2], cbar_size=0.2):
    """
    A simple function to plot a slice of the system in the xy-plane with a
    specific colormap (HLS by default)

    omf_file        :: Path to omf file
    ax              :: Matplotlib axes object where to plot if specified
    savefig         :: Filename of a figure to be saved with a format allowed
                       by matplotlib
    cmap            :; hls or any colormap accepted by matplotlib
    cbar            :: show a wheel colorbar for the hls cmap or a linear
                       colorbar if other cmap is used
    dpi             :: figure dpi
    cbar_offsets    :: offset x-y positions for the colorbar (hls only)
    cbar_size       :: size relative to plot axes fraction (only for hls color)
    """

    data = MagnetisationData(omf_file)
    data.generate_field()
    data.generate_coordinates()

    if not ax:
        f = plt.figure()
        ax = f.add_subplot(111)

    if cmap == 'hls':
        spin_data = plot_tools.generate_colours(data.field[:, :])
        spin_data[data.field_norm < 1e-10] = [1., 1., 1.]
        spin_data = spin_data.reshape(-1, data.nx, 3)

        p = ax.imshow(spin_data, origin='lower', interpolation='None',
                      vmin=0, vmax=2 * np.pi,
                      extent=[data.xmin * 1e9, data.xmax * 1e9,
                              data.ymin * 1e9, data.ymax * 1e9]
                      )
        if cbar:
            box = ax.get_position()
            axColor = plt.axes([box.x1 + cbar_offsets[0],
                                box.y1 + cbar_offsets[1],
                                cbar_size, cbar_size], projection='polar')
            daz = 361  # Number of discretised values i.e. N of colours
            azimuths = np.linspace(0, 360, daz)
            dzn = 30   # Zeniths or "rings" from the centre towards the bndry
            zeniths = np.arange(0, dzn, 1)

            rgb = np.ones((dzn * daz, 3))
            # Set the HLS hue value from 0 to 2 PI from the azimuth values
            # We tile the circle "dz" times:
            #   [0 ... 2PI] -> [0...2PI 0 .. 2PI ...]
            rgb[:, 0] = np.tile(np.radians(azimuths), dzn)
            # For every circle (daz values) we increase the Light value
            # from 1 to 0, i.e. from white to black, dz times:
            #               |--dzn--|
            #  [1 .. 0] -> [1 1 ... 1 0.8 0.8 ... 0.8 0.6 ... 0 0 ... 0]
            greys = np.zeros(dzn)
            # Last 2 rings are left completely black:
            greys[:-1] = np.linspace(1, 0, dzn - 1)
            rgb[:, 1] = np.repeat(greys, daz)
            # Now we convert every row in HLS to RGB values
            rgb = np.apply_along_axis(plot_tools.convert_to_RGB, 1, rgb)

            # And plot in the polar axes:
            axColor.pcolormesh(np.radians(azimuths), zeniths,
                               # only necessary as required n of args:
                               np.zeros((dzn, daz)),
                               # cmap=plt.cm.hsv
                               color=rgb,
                               shading='auto'
                               )
            axColor.set_yticks([])
            # axColor.set_xticks([0, np.pi * 0.5, np.pi, 1.5 * np.pi])
            axColor.set_thetagrids([0, 90, 180, 270])
            axColor.tick_params(axis='x', pad=0)
            axColor.set_xticklabels([r'$0$', r'$\pi/2$',
                                     r'$\pi$', r'$3\pi/2$'],
                                    # fontsize=18
                                    )
            axColor.text(0.5, 0.5, r'$\vec{m}$',
                         horizontalalignment='center',
                         verticalalignment='center',
                         transform=axColor.transAxes,
                         # fontsize=20
                         )

    else:
        # Spin data in a grid
        spin_z = data.field[:, 2].reshape(-1, data.nx)

        ax.imshow(spin_z, origin='lower', cmap=cmap, interpolation='None',
                  vmin=-1, vmax=1,
                  extent=[data.xmin * 1e9, data.xmax * 1e9,
                          data.ymin * 1e9, data.ymax * 1e9]
                  )

        if cbar:
            plt.colorbar()

    ax.set_ylabel(r'y (nm)')
    ax.set_xlabel(r'x (nm)')

    if savefig:
        plt.savefig(savefig, bbox_inches='tight', dpi=dpi)

    plt.show()


def plot_charge_density(omf_file, ax=None, savefig=None, dpi=150,
                        plane='xy', index=0):
    """
    Testing
    Plot the sk number density
    """

    data = MagnetisationData(omf_file)
    data.generate_field()
    data.generate_coordinates()
    data.compute_sk_number(plane=plane, index=index)

    if not ax:
        f, ax = plt.subplots()

    charge = data.sk_number
    vmax = np.max(np.abs(charge))

    charge.reshape(-1,)[data.field_norm < 1e-10] = np.nan

    ax.imshow(charge, origin='lower', cmap='RdYlBu', interpolation='None',
              vmin=-vmax, vmax=vmax,
              extent=[data.xmin * 1e9, data.xmax * 1e9,
                      data.ymin * 1e9, data.ymax * 1e9]
              )
    ax.set_ylabel(r'y (nm)')
    ax.set_xlabel(r'x (nm)')
    plt.colorbar()

    if savefig:
        plt.savefig(savefig, bbox_inches='tight', dpi=dpi)
