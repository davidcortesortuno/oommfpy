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
# rcParams.update({'figure.autolayout': True})


# -----------------------------------------------------------------------------


def plot_omf_slices(input_omf_file, z=False, quiver=False, hls=False):
    """
    Generates an interactive visualisation of the system by showing slices
    in the xy-plane with varying z coordinate
    """
    data = MagnetisationData(input_omf_file)
    data.generate_field()
    data.generate_coordinates()

    if not z:
        z = data.zs.min()

    dz = data.dz * 1e9  # To get in nm unit. We might update this in the future
    if z > data.zs.max() + 0.5 * dz or z < data.zs.min() - 0.5 * dz:
        raise Exception("z outside range of data")
    else:
        z_tree = ss.cKDTree(data.zs.reshape(len(data.zs), 1))
        z_data = data.zs[z_tree.query([z])[1]]

    f, ax = plt.subplots(1, 1,
                         # constrained_layout=True
                         )

    # Filters to plot a slice perp to the z direction
    z_fltr = data.z == z_data
    # Filter to remove points with nzero magnetisation
    Ms_fltr = ((data.mx[z_fltr] ** 2 +
                data.my[z_fltr] ** 2 +
                data.mz[z_fltr] ** 2) < 1e-6
               )

    if not hls:
        m_data = data.field_z[z_fltr]
        m_data[Ms_fltr] = np.nan
        m_data.shape = (-1, len(data.ys))
    else:
        m_data = plot_tools.generate_colours(data.field[z_fltr].reshape(-1, 3))
        m_data[Ms_fltr] = [0.9, 0.9, 0.9]
        m_data.shape = (len(data.xs), len(data.ys), 3)

    # Plot as an image with pixels coloured by mz
    p = ax.imshow(m_data,
                  cmap='RdYlBu', vmin=-1, vmax=1,
                  extent=[data.xs.min() - 0.5 * data.dx * 1e9,
                          data.xs.max() + 0.5 * data.dx * 1e9,
                          data.ys.min() - 0.5 * data.dy * 1e9,
                          data.ys.max() + 0.5 * data.dy * 1e9
                          ],
                  # aspect='auto',
                  origin='lower'
                  )

    # pq = ax.quiver(data.x[z_fltr],
    #                data.y[z_fltr],
    #                data.mx[z_fltr],
    #                data.my[z_fltr],
    #                scale_units='xy', scale=0.3,
    #                color='k'
    #                )

    ax.set_ylabel(r'$y$  (nm)')
    ax.set_xlabel(r'$x$  (nm)')

    # colour bar ..............................................................
    # box = ax.get_position()
    # axCb = plt.axes([box.x1 + 0.01, box.y0, 0.03, box.height])
    # cb = matplotlib.colorbar.ColorbarBase(axCb, 'RdYlBu', orientation="vertical",
    #                                       ticks=[-1, 0, 1],
    #                                       norm=matplotlib.colors.Normalize(vmin=-1,
    #                                                                        vmax=1))

    # Slider ..................................................................
    # Plot a bar to change the z-slice being plotted

    # Create a new set of axes below the plot with the plot width
    divider = make_axes_locatable(ax)
    axSl = divider.append_axes("bottom", "3%", pad="15%")

    # Not working with tight layout:
    # Get the axes positions:
    # box = ax.get_position()
    # axSl = plt.axes([box.x0, box.y0 - 0.12, box.width, 0.02],
    #                 facecolor='lightgoldenrodyellow')

    # Create the Slider widget
    szslice = Slider(axSl, 'z-slice', data.zs.min(), data.zs.max(),
                     valinit=z)

    # This function defines how the data is updated in the plot
    # We find a new z-value and apply the filters again
    def update(z_val):
        z_data = data.zs[z_tree.query([z_val])[1]]
        z_fltr = data.z == z_data
        Ms_fltr = ((data.field_x[z_fltr] ** 2 +
                    data.field_y[z_fltr] ** 2 +
                    data.field_z[z_fltr] ** 2) < 1e-6
                   )
        mz_data = data.field_z[z_fltr]
        mz_data[Ms_fltr] = np.nan
        mz_data.shape = (-1, len(data.xs))
        # Update the plot:
        p.set_data(mz_data)
        f.canvas.draw_idle()

    # Define the update rule for the Slider widget:
    szslice.on_changed(update)

    # .........................................................................

    plt.tight_layout()
    plt.show()


# -----------------------------------------------------------------------------
# CLI functions

@click.command()
@click.option('-i', '--input_omf_file', type=str,
              help='Path to OMF file', required=True)
@click.option('--z', type=float, help='z coordinate')
@click.option('--hls', is_flag=True, help='Plot with HLS colourmap')
def plot_omf_slices_cli(input_omf_file, z, hls):
    plot_omf_slices(input_omf_file, z=z, hls=hls)


# -----------------------------------------------------------------------------


if __name__ == '__main__':
    plot_omf_slices_cli()
