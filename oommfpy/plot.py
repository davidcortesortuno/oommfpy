from . import OOMMFData
import colorsys
import numpy as np
import click
import scipy.spatial as ss
import matplotlib.pyplot as plt


# -----------------------------------------------------------------------------
# Utilities to generate a HSL colourmap from the magnetisation field data

def convert_to_RGB(hls_color):
    return np.array(colorsys.hls_to_rgb(hls_color[0] / (2 * np.pi),
                                        hls_color[1],
                                        hls_color[2]))


def generate_colours(field_data, colour_model='rgb'):
    """
    field_data      ::  (n, 3) array
    """
    hls = np.ones_like(field_data)
    hls[:, 0] = np.arctan2(field_data[:, 1],
                           field_data[:, 0]
                           )
    hls[:, 0][hls[:, 0] < 0] = hls[:, 0][hls[:, 0] < 0] + 2 * np.pi
    hls[:, 1] = 0.5 * (field_data[:, 2] + 1)

    if colour_model == 'rgb':
        rgbs = np.apply_along_axis(convert_to_RGB, 1, hls)
        return rgbs

    elif colour_model == 'hls':
        return hls

    else:
        raise Exception('Specify a valid colour model: rgb or hls')


# -----------------------------------------------------------------------------


@click.command()
@click.option('-i', '--input_omf_file', type=str,
              help='Path to OMF file', required=True)
@click.option('--z', type=float, help='z coordinate')
def plot_omf(input_omf_file, z=False):
    data = OOMMFData(input_omf_file)
    data.read_m()
    data.set_coordinates()

    if not z:
        z = data.zs.min()

    dz = data.dz * 1e9  # To get in nm unit. We might update this in the future
    if z > data.zs.max() + 0.5 * dz or z < data.zs.min() - 0.5 * dz:
        raise Exception("z outside range of data")
    else:
        z_tree = ss.cKDTree(data.zs.reshape(len(data.zs), 1))
        z_data = data.zs[z_tree.query([z])[1]]

    f = plt.figure()
    ax = f.add_subplot(111)

    z_fltr = data.z == z_data
    Ms_fltr = ((data.mx[z_fltr] ** 2 +
                data.my[z_fltr] ** 2 +
                data.mz[z_fltr] ** 2) < 1e-6
               )
    # ax.scatter(data.x[z_fltr][Ms_fltr],
    #            data.y[z_fltr][Ms_fltr],
    #            c=data.mz[z_fltr][Ms_fltr],
    #            cmap='RdYlBu', vmin=-1, vmax=1,
    #            marker='s'
    #            )
    mz_data = data.mz[z_fltr]
    mz_data[Ms_fltr] = np.nan
    mz_data.shape = (-1, len(data.xs))
    ax.imshow(mz_data,
              cmap='RdYlBu', vmin=-1, vmax=1,
              extent=[data.xs.min() - 0.5 * data.dx * 1e9,
                      data.xs.max() + 0.5 * data.dx * 1e9,
                      data.ys.min() - 0.5 * data.dy * 1e9,
                      data.ys.max() + 0.5 * data.dy * 1e9
                      ]
              )
    plt.show()
