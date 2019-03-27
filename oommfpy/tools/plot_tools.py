import colorsys
import numpy as np


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
