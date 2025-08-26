import colorsys
import hsluv
import numpy as np


# -----------------------------------------------------------------------------
# Utilities to generate a HSL colourmap from the magnetisation field data
#
def convert_HSLUV_to_RGB(hls_color):
    """
    From hsluv lib:
    # hue is a float between 0 and 360, saturation and lightness are floats between 0 and 100. 
    """
    # NOTE: saturation is in the [2] column
    return np.array(hsluv.hsluv_to_rgb([hls_color[0] * 180 / np.pi,
                                        hls_color[2] * 100.,
                                        hls_color[1] * 100.]))


def convert_HLS_to_RGB(hls_color):
    return np.array(colorsys.hls_to_rgb(hls_color[0] / (2 * np.pi),
                                        hls_color[1],
                                        hls_color[2]))


def generate_colours(field_data, colour_model='hls'):
    """
    Convert field data to HLS or HSLUV colours, in RGB format
    This can be used in imshow

    field_data      ::  (n, 3) array
    """
    hls = np.ones_like(field_data)
    hls[:, 0] = np.arctan2(field_data[:, 1], field_data[:, 0])
    hls[:, 0][hls[:, 0] < 0] = hls[:, 0][hls[:, 0] < 0] + 2 * np.pi
    hls[:, 1] = 0.5 * (field_data[:, 2] + 1)

    if colour_model == 'hls':
        rgbs = np.apply_along_axis(convert_HLS_to_RGB, 1, hls)
        return rgbs

    elif colour_model == 'hsluv':
        rgbs = np.apply_along_axis(convert_HSLUV_to_RGB, 1, hls)
        return rgbs

    else:
        raise Exception('Specify a valid colour model: hls or hsluv')


# -----------------------------------------------------------------------------
