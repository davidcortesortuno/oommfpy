# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.6.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

import oommfpy as op
import numpy as np
import matplotlib.pyplot as plt

# ## Basic functionality

# The most useful method in this library is the `OOMMFData` class, which can read a `omf` file in any format (binary, structured, etc.). From this class, the magnetisation and coordinates can be computed:

data = op.MagnetisationData('../../test/omfs/isolated_sk_Cnv_n_060-Oxs_MinDriver-Magnetization-00-0000498.omf')
data.generate_field()
data.generate_coordinates()

# After reading the coordinates we have the `data.x`, `data.y` and `data.z` arrays with the coordinates. Unique values are stored in the `xs`, `ys` and `zs` arrays:

data.xs

# Similarly for the magnetisation, where the field (spin orientations) is stored in the `(n,3)`-array `data.m` and the magnetisation magnitude at every mesh site is stored in `data.Ms` (or `data.field_norm`):

data.mx

# We can compute the skyrmion number for a given plane and slice:

data.compute_sk_number(plane='xy', index=0)

# The index means computing the slice at `z=`:

data.zs[0]

# We can use this information to generate a plot:

# +
fltr = np.logical_and(data.z == data.zs[0],
                      data.Ms > 1e-4)

f = plt.figure(figsize=(6, 6))
# More effective to use imshow to avoid scaling the
# scattered data point size
plt.scatter(data.x[fltr], data.y[fltr],
            c=data.mz[fltr],
            cmap='RdYlBu',
            marker='s'
            )
plt.show()
# -

# ## Plot Functions

# The library also offers plot tools that can be used to quickly analyse a system:

from oommfpy.tools import plot_omf
from oommfpy.tools.plot_omf import plot_charge_density

# We can also choose not to pass an axes object
f, ax = plt.subplots()
plot_omf('../../test/omfs/isolated_sk_Cnv_n_060-Oxs_MinDriver-Magnetization-00-0000498.omf',
         # cbar=False,  # uncomment if colorbar is not required
         ax=ax,  # comment if not passing an axes object
         cbar_offsets=[-0.25, -0.2], cbar_size=0.2)

# And the charge density:

plot_charge_density('../../test/omfs/isolated_sk_Cnv_n_060-Oxs_MinDriver-Magnetization-00-0000498.omf',
                    plane='xy', index=0)
