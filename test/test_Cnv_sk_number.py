# import numpy as np
import oommfpy.oommfpy as op
import re
import subprocess
import shutil
from pathlib import Path
import pytest
this_dir = Path(__file__).resolve().parent


def generate_omfs():

    """
    Simulate a skyrmion reducing the mesh discretisation
    using the OOMMF script: isolated_sk_DMI_Cnv.mif

    It is necessary to have OOMMF installed with the latest DMI
    modules. Recommended to install JOOMMF's OOMMF conda package
    """

    OMF_DIR = Path(this_dir / 'omfs/')
    if OMF_DIR.exists():
        shutil.rmtree(OMF_DIR)
    OMF_DIR.mkdir(exist_ok=True)

    for n in range(20, 101, 20):
        SIM_NAME = 'omfs/isolated_sk_Cnv_n_{:03d}'.format(n)
        SCRIPT = Path(this_dir / 'isolated_sk_DMI_Cnv.mif')

        job = ('oommf boxsi -threads 2 '
               '-parameters "NX {0} '
               'BASENAME {1}" '
               '"{2}"'.format(n, SIM_NAME, SCRIPT)
               )
        print(job)

        subprocess.call(job, shell=True)


# -----------------------------------------------------------------------------

# Function to sort OMF files according to number of mesh sites n from filename
def get_n(f):
    return int(re.search(r'(?<=n_)\d+(?=-Oxs)', str(f)).group(0))


@pytest.mark.skip(reason="Requires better comparison / asserts")
def test_sk_number_vs_mesh_discretisation():
    """
    Compute the sk number at the bottom of a disk with a skyrmion
    (the function computes the sk number in the xy plane). The system has a
    DMI with Cnv symmetry (interfacial)

    The test repeats the calculation for samples with increasing mesh
    resolutions (smaller mesh spacings) and prints the result for every mesh
    discretisation. Sk number values should converge to 1 for finer meshes

    Samples with different mesh sizes are computed through the
    isolated_sk_DMI_Cnv.mif OOMMF script and are called from the
    generate_omfs() function
    """

    generate_omfs()

    # The files in the right order
    _files = Path(this_dir / 'omfs/').glob('*.omf')
    _files = sorted(_files, key=get_n)

    mesh_length = 60  # nm

    print('Sk Number as a function of discretisation')
    print('FD: finite differences     SL: spin lattice')
    for FILE in _files:
        n = get_n(FILE)
        oommf_data = op.MagnetisationData(FILE)
        oommf_data.generate_field()

        print('dx = {:.2f} nm --> Q_FD = {}'.format(
            mesh_length / n,
            oommf_data.compute_sk_number(plane='xy', index=0,
                                         method='finite_differences')))

        print('dx = {:.2f} nm --> Q_SL = {}'.format(
            mesh_length / n,
            oommf_data.compute_sk_number(plane='xy', index=0,
                                         method='spin_lattice')))
        print('---')


if __name__ == '__main__':
    test_sk_number_vs_mesh_discretisation()
