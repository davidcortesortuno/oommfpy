import oommfpy
import numpy as np
import glob


# -----------------------------------------------------------------------------

# Expected Values from the OOMMF simulation
coords_ref = np.array([[0.5, 0.5, 0.5],
                       [1.5, 0.5, 0.5],
                       [2.5, 0.5, 0.5],
                       [0.5, 1.5, 0.5],
                       [1.5, 1.5, 0.5],
                       [2.5, 1.5, 0.5],
                       [0.5, 2.5, 0.5],
                       [1.5, 2.5, 0.5],
                       [2.5, 2.5, 0.5],
                       ])

spins_ref = np.array([[0.0, 0.0, 1.0],
                      [0.0, 0.0, 1.0],
                      [0.0, 0.0, 1.0],
                      [0.0, 1.0, 0.0],
                      [0.0, 1.0, 0.0],
                      [0.0, 1.0, 0.0],
                      [1.0, 0.0, 0.0],
                      [1.0, 0.0, 0.0],
                      [1.0, 0.0, 0.0],
                      ])


# -----------------------------------------------------------------------------


def test_all_formats():
    for omf_file in glob.glob('test_*.omf'):
        print(omf_file)

        data = oommfpy.OOMMFData(omf_file)
        # print(data.data_format)
        data.read_m()

        # Numpy assertion
        np.testing.assert_array_equal(data.mx, spins_ref[:, 0])
        np.testing.assert_array_equal(data.my, spins_ref[:, 1])
        np.testing.assert_array_equal(data.mz, spins_ref[:, 2])

        # The base values for the coordinates, in nm units
        # It seems irreg binary4 format is the less accurate
        assert np.abs(data.xbase * 1e9 - 0.5) < 1e-7
        assert np.abs(data.ybase * 1e9 - 0.5) < 1e-7
        assert np.abs(data.zbase * 1e9 - 0.5) < 1e-7


if __name__ == "__main__":
    test_all_formats()
