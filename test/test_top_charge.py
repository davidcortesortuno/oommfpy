import oommfpy as omp
import numpy as np
from pathlib import Path


def test_sknum_spin_lattice():
    """
    In this test we have five spins:

            SPINS:                       REFERENCE:

               |  | |                      ^ +y
               |  v | <--                  |
           ----|----|----                  |
            ⦻  |--> |  ⦿                   ⦿-----> +x
               |    |                      +z

    And we need to compute the skyrmion number using the spin lattice method
    For the lower left site we have 4 triangles in counter-clockwise order:
       (+x +z -y)  (-y -z +x)  (+x 0 0)  (+x 0 0)
    The zeros appear at the boundary (no material)
    """
    parent = Path(__file__).resolve().parent

    data = omp.MagnetisationData(parent / './test_top_charge.omf')
    data.generate_field()
    sk_num = data.compute_sk_number(method='spin_lattice')
    print(data.sk_number)

    # We will compare this with the semi-analytical results

    # For the spin next to the lower left corner the triangles that produces a
    # non-zero charge is the first and second triangle: (+x +z -y)  (-y -z +x)
    # The spin lattice formula is:
    #       q_ijk = arctan2(2 * mi * (mj X mk) / (1 + mi * mj + ...))
    # Because the 1st triangle does not have a zero spin in the opposite
    # diagonal corner, its weight is 0.5
    denom = 1.0  # the dot products all cancel to zero
    num = 0.5 * 2 * np.arctan(1.)
    expected_sknum_low_left = num / denom
    # The 2nd triangle does have a zero spin in the opposite
    # diagonal corner, so its weight is 1.0
    denom = 1.0  # the dot products all cancel to zero
    num = 1.0 * 2 * np.arctan(1.)
    expected_sknum_low_left += (num / denom)

    assert abs(data.sk_number[0, 1] - expected_sknum_low_left) < 1e-8

    # we repeat the same for the other triangles, most of them give a charge
    # of zero because the cross products are zero. The other non-zero is the
    # one at the top right with the triangle: (-x -y +z)
    # such that mi * mj X mk = -x * (-x) = 1
    denom = 1.0  # the dot products all cancel to zero
    num = 0.5 * 2 * np.arctan(1.)  # again the weight is 0.5
    expected_sknum_top_right = num / denom
    assert abs(data.sk_number[1, 2] - expected_sknum_top_right) < 1e-8

    # Finally we can check the total charge
    fpi = 4 * np.pi
    exp_tot = (expected_sknum_low_left + expected_sknum_top_right) / fpi
    assert abs(sk_num - exp_tot) < 1e-8
    print(exp_tot)


if __name__ == "__main__":

    test_sknum_spin_lattice()
