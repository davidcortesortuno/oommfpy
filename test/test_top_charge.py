import oommfpy as omp
import numpy as np

# TODO: check the top charge for a sample with more complex corners

def test_sknum_spin_lattice():
    """
    In this test we have four spins:

     SPINS:                       REFERENCE:

          | |                      ^ +y
          v | <--                  |
        ----|----                  |
        --> |  ⦿                   ⦿-----> +x
            |                      +z

    And we need to compute the skyrmion number using the spin lattice method
    For the lower left site we have 4 triangles in counter-clockwise order:
       (+x +z -y)  (-y 0 +x)  (+x 0 0)  (+x 0 0)
    The zeros appear at the boundary (no material)
    """
    data = omp.MagnetisationData('./test_top_charge.omf')
    data.generate_field()
    sk_num = data.compute_sk_number(method='spin_lattice')

    # We will compare this with the semi-analytical results
    # For the lower left corner the only triangle that produces a non-zero
    # charge is the first triangle: (+x +z -y)
    # The spin lattice formula is:
    #       q_ijk = arctan2(2 * mi * (mj X mk) / (1 + mi * mj + ...))
    # Because this triangle does not have a zero spin in the opposite
    # diagonal corner, its weight is 0.5
    denom = 1.0  # the dot products all cancel to zero
    num = 0.5 * 2 * np.arctan(1.)
    expected_sknum_low_left = num / denom

    assert data.sk_number[0, 0] == expected_sknum_low_left

    # we can repeat the same for the other triangles, two of them give a charge
    # of zero because the cross products are zero. The other non-zero is the
    # one at the top right with the triangle: (-x -y +z)
    # such that mi * mj X mk = -x * (-x) = 1
    denom = 1.0  # the dot products all cancel to zero
    num = 0.5 * 2 * np.arctan(1.)  # again the weight is 0.5
    expected_sknum_top_right = num / denom
    assert data.sk_number[1, 1] == expected_sknum_top_right

    # Finally we can check the total charge
    fpi = 4 * np.pi
    assert sk_num == (expected_sknum_low_left + expected_sknum_top_right) / fpi


if __name__ == "__main__":

    test_sknum_spin_lattice()
