import unittest
import warnings

import numpy as np

from potpatch.objects import AtomConfig, Lattice, MaterialSystemInfo, VR
from potpatch.patch import inspect_ingredient
from scripts.patch_alloy import inspect_ingredient as inspect_alloy_ingredient


class TestSupercellSizeInference(unittest.TestCase):
    def make_info(self, lattice_matrix, mesh_shape):
        lattice = Lattice(np.array(lattice_matrix, dtype=float), "angstrom")
        return MaterialSystemInfo(
            atomconfig=AtomConfig(lattice=lattice),
            vr=VR(lattice=lattice, mesh=np.zeros(mesh_shape)),
        )

    def test_inspect_returns_lattice_derived_size_with_compatible_mesh(self):
        bulk = self.make_info(np.eye(3), (8, 8, 8))
        supercell = self.make_info(np.diag((4, 3, 2)), (32, 24, 16))

        size = inspect_ingredient(supercell, bulk)

        np.testing.assert_array_equal(size, np.array((4, 3, 2)))

    def test_inspect_rejects_mesh_incompatible_with_lattice_size(self):
        bulk = self.make_info(np.eye(3), (8, 8, 8))
        supercell = self.make_info(np.diag((4, 3, 2)), (16, 16, 16))

        with self.assertRaisesRegex(
                ValueError, "VR mesh magnification does not match lattice"):
            inspect_ingredient(supercell, bulk)

    def test_inspect_rejects_noninteger_lattice_size(self):
        bulk = self.make_info(np.eye(3), (8, 8, 8))
        supercell = self.make_info(np.diag((1.5, 2, 2)), (16, 16, 16))

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with self.assertRaisesRegex(ValueError, "lattices is not integer"):
                inspect_ingredient(supercell, bulk)

    def test_alloy_uses_fractional_lattice_ratio(self):
        bulk = self.make_info(np.eye(3), (8, 8, 8))
        supercell = self.make_info(np.eye(3) * 0.5, (8, 8, 8))
        supercell.charge = 0
        supercell.epsilon = np.eye(3)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            size = inspect_alloy_ingredient(supercell, bulk)

        np.testing.assert_allclose(size, (0.5, 0.5, 0.5))


if __name__ == "__main__":
    unittest.main()
