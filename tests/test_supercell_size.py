import unittest
import warnings

import numpy as np

from potpatch.objects import Atom, Lattice, MaterialSystemInfo, VR
from potpatch.patch import (inspect_ingredient, patch_vr,
                            _resample_periodic_mesh)
from scripts.patch_alloy import inspect_ingredient as inspect_alloy_ingredient


class TestSupercellSizeInference(unittest.TestCase):
    def make_info(self, lattice_matrix, mesh_shape):
        lattice = Lattice(np.array(lattice_matrix, dtype=float), "angstrom")
        return MaterialSystemInfo(
            atomconfig=Atom(lattice=lattice),
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

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            with self.assertRaisesRegex(ValueError, "exceeds 5%"):
                inspect_ingredient(supercell, bulk)
        self.assertTrue(any("exceeds 5%" in str(w.message)
                            for w in caught))
        self.assertEqual(supercell.vr.mesh.shape, (16, 16, 16))

    def test_inspect_resamples_multi_point_difference_before_patching(self):
        bulk = self.make_info(np.eye(3), (32, 32, 32))
        supercell = self.make_info(np.diag((2, 1, 1)), (66, 31, 32))

        def potential(shape, x_cycles):
            x, y, z = np.meshgrid(
                *(np.arange(n) / n for n in shape), indexing="ij")
            return (1.5 + np.cos(2 * np.pi * x_cycles * x)
                    + 0.2 * np.sin(2 * np.pi * y)
                    + 0.3 * np.cos(2 * np.pi * z))

        bulk.vr.mesh = potential((32, 32, 32), 1)
        source_mesh = potential((66, 31, 32), 2)
        supercell.vr.mesh = source_mesh
        original = source_mesh.copy()

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            size = inspect_ingredient(supercell, bulk)

        np.testing.assert_array_equal(size, (2, 1, 1))
        np.testing.assert_allclose(supercell.vr.mesh,
                                   np.tile(bulk.vr.mesh, (2, 1, 1)), atol=1e-12)
        np.testing.assert_array_equal(source_mesh, original)
        self.assertFalse(caught)
        patched = patch_vr(supercell.vr, bulk.vr, size, (4, 1, 1))
        self.assertEqual(patched.mesh.shape, (128, 32, 32))

    def test_inspect_rejects_spacing_difference_over_five_percent(self):
        bulk = self.make_info(np.eye(3), (32, 32, 32))
        supercell = self.make_info(np.diag((2, 1, 1)), (68, 32, 32))
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            with self.assertRaisesRegex(ValueError, "exceeds 5%"):
                inspect_ingredient(supercell, bulk)
        self.assertEqual(supercell.vr.mesh.shape, (68, 32, 32))
        self.assertTrue(any("spacing difference" in str(w.message)
                            for w in caught))

    def test_patch_vr_rejects_unreconciled_mesh(self):
        bulk = self.make_info(np.eye(3), (8, 8, 8))
        supercell = self.make_info(np.diag((2, 1, 1)), (17, 8, 8))
        with self.assertRaisesRegex(ValueError, "call inspect_ingredient first"):
            patch_vr(supercell.vr, bulk.vr, (2, 1, 1), (4, 1, 1))

    def test_fractional_alloy_patch_keeps_existing_grid_behavior(self):
        bulk = self.make_info(np.eye(3), (8, 8, 8))
        alloy = self.make_info(np.eye(3) * 0.5, (8, 8, 8))
        alloy.vr.mesh.fill(2.0)
        result = patch_vr(alloy.vr, bulk.vr, (0.5, 0.5, 0.5), (1, 1, 1))
        np.testing.assert_array_equal(result.mesh, alloy.vr.mesh)

    def test_fourier_resampling_preserves_nyquist_mode_and_mean(self):
        x = np.arange(8) / 8
        source = np.cos(2 * np.pi * 4 * x)[:, None, None]
        source = np.broadcast_to(source, (8, 3, 3)) + 2.0
        output = _resample_periodic_mesh(source, (9, 3, 3))
        expected = np.cos(2 * np.pi * 4 * np.arange(9) / 9)
        np.testing.assert_allclose(output[:, 0, 0], expected + 2.0,
                                   atol=1e-12)
        self.assertAlmostEqual(output.mean(), source.mean())

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
