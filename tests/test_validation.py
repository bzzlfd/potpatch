"""Regression checks for explicit material-system validation."""

import unittest
from unittest.mock import patch

import numpy as np

from potpatch.objects import Atom, Lattice, MaterialSystemInfo, VR
from potpatch.patch import inspect_ingredient
from potpatch.validation import PATCH_REQUIRED_FIELDS


def lattice(scale=1.0):
    return Lattice(np.eye(3) * scale, "angstrom")


def atomconfig(scale=1.0):
    return Atom(
        natoms=1, lattice=lattice(scale), itypes=np.array([14]),
        positions=np.array([[0.0, 0.0, 0.0]]),
        moves=np.array([[1, 1, 1]]),
    )


def vr(scale=1.0):
    return VR(lattice=lattice(scale), mesh=np.zeros((2, 2, 2)))


class MaterialValidationTests(unittest.TestCase):
    def test_conflict_is_reported_without_replacing_either_lattice(self):
        atoms = atomconfig()
        potential = vr(2.0)
        info = MaterialSystemInfo(atomconfig=atoms, vr=potential)

        report = info.validate()

        self.assertIs(info.atomconfig, atoms)
        self.assertIs(info.vr, potential)
        np.testing.assert_allclose(atoms.lattice.AL, np.eye(3))
        np.testing.assert_allclose(potential.lattice.AL, np.eye(3) * 2)
        self.assertEqual(
            next(c.status for c in report.checks
                 if c.name == "atomconfig.lattice == vr.lattice"), "fail"
        )
        with self.assertRaisesRegex(ValueError, "atomconfig.lattice == vr.lattice"):
            report.raise_for_errors()

    def test_direct_child_update_is_allowed_and_rechecked(self):
        info = MaterialSystemInfo(atomconfig=atomconfig(), vr=vr())
        earlier = info.validate()
        self.assertTrue(earlier.passed)

        info.vr.lattice = lattice(3.0)

        self.assertIs(info.last_validation, earlier)
        self.assertFalse(info.validate().passed)
        np.testing.assert_allclose(info.atomconfig.lattice.AL, np.eye(3))
        np.testing.assert_allclose(info.vr.lattice.AL, np.eye(3) * 3)

    def test_in_place_array_edit_is_found_on_explicit_recheck(self):
        info = MaterialSystemInfo(atomconfig=atomconfig(), vr=vr())
        info.validate()

        info.vr.lattice.AL[0, 0] = 1.5

        self.assertTrue(info.last_validation.passed)
        self.assertFalse(info.validate().passed)

    def test_explicit_reference_does_not_rewrite_children(self):
        info = MaterialSystemInfo(atomconfig=atomconfig(), vr=vr())
        info.lattice = lattice(4.0)

        np.testing.assert_allclose(info.atomconfig.lattice.AL, np.eye(3))
        np.testing.assert_allclose(info.vr.lattice.AL, np.eye(3))
        self.assertFalse(info.validate().passed)

    def test_constructor_lattices_do_not_share_mutable_arrays(self):
        common = lattice()
        info = MaterialSystemInfo(
            lattice=common, atomconfig=Atom(lattice=common),
            vr=VR(lattice=common, mesh=np.zeros((2, 2, 2))),
        )
        info.vr.lattice.AL[0, 0] = 2.0

        np.testing.assert_allclose(common.AL, np.eye(3))
        np.testing.assert_allclose(info.atomconfig.lattice.AL, np.eye(3))
        np.testing.assert_allclose(info.lattice.AL, np.eye(3))
        self.assertFalse(info.validate().passed)

    def test_in_place_lattice_scaling_can_be_revalidated(self):
        info = MaterialSystemInfo(atomconfig=atomconfig(), vr=vr())
        info.vr.lattice *= (2, 1, 1)

        self.assertFalse(info.validate().passed)
        np.testing.assert_allclose(
            info.vr.lattice.AL, np.diag([2.0, 1.0, 1.0])
        )

    def test_missing_data_is_distinct_from_conflict(self):
        empty = MaterialSystemInfo()
        self.assertFalse(empty.validate().passed)
        self.assertEqual(empty.last_validation.checks, ())

        empty.atomconfig.lattice = lattice()
        report = empty.validate()
        self.assertEqual(report.failures, ())
        self.assertTrue(any(c.status == "insufficient" for c in report.checks))
        report.raise_for_errors()
        with self.assertRaisesRegex(ValueError, "required.vr.mesh"):
            empty.validate(required_paths=PATCH_REQUIRED_FIELDS).raise_for_errors()

    def test_shape_conflict_and_write_gate(self):
        atoms = atomconfig()
        atoms.natoms = 2
        info = MaterialSystemInfo(atomconfig=atoms)
        self.assertTrue(any(
            c.status == "fail" and c.name == "atomconfig.natoms == len(positions)"
            for c in info.validate().checks
        ))
        with patch("builtins.open") as open_file:
            with self.assertRaises(ValueError):
                atoms.write("atom.config")
            open_file.assert_not_called()

    def test_complete_info_is_ready_for_patch(self):
        info = MaterialSystemInfo(atomconfig=atomconfig(), vr=vr())
        report = info.validate(required_paths=PATCH_REQUIRED_FIELDS)
        report.raise_for_errors()
        self.assertTrue(report.passed)

    def test_vr_write_rejects_missing_mesh_before_opening(self):
        potential = VR(lattice=lattice())
        with patch("builtins.open") as open_file:
            with self.assertRaisesRegex(ValueError, "required.vr.mesh"):
                potential.write("OUT.VR")
            open_file.assert_not_called()

    def test_patch_inspection_rejects_conflicting_sources(self):
        bulk = MaterialSystemInfo(atomconfig=atomconfig(), vr=vr())
        supcl = MaterialSystemInfo(atomconfig=atomconfig(2.0), vr=vr(3.0))

        with self.assertRaisesRegex(ValueError, "atomconfig.lattice == vr.lattice"):
            inspect_ingredient(supcl, bulk)

    def test_frozen_atom_inspection_requires_positions(self):
        bulk = MaterialSystemInfo(atomconfig=atomconfig(), vr=vr())
        supcl = MaterialSystemInfo(atomconfig=atomconfig(), vr=vr())
        supcl.atomconfig.positions = None

        with self.assertRaisesRegex(ValueError, "required.atomconfig.positions"):
            inspect_ingredient(supcl, bulk, frozen_confirm=1.0)


if __name__ == "__main__":
    unittest.main()
