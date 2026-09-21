"""The short Atom name denotes a complete structure and preserves old imports."""

import unittest

import numpy as np

from potpatch import Atom, AtomConfig, Lattice, MaterialSystemInfo
from potpatch.objects import Atom as ObjectAtom


class AtomNameTests(unittest.TestCase):
    def test_public_name_and_compatibility_alias(self):
        self.assertIs(Atom, ObjectAtom)
        self.assertIs(AtomConfig, Atom)
        self.assertEqual(Atom.__name__, "Atom")
        self.assertIsInstance(MaterialSystemInfo().atomconfig, Atom)

    def test_atom_represents_entire_crystal_structure(self):
        atoms = Atom(
            natoms=2,
            lattice=Lattice(np.eye(3), "angstrom"),
            itypes=np.array([14, 8]),
            positions=np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
            moves=np.ones((2, 3), dtype=int),
        )

        supercell = atoms * (2, 1, 1)

        self.assertIsInstance(supercell, Atom)
        self.assertEqual(supercell.natoms, 4)
        np.testing.assert_array_equal(supercell.itypes, [14, 8, 14, 8])
        np.testing.assert_allclose(
            supercell.lattice.in_unit("angstrom"), np.diag([2.0, 1.0, 1.0])
        )


if __name__ == "__main__":
    unittest.main()
