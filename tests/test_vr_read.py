import tempfile
import unittest
import warnings

import numpy as np

from potpatch.objects import Lattice, VR
from potpatch.utils import write_fortran_binary_block
from potpatch.datatype import INTEGER_OUT


def vr() -> VR:
    lattice = Lattice(np.diag((5.2, 5.2, 5.2)), "angstrom")
    mesh = np.arange(2*2*2, dtype=float).reshape(2, 2, 2)
    return VR(lattice=lattice, mesh=mesh)


def write_vr_with_nstate(filename, potential: VR, nstate: int | None):
    """
    write an OUT.VR(PWmat fmt, nnodes=1) whose first record has 4 integers,
    or 5 integers(`nstate` appended) when `nstate` is not None.
    if nstate > 1, `nstate` copies of the mesh data are written.
    """
    with open(filename, "bw") as io:
        header = [*potential.mesh.shape, 1] \
            + ([] if nstate is None else [nstate])
        write_fortran_binary_block(io, np.array(header, dtype=INTEGER_OUT))

        AL = potential.lattice.in_unit("angstrom")
        write_fortran_binary_block(io, np.reshape(AL, AL.size))

        mesh1d = np.reshape(potential.mesh, potential.mesh.size)
        for _ in range(1 if nstate is None else nstate):
            write_fortran_binary_block(io, mesh1d)


class TestVRReadHeader(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def path(self, name="OUT.VR"):
        import os
        return os.path.join(self.tmp.name, name)

    def assert_readback_matches(self, read_vr: VR):
        expected = vr()
        np.testing.assert_allclose(read_vr.mesh, expected.mesh)
        np.testing.assert_allclose(
            read_vr.lattice.AL_AU, expected.lattice.AL_AU)
        self.assertEqual(read_vr.nnodes, 1)

    def test_read_4_int_header(self):
        path = self.path()
        vr().write(path, nnodes=1)

        read_vr = VR(filename=path)

        self.assert_readback_matches(read_vr)
        self.assertEqual(read_vr.nstate, 1)

    def test_read_5_int_header(self):
        path = self.path()
        write_vr_with_nstate(path, vr(), nstate=1)

        read_vr = VR(filename=path)

        self.assert_readback_matches(read_vr)
        self.assertEqual(read_vr.nstate, 1)

    def test_read_5_int_header_multi_state_warns_and_keeps_first(self):
        path = self.path()
        write_vr_with_nstate(path, vr(), nstate=2)

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            read_vr = VR(filename=path)

        self.assertTrue(any("reads only the first" in str(w.message)
                            for w in caught))
        self.assert_readback_matches(read_vr)
        self.assertEqual(read_vr.nstate, 2)

    def test_read_unexpected_int_count_raises(self):
        path = self.path()
        with open(path, "bw") as io:
            write_fortran_binary_block(io, np.array([2, 2, 2], dtype=INTEGER_OUT))

        with self.assertRaisesRegex(ValueError, "3 integers"):
            VR(filename=path)


if __name__ == "__main__":
    unittest.main()
