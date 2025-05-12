from fractions import Fraction
from typing import Tuple
from os.path import join, abspath
from pathlib import Path

from numpy import prod, array, zeros

from potpatch import (AtomConfig, modify_supercell)


def in_imbox(pos, planl, planr):
    if all(planl <= pos) and all(pos < planr):
        return True
    else:
        return False


def make_supercell(bulkAtom: AtomConfig, 
                   m123: Tuple[Fraction, Fraction, Fraction],
                   vr_mesh_n123 = None, 
                   ) -> AtomConfig:
    m123 = array(m123, dtype=Fraction)
    fm123 = m123.astype(float)
    natoms = bulkAtom.natoms * prod(m123)
    assert natoms.denominator == 1, f"natoms({natoms}) should be an integer"
    natoms = int(natoms)

    lattice = bulkAtom.lattice * fm123
    atoms_itype = zeros(natoms, dtype=int)
    atoms_position = zeros((natoms, 3), dtype=float)
    atoms_move = zeros((natoms, 3), dtype=int)

    cnt = 0
    planl, planr = -fm123/2, +fm123/2
    bulkAtom.revise_atomsposition()
    for ia in range(bulkAtom.natoms):
        pos = bulkAtom.positions[ia]
        pos[pos >= 0.5] -= 1.0
        if in_imbox(pos, planl, planr):
            atoms_itype[cnt] = bulkAtom.itypes[ia]
            atoms_position[cnt] = pos
            atoms_move[cnt] = bulkAtom.moves[ia]
            cnt += 1
    assert cnt == natoms, f"cnt({cnt}) should equal natoms({natoms})"

    atoms_position /= fm123
    supclAtom = AtomConfig(natoms=natoms, lattice=lattice, 
                           itypes=atoms_itype, positions=atoms_position, 
                           moves=atoms_move)
    supclAtom.revise_atomsposition()

    return supclAtom


if __name__ == '__main__':
    rootdir = Path(__file__).parent
    bulk = ""
    supcl = ""
    debug = 56  # set frozen atoms to 56(Ba)

    bulkAC = AtomConfig(filename=join(rootdir, bulk))
    supclAC = make_supercell(
                        bulkAtom=bulkAC, 
                        m123=(Fraction("1/4"), Fraction(1, 4), Fraction(0.25)))
    modify_supercell(supclAC, 
                     1.0, 
                     debug)  # for debug
    supclAC.write(filename=join(rootdir, supcl), 
                        comment=f"from {abspath(rootdir)}/{bulk}", 
                        fmt="PWmat")
