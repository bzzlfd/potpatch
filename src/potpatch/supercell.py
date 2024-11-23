from itertools import product
from warnings import warn

import numpy as np
import numpy.typing as npt

from potpatch.objects import AtomConfig, Lattice
from potpatch.datatype import REAL_8, INTEGER


def which_lattice_is_bulk(lat1: Lattice, lat2: Lattice):
    AL_1 = lat1.in_unit("angstrom")
    AL_2 = lat2.in_unit("angstrom")
    
    magM = AL_2 @ np.linalg.inv(AL_1)        
    if all(magM[i, i] >= 1 for i in range(3)): 
        whichLatticeIsBulk = 1 
    elif all(magM[i, i] <= 1 for i in range(3)):
        whichLatticeIsBulk = 2
    else:
        assert False, "can't infer which lattice is bulk"

    return whichLatticeIsBulk


def infer_supercell_size(bulk: Lattice, supcl: Lattice):
    """
    infer the size of supercell from bulk and supercell *lattices* only
    returns:
        mag: integer array, size of supercell
        mag_f: float array, size of supercell
    """
    AL_b = bulk.in_unit("angstrom")
    AL_s = supcl.in_unit("angstrom")
    magM = AL_s @ np.linalg.inv(AL_b)
    
    mag = np.array([round(magM[i, i], 0) for i in range(3)], dtype=INTEGER)
    mag_f = np.array([magM[i, i] for i in range(3)], dtype=REAL_8)
    
    for i, j in product(range(3), repeat=2):
        if i != j:
            assert abs(magM[i, j]) < 1e-6, \
                f"AL_s / AL_b = {magM}"
        if i == j and abs(magM[i, i] - round(magM[i, i])) > 1e-6:
            warn(f"warning magM[i, i]{mag_f} is not integer", stacklevel=3)
    return mag, mag_f


def make_supercell(bulkAtom: AtomConfig, m123) -> AtomConfig:
    "``m123`` should be a tuple contains 3 integer numbers, e.g. ``(4,4,4)``"
    supclAtom = bulkAtom * m123
    supclAtom.revise_atomsposition()
    return supclAtom


def modify_supercell(supclAtom: AtomConfig,
                     frozen_edge_width=1.0, 
                     mv2outsider=None) -> AtomConfig:
    """
    for supercell RELAX job: freeze atoms whose positions are closed to
    the pathcing edge(interface between supercell and bulk)
    frozen_edge_width: in angstrom
    """
    rtac = supclAtom * (1, 1, 1)  # make a copy

    for (i, pos) in enumerate(rtac.positions):
        rtac.moves[i] = np.array([1, 1, 1])
    
    count = 0
    for (i, pos) in enumerate(rtac.positions):
        if closed_to_edge(rtac.lattice, pos, frozen_edge_width):
            rtac.moves[i] = np.array([0, 0, 0])
            if mv2outsider is not None:
                assert int(mv2outsider) == mv2outsider
                rtac.itypes[i] = INTEGER(mv2outsider)
            count += 1
    rtac._natoms_frozen = count

    rtac.revise_atomsposition()
    rtac.sort_atomposition()
    print("...", f"modify_supercell(): {count} atoms frozen")
    return rtac


def closed_to_edge(lattice: Lattice, position, frozen_edge_width):
    position = np.array(position)
    position %= 1
    indx = position >= 0.5
    position[indx] -= 1.0

    AL = lattice.in_unit("angstrom")
    iAL = (np.linalg.inv(AL))
    iuAL = (iAL / np.linalg.norm(iAL, axis=0)).T
    L = np.array([np.abs(np.dot(AL[i], iuAL[i])) for i in range(3)])
    # print(L*position)
    if any(position*L <= -0.5*L+frozen_edge_width) \
            or any(position*L >= 0.5*L-frozen_edge_width):
        return True
    else:
        return False
