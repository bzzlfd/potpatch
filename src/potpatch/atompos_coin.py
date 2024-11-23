"atompos_coincidence"

import warnings
from itertools import product

import numpy as np

from potpatch.datatype import INTEGER
from potpatch.objects import AtomConfig
from potpatch.supercell import (make_supercell, closed_to_edge, 
                                infer_supercell_size)
from potpatch.utils import timing


@timing()
def check_atompos_consistency(bulk: AtomConfig, supcl: AtomConfig, 
                              tol: float = 1e-6, 
                              frozen_range: float = np.inf,  # in angstrom, close to edge
                              supcl_size: tuple = None,  
                              ):
    """
    TODO: distinguish between different atom types (as an option) 
    """
    _, distances, nwarn = bulk_order_mapto_supcl(bulk, supcl, warn_tol=tol)

    max_dist = 0
    for i, pos in enumerate(supcl.positions):
        if closed_to_edge(supcl.lattice, pos, frozen_range):
            distance = distances[i]
            if distance > tol:
                max_dist = max(max_dist, distance)

    return nwarn, max_dist


def bulk_order_mapto_supcl(bulk: AtomConfig, supcl: AtomConfig, 
                           warn_tol: float = 1e-6):
    """
    return a list of indices that map the closest positions of bulk 
    to the positions of supcl

    compatible with `bulk size` > `supcl size` for alloy-potental-patching 
    support

    *order refers to the supcl*: the invariance of the order of the 
    supcl.positions is important
    """
    mag, mag_f = infer_supercell_size(bulk.lattice, supcl.lattice)
    if all(mag_f <= 1):
        print("... bulk_order_mapto_supcl(): case `bulk size` > `supcl size`") 
        repeat = (1, 1, 1)
    elif all(mag_f >= 1):
        assert all(mag == mag_f), \
            "... bulk_order_mapto_supcl(): case `bulk size` < `supcl size`, "\
            f"but `mag`{mag_f} not integer"
        repeat = mag
    else:
        raise ValueError("... bulk_order_mapto_supcl(): relationship between "
                         "`bulk size` and `supcl size` is neither "
                         "`all(bulk_size > supcl_size)` OR "
                         "`all(bulk_size < supcl_size) and `mag` is integer`")

    order = np.zeros(supcl.natoms, dtype=INTEGER)
    distances = np.zeros(supcl.natoms)

    AL_bulk = bulk.lattice.in_unit("angstrom")
    AL_supcl = supcl.lattice.in_unit("angstrom")

    bulk.revise_atomsposition()
    supcl.revise_atomsposition()
    supcl.positions[supcl.positions >= 0.5] -= 1.0
    # |
    nwarn = 0
    shifts = [np.array(shf) for shf in product(
                                    range(-repeat[0], repeat[0]+1), 
                                    range(-repeat[1], repeat[1]+1), 
                                    range(-repeat[2], repeat[2]+1))]  # 8, 3
    # |
    pos_bulk_trans_shifted = (bulk.positions[:, None, :] + shifts) @ AL_bulk
    for i, pos_supcl in enumerate(supcl.positions):
        pos_supcl_trans = pos_supcl @ AL_supcl
        distances = np.linalg.norm(
            # Float[natoms, 8, 3] .- Float[3]
            pos_bulk_trans_shifted - pos_supcl_trans, axis=2)
        
        argmin, len_shifts = np.argmin(distances), len(shifts)
        i_bulk, j_bulk = argmin // len_shifts,  argmin % len_shifts
        order[i] = j_bulk
        distances[i] = np.min(distances)

        if distances[i] > warn_tol:
            nearest_pos_bulk = bulk.positions[i_bulk] + shifts[j_bulk]
            _pos_supcl, _pos_bulk = pos_supcl[:], nearest_pos_bulk[:]
            _pos_supcl[_pos_supcl < 0.0] += 1.0
            _pos_bulk[_pos_bulk < 0.0] += 1.0
            warnings.warn(
                f"Atom{i} {_pos_supcl} (from supcl) and "
                f"Atom{j_bulk} {_pos_bulk} (from bulk) don't coincide. "
                f"Discrepancy: {distances[i]:.6f} angstrom")
            nwarn += 1

    return order, distances, nwarn
