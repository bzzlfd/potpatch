import warnings
from itertools import product

import numpy as np

from potpatch.objects import AtomConfig
from potpatch.supercell import (make_supercell, closed_to_edge, 
                                infer_supercell_size)
from potpatch.utils import timing


@timing()
def check_atompos_consistency(bulk: AtomConfig, supcl: AtomConfig, 
                              tol: float = 1e-6, 
                              frozen_range: float = np.inf  # in angstrom, close to edge
                              ):
    """
    TODO: distinguish between different atom types (as an option) 
    """
    supcl_size = infer_supercell_size(bulk.lattice, supcl.lattice)
    print(f"Supercell size: {supcl_size} (guess)")

    supcl_1 = make_supercell(bulk, supcl_size)
    AL = supcl.lattice.in_unit("angstrom")
    nwarn, max_dist = 0, 0
    shifts = [np.array(shf) for shf in product(range(-1, 2), repeat=3)]
    
    pos2_transformed_shifted = (supcl.positions[:, None, :] + shifts) @ AL 
    for pos1 in supcl_1.positions:
        if closed_to_edge(supcl_1.lattice, pos1, frozen_range):
            pos1_transformed = pos1 @ AL  

            distances = np.linalg.norm(
                # Float[natoms, 27, 3] .- Float[3]
                pos2_transformed_shifted - pos1_transformed, 
                axis=2)
            min_distance = np.min(distances)
            
            if min_distance > tol:
                argmin, len_shifts = np.argmin(distances), len(shifts)
                i, j               = argmin % len_shifts, argmin // len_shifts
                nearest_pos2       = shifts[i] + supcl.positions[j]
                
                warnings.warn(
                    f"{pos1} (from bulk make_supercell) and "
                    f"{nearest_pos2} (from supcl) don't coincide. "
                    f"Discrepancy: {min_distance:.6f} angstrom")
                nwarn += 1
                max_dist = max(max_dist, min_distance)

    return nwarn, max_dist
