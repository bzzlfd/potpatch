import warnings
from itertools import product

import numpy as np

from potpatch.objects import AtomConfig
from potpatch.supercell import (modify_supercell, make_supercell, 
                                closed_to_edge, infer_supercell_size)
from potpatch.utils import timing


@timing()
def check_atompos_consistency(bulk: AtomConfig, supcl: AtomConfig, 
                              tol: float = 1e-6, 
                              frozen_range: float = np.inf  # in angstrom, close to edge
                              ):
    supcl_size = infer_supercell_size(bulk.lattice, supcl.lattice)
    print(f"Supercell size: {supcl_size} (guess)")

    supcl_1 = make_supercell(bulk, supcl_size)
    supcl_2 = supcl
    AL = supcl.lattice.in_unit("angstrom")
    nwarn, max_dist = 0, 0
    shifts = [np.array(shf) for shf in product(range(-1, 2), repeat=3)]
    
    for pos1 in supcl_1.positions:
        if closed_to_edge(supcl_1.lattice, pos1, frozen_range):
            # 使用矢量化计算所有 `supcl_2` 中原子位置的平移距离
            pos1_transformed = pos1 @ AL  # 先将 `pos1` 转换到笛卡尔坐标系
            pos2_transformed = supcl.positions @ AL  # 转换整个 `supcl` 的原子坐标

            # 计算 `pos1` 到 `supcl_2` 所有原子及平移的距离
            distances = np.linalg.norm(pos2_transformed[:, None, :] + shifts - pos1_transformed, axis=2)
            min_distance = np.min(distances)  # 找到最小距离
            
            # 如果最小距离大于容差，则记录警告信息
            if min_distance > tol:
                nearest_shift = shifts[np.argmin(distances) // len(supcl.positions)]
                nearest_pos2 = supcl.positions[np.argmin(distances) % len(supcl.positions)] + nearest_shift
                nearest_pos2_transformed = nearest_pos2 @ AL
                
                warnings.warn(
                    f"{pos1} (from bulk make_supercell) and {nearest_pos2_transformed} (from supcl) don't coincide. "
                    f"Discrepancy: {min_distance:.6f} angstrom")
                nwarn += 1
                max_dist = max(max_dist, min_distance)

    return nwarn, max_dist
