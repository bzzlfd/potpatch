import os

import numpy as np

from potpatch import VATOM, AtomConfig, Lattice
from potpatch.atompos_coin import bulk_order_mapto_supcl
from potpatch.supercell import infer_supercell_size
from potpatch.diff_vatom import write_diffvatom


def read_vatom(bulk_dir, supcl_dir):
    # read atom.config  for `lattice` and `order`
    bulk_ac = os.path.join(bulk_dir, "atom.config")
    supcl_ac = os.path.join(supcl_dir, "atom.config")
    bulk_ac = AtomConfig().read(bulk_ac)
    supcl_ac = AtomConfig().read(supcl_ac)
    bulk_lattice = bulk_ac.lattice
    supcl_lattice = supcl_ac.lattice

    # read OUT.VATOM
    bulk_vatom = os.path.join(bulk_dir, "OUT.VATOM")
    supcl_vatom = os.path.join(supcl_dir, "OUT.VATOM")
    bulk_vatom = VATOM().read_vatom(bulk_vatom)
    supcl_vatom = VATOM().read_vatom(supcl_vatom)

    # 
    order, _, _ = bulk_order_mapto_supcl(bulk_ac, supcl_ac, warn_tol=np.inf)
    r          = np.zeros(supcl_vatom.natoms, dtype=np.float64)
    ξ          = np.zeros(supcl_vatom.natoms, dtype=np.float64)
    dv_bulk    = np.zeros(supcl_vatom.natoms, dtype=np.float64)
    dv_supcl   = np.zeros(supcl_vatom.natoms, dtype=np.float64)

    for i in range(supcl_vatom.natoms):
        dv_bulk[i] = bulk_vatom.vatoms_au[order[i]]
        dv_supcl[i] = supcl_vatom.vatoms_au[i]

        pos = supcl_vatom.positions[i]
        pos = pos % 1.0
        pos[pos >= 0.5] -= 1.0
        r[i] = np.linalg.norm(pos @ supcl_lattice.in_unit("atomic unit"))
        ξ[i] = np.nan

    return r, ξ, dv_bulk, dv_supcl, bulk_vatom, supcl_vatom, order


if __name__ == "__main__":
    bulk_dir = "bulk"
    supcl_dir = "supcl"
    supcl_corr_dir = "supcl_corr"
    charge = 0

    epsilon = np.full((3, 3), np.nan, dtype=np.float64)

    r, ξ, dv_bulk, dv_supcl, bulk_vatom, supcl_vatom, order = \
        read_vatom(bulk_dir, supcl_dir)
    _, _, _,       dv_supcl_corr, _,     _,           _ = \
        read_vatom(bulk_dir, supcl_corr_dir)

    write_diffvatom("diff_vatom.dat", 
                    bulk_vatom, supcl_vatom, 
                    order, epsilon, charge,
                    r, ξ, 
                    dv_bulk, dv_supcl, dv_supcl_corr, 
                    )
