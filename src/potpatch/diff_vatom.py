"""
a script for comparing the difference between the supercell.vatom and 
the bulk.vatom

*closest atom* match but not 1-1 map

using Gaussian weighted averages, which is a different implementation 
from the PWmat `OUT.VATOM`
"""
from os.path import exists
from itertools import product

import numpy as np
from numba import jit, guvectorize, prange

from potpatch.objects import Lattice, VR, AtomConfig, MaterialSystemInfo
from potpatch.supercell import infer_supercell_size
from potpatch.datatype import INTEGER
from potpatch.constant import HA, BOHR
from potpatch.utils import timing


@timing()
def diff_vatom(bulk: MaterialSystemInfo, supcl: MaterialSystemInfo, 
               sigma: float = 1.):
    size = infer_supercell_size(bulk.lattice, supcl.lattice)
    ac_bulk  = bulk.atomconfig * size
    ac_supcl = supcl.atomconfig
    vr_bulk  = bulk.vr * size
    vr_supcl = supcl.vr
    lattice = supcl.lattice
    n123 = supcl.vr.n123

    r = np.zeros(ac_supcl.natoms, dtype=np.float64)
    AL = lattice.in_unit("atomic unit")
    ac_supcl.revise_atomsposition()
    ac_supcl.positions[ac_supcl.positions >= 0.5] -= 1.0
    for i in range(ac_supcl.natoms):
        shifts = [np.array(shf) for shf in product(range(-1, 1), repeat=3)]
        coords = (ac_supcl.positions[i, :] + shifts) @ AL
        r[i] = np.min(np.linalg.norm(coords, axis=1))

    ac_bulk = closest_resort(AL, ac_supcl, ac_bulk)

    iAL = (np.linalg.inv(AL))
    iuAL = (iAL / np.linalg.norm(iAL, axis=0)).T
    L = np.array([np.abs(np.dot(AL[i], iuAL[i])) for i in range(3)])
    window = INTEGER(6 * sigma / L * n123) + 1
    if any(window > n123):
        print(f"Warning: window size({window=}) exceeds supercell size{n123=}")
    window[window > n123] = n123[window > n123]  # TODO edge case
    w1, w2, w3 = window // 2
    window_lb = - np.array([w1, w2, w3])
    window_ur = np.array([w1, w2, w3])

    diff_vatom = np.zeros((ac_supcl.natoms, 3), dtype=np.float64)
    ac_supcl.revise_atomsposition()
    ac_supcl.positions[ac_supcl.positions >= 0.5] -= 1.0
    ac_bulk.revise_atomsposition()
    ac_bulk.positions[ac_bulk.positions >= 0.5] -= 1.0
    for i in prange(ac_supcl.natoms):
        pos = ac_bulk.positions[i]
        center = INTEGER(np.ceil(pos * n123))
        diff_vatom[i, 0] = \
            gaussian_integrate(AL, vr_bulk.mesh, 
                               center+window_lb, center+window_ur, pos, sigma)
        pos = ac_supcl.positions[i]
        center = INTEGER(np.ceil(pos * n123))
        diff_vatom[i, 1] = \
            gaussian_integrate(AL, vr_supcl.mesh, 
                               center+window_lb, center+window_ur, pos, sigma)

    diff_vatom[:, 2] = diff_vatom[:, 1] - diff_vatom[:, 0]

    return r, diff_vatom, ac_bulk, ac_supcl


def closest_resort(AL, supcl: AtomConfig, bulk: AtomConfig):
    """may not resort
    operate `bulk`
    """
    supcl.revise_atomsposition()
    supcl.positions[supcl.positions >= 0.5] -= 1.0
    bulk.revise_atomsposition()

    order = np.zeros(len(supcl.positions), dtype=np.int32)
    shifts = [np.array(shf) for shf in product(range(-1, 1), repeat=3)]
    
    pos2_transformed_shifted = (bulk.positions[:, None, :] + shifts) @ AL 
    for i, pos1 in enumerate(supcl.positions):
        pos1_transformed = pos1 @ AL  

        distances = np.linalg.norm(
            # Float[8, 3] .- Float[3]
            pos2_transformed_shifted - pos1_transformed, 
            axis=2)
        
        argmin, len_shifts = np.argmin(distances), len(shifts)
        order[i]           = argmin // len_shifts 

    bulk.revise_atomsposition()
    bulk.positions[bulk.positions >= 0.5] -= 1.0
    bulk.positions = bulk.positions[order]
    bulk.itypes = bulk.itypes[order]

    return bulk


def gaussian_integrate(AL, mesh: np.ndarray, window_lb, window_ur, pos, sigma):
    Ω = np.abs(np.linalg.det(AL))
    n123 = np.array(mesh.shape)
    dΩ = Ω / np.prod(n123)
    
    mesh_indices = np.mgrid[window_lb[0]:window_ur[0], 
                            window_lb[1]:window_ur[1], 
                            window_lb[2]:window_ur[2]]  #

    flat_indices = mesh_indices.reshape(3, -1).T
    
    mesh_pos = (flat_indices / n123 - pos) @ AL.T  
    r2 = np.sum(mesh_pos ** 2, axis=1)  
    
    gaussian_vals = np.exp(-r2 / (2 * sigma ** 2))
    
    vr_sum = np.sum(mesh[flat_indices[:, 0], 
                         flat_indices[:, 1], 
                         flat_indices[:, 2]] 
                    * gaussian_vals)  #
    gaussian_sum = np.sum(gaussian_vals)
    
    vr_avg = vr_sum / gaussian_sum
    
    # gaussian_sum *= (np.sqrt(2 * np.pi) * sigma)**-3 * dΩ
    # if (debug:= True) or gaussian_sum < 0.9:
    #     print(f"{i=}, {j=}, {k=}, guassian_sum={gaussian_sum}")
    
    return vr_avg


def write_diffvatom(filename: str, supcl: AtomConfig, bulk: AtomConfig,
                    r: np.ndarray, diff_vatom: np.ndarray):
    """
    natoms
    ...
    # bulk.ac.itypes, bulk.ac.positions, bulk.vatom, supcl.ac.itypes, supcl.ac.positions, supcl.vatom, diff_vatom,
    # metadata
    """
    if exists(filename):
        pass  # TODO overwrite or rename

    r           = r * BOHR
    diff_vatom  = diff_vatom * HA

    supcl.revise_atomsposition()
    bulk.revise_atomsposition()
    with open(filename, "w") as f:
        f.write(f" {supcl.natoms} atoms\n")
        for i in range(supcl.natoms):
            f.write(f"{r[i]:15.8f}")
            f.write(f"{bulk.itypes[i]:20d} {bulk.positions[i, 0]:15.8f} {bulk.positions[i, 1]:15.8f} {bulk.positions[i, 2]:15.8f} {diff_vatom[i, 0]:15.8f}")
            f.write(f"{supcl.itypes[i]:20d} {supcl.positions[i, 0]:15.8f} {supcl.positions[i, 1]:15.8f} {supcl.positions[i, 2]:15.8f} {diff_vatom[i, 1]:15.8f}")
            f.write(f"{diff_vatom[i, 2]:20.8f}\n")
        f.write("# r, bulk.ac.itypes, bulk.ac.positions, bulk.vatom, supcl.ac.itypes, supcl.ac.positions, supcl.vatom, diff_vatom\n")
