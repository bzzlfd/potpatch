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
from numpy import sum, einsum, prod, abs, sqrt, exp, pi
from numpy.linalg import norm, inv, det, eigvalsh
from numba import jit, guvectorize, prange

from potpatch.objects import Lattice, VR, AtomConfig, MaterialSystemInfo
from potpatch.atompos_coin import bulk_order_mapto_supcl
from potpatch.datatype import INTEGER
from potpatch.constant import HA, BOHR
from potpatch.utils import timing


@timing()
def diff_vatom(bulk: MaterialSystemInfo, supcl: MaterialSystemInfo, 
               sigma: float = 1.):
    ac_bulk  = bulk.atomconfig
    ac_supcl = supcl.atomconfig
    vr_bulk  = bulk.vr
    vr_supcl = supcl.vr
    epsilon = supcl.epsilon
    n123 = supcl.vr.n123

    # gaussian integration window
    AL = supcl.lattice.in_unit("atomic unit")
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

    order, _, _ = bulk_order_mapto_supcl(ac_supcl, ac_bulk, warn_tol=np.inf)
    r          = np.zeros(ac_supcl.natoms, dtype=np.float64)
    ξ          = np.zeros(ac_supcl.natoms, dtype=np.float64)
    diff_vatom = np.zeros((ac_supcl.natoms, 3), dtype=np.float64)
    for i in prange(ac_supcl.natoms):
        # bulk
        pos = ac_bulk.positions[order[i]]
        center = INTEGER(np.ceil(pos * n123))
        _, _, diff_vatom[i, 0] = \
            gaussian_integrate(bulk.lattice.in_unit("atomic unit"),
                               vr_bulk.mesh, 
                               center+window_lb, center+window_ur, 
                               pos, epsilon, sigma)
        # supercell, r, ξ
        pos = ac_supcl.positions[i]
        center = INTEGER(np.ceil(pos * n123))
        r[i], ξ[i], diff_vatom[i, 1] = \
            gaussian_integrate(supcl.lattice.in_unit("atomic unit"),
                               vr_supcl.mesh, 
                               center+window_lb, center+window_ur, 
                               pos, epsilon, sigma)

    diff_vatom[:, 2] = diff_vatom[:, 1] - diff_vatom[:, 0]

    return r, ξ, diff_vatom, ac_bulk, ac_supcl


def gaussian_integrate(AL, mesh: np.ndarray, window_lb, window_ur, 
                       pos, epsilon, sigma):
    """
    pos in [-0.5, 0.5), corresponding to the image charge correction

    n(ξ) = 1/√(ε₁ε₂ε₃) √(p/π)³ exp(-p ξ²) 
        ξ    = √(r₁²/ε₁ + r₂²/ε₂ + r₃²/ε₃)
        p    = εₘₐₓ/(2σ²)
    """
    Ω = abs(det(AL))
    n123 = np.array(mesh.shape)
    dΩ = Ω / prod(n123)
    
    inv_eps = inv(epsilon)
    eps_ii = eigvalsh(epsilon)
    
    mesh_indices = np.mgrid[window_lb[0]:window_ur[0], 
                            window_lb[1]:window_ur[1], 
                            window_lb[2]:window_ur[2]]  #

    flat_indices = mesh_indices.reshape(3, -1).T  # [m, 3]
    flat_indices = flat_indices % n123 

    r = pos @ AL.T
    ξ, r = sqrt(r.T @ inv_eps @ r), norm(r)
    
    rT = (flat_indices / n123 - pos) @ AL.T  # [m, 3]
    ξ2 = einsum("ij,ij->i", rT @ inv_eps, rT)
    p = max(eps_ii) / (2 * sigma ** 2)
    
    gaussian_vals = np.exp(- p * ξ2)
    
    vr_sum = sum(mesh[flat_indices[:, 0], 
                      flat_indices[:, 1], 
                      flat_indices[:, 2]] * gaussian_vals)
    gaussian_sum = sum(gaussian_vals)
    
    vr_avg = vr_sum / gaussian_sum
    
    gaussian_sum *= 1/sqrt(prod(eps_ii)) * sqrt(p / pi)**3 * dΩ
    if (debug := True) and gaussian_sum < 0.9:
        print(f"{window_lb=}, {window_ur=}, guassian_sum={gaussian_sum}")
    
    return r, ξ, vr_avg


def write_diffvatom(filename: str, supcl: AtomConfig, bulk: AtomConfig, 
                    epsilon: np.ndarray, 
                    r: np.ndarray, ξ: np.ndarray, diff_vatom: np.ndarray):
    """
    natoms
    ...
    # bulk.ac.itypes, bulk.ac.positions, bulk.vatom, supcl.ac.itypes, supcl.ac.positions, supcl.vatom, diff_vatom,
    # metadata
    """
    if exists(filename):
        pass  # TODO overwrite or rename

    ξ[ξ < 1e-10] = 1e-10  # avoid divide by zero
    pointv      = 1 / ξ / sqrt(prod(eigvalsh(epsilon)))

    r           = r * BOHR
    ξ           = ξ * BOHR
    pointv      = pointv * HA        
    diff_vatom  = diff_vatom * HA

    supcl.revise_atomsposition()
    bulk.revise_atomsposition()
    with open(filename, "w") as f:
        f.write(f" {supcl.natoms} atoms\n")
        f.write(" EPSILON\n")
        for i in range(3):
            f.write(f"{epsilon[i, 0]:15.8f} {epsilon[i, 1]:15.8f} {epsilon[i, 2]:15.8f}\n")
        f.write(" diff_VATOM  # r, ξ, 1/√(ε₁ε₂ε₃)/ξ, "
                "bulk.ac.itypes, bulk.ac.positions[0:3], bulk.vatom, "
                "supcl.ac.itypes, supcl.ac.positions[0:3], supcl.vatom, "
                "diff_vatom\n")
        for i in range(supcl.natoms):
            f.write(f"{r[i]:15.8f} {ξ[i]:15.8f} {pointv[i]:15.8f} ")
            f.write(f"{bulk.itypes[i]:20d} {bulk.positions[i, 0]:15.8f} {bulk.positions[i, 1]:15.8f} {bulk.positions[i, 2]:15.8f} {diff_vatom[i, 0]:15.8f}")
            f.write(f"{supcl.itypes[i]:20d} {supcl.positions[i, 0]:15.8f} {supcl.positions[i, 1]:15.8f} {supcl.positions[i, 2]:15.8f} {diff_vatom[i, 1]:15.8f}")
            f.write(f"{diff_vatom[i, 2]:20.8f}\n")
