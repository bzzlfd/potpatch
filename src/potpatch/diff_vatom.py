"""
a script for comparing the difference between the supercell.vatom and 
the bulk.vatom

*closest atom* match but not 1-1 map

using Gaussian weighted averages, which is a different implementation 
from the PWmat `OUT.VATOM`

in guassian_integrate and get_window, we use the same ε scaling
"""
from os.path import exists
from itertools import product

import numpy as np
from numpy import sum, einsum, prod, abs, sqrt, exp, pi
from numpy.linalg import norm, inv, det, eigvalsh
from numpy.typing import NDArray
from numba import jit, guvectorize, prange

from potpatch.objects import Lattice, VR, AtomConfig, VATOM, MaterialSystemInfo
from potpatch.atompos_coin import bulk_order_mapto_supcl
from potpatch.supercell import infer_supercell_size
from potpatch.datatype import INTEGER, INTEGER_8
from potpatch.constant import HA, BOHR
from potpatch.utils import timing


@timing()
def diff_vatom(bulk: MaterialSystemInfo, supcl: MaterialSystemInfo, 
               sigma: float = 1.):
    "calulate and return in atomic unit: r, ξ, diff_vatom "
    ac_bulk  = bulk.atomconfig
    ac_supcl = supcl.atomconfig
    vr_bulk  = bulk.vr
    vr_supcl = supcl.vr
    epsilon = supcl.epsilon
    n123 = supcl.vr.n123

    # mag
    _, mag_f = infer_supercell_size(bulk.lattice, supcl.lattice)

    # gaussian integration window
    window = get_window(supcl.lattice.in_unit("atomic unit"), 
                        n123, epsilon, sigma)
    if debug := False:
        window = np.ones(3, dtype=np.int64) * 2

    order, _, _ = bulk_order_mapto_supcl(ac_bulk, ac_supcl, warn_tol=np.inf)
    r          = np.zeros(ac_supcl.natoms, dtype=np.float64)
    ξ          = np.zeros(ac_supcl.natoms, dtype=np.float64)
    dv_bulk    = np.zeros(ac_supcl.natoms, dtype=np.float64)
    dv_supcl   = np.zeros(ac_supcl.natoms, dtype=np.float64)
    for i in prange(ac_supcl.natoms):
        # bulk
        pos = ac_bulk.positions[order[i]]
        _, _, dv_bulk[i] = \
            gaussian_integrate(bulk.lattice.in_unit("atomic unit"),
                               vr_bulk.mesh, window, pos, epsilon, sigma)
        # supercell, r, ξ
        pos = ac_supcl.positions[i]
        r[i], ξ[i], dv_supcl[i] = \
            gaussian_integrate(supcl.lattice.in_unit("atomic unit"),
                               vr_supcl.mesh, window, pos, epsilon, sigma)

    return r, ξ, dv_bulk, dv_supcl, ac_bulk, ac_supcl, order


def get_window(AL, n123, epsilon, sigma):
    """
    epsilon supports cos_theta, scaling is not needed here
    """
    iAL = (np.linalg.inv(AL))
    iuAL = (iAL / np.linalg.norm(iAL, axis=0)).T

    L, cos_theta = np.zeros(3), np.zeros(3)
    for i in range(3):
        a, p = AL[i], iuAL[i]
        cos_theta[i] = p.T @ epsilon @ p / norm(p) / norm(epsilon @ p)
        L[i] = np.dot(p, a)
    if debug := False:
        print(f"{L=}, {cos_theta=}")
    L = L / cos_theta

    window = INTEGER(6 * sigma / L * n123) + 1

    return window


# from potpatch.utils import gen_counter
# cnt = gen_counter()
def gaussian_integrate(AL, mesh: np.ndarray, window,
                       pos, epsilon, sigma):
    """
    pos in [-0.5, 0.5), corresponding to the image charge correction

    n(ξ) = 1/√(ε₁ε₂ε₃) √(p/π)³ exp(-p ξ²) 
        ξ    = √(r₁²/ε₁ + r₂²/ε₂ + r₃²/ε₃)
        p    = εₘₐₓ/(2σ²)
    
    ε is scaled by the `min(ε)`
    """
    Ω = abs(det(AL))
    n123 = np.array(mesh.shape)
    dΩ = Ω / prod(n123)

    # epsilon
    inv_eps = inv(epsilon)
    eps_ii = eigvalsh(epsilon)
    # |
    r = pos @ AL
    ξ, r = sqrt(r.T @ inv_eps @ r), norm(r)
    # |
    meps_ii = min(eps_ii)
    # |
    inv_eps = inv_eps * meps_ii
    eps_ii = eps_ii / meps_ii
    
    # indeces
    center = INTEGER_8(np.ceil(pos * n123))
    window_lb, window_ur = center - window // 2, center + window // 2
    mesh_indices = np.mgrid[window_lb[0]:window_ur[0], 
                            window_lb[1]:window_ur[1], 
                            window_lb[2]:window_ur[2]]  #
    flat_indices = mesh_indices.reshape(3, -1).T  # [m, 3]
    
    rT = (flat_indices / n123 - pos) @ AL  # [m, 3]
    ξ2 = einsum("ij,ij->i", rT @ inv_eps, rT)
    p = 1 / (2 * sigma ** 2)
    
    gaussian_vals = np.exp(- p * ξ2)
    
    flat_indices = flat_indices % n123 
    vr_sum = sum(mesh[flat_indices[:, 0], 
                      flat_indices[:, 1], 
                      flat_indices[:, 2]] * gaussian_vals)
    gaussian_sum = sum(gaussian_vals)
    
    vr_avg = vr_sum / gaussian_sum
    
    gaussian_sum *= 1/sqrt(prod(eps_ii)) * sqrt(p / pi)**3 * dΩ
    if (verbose := True) and not (0.9 < gaussian_sum < 1.0):
        print(f"{gaussian_sum:.6f}->(guassian_sum) {n123 = } \n"
              # f"\t{eps_ii=}, {inv_eps=}\n"
              f"\t{window_lb=}, {window_ur=} \n"
              f"\t{pos@AL=}")
        if debug := False:
            print(f"{pos=}")
            print(f"{rT[:]}")
            print(f"{flat_indices[:]}")
            print(f"{flat_indices/n123[:]}")
            print(f"{(flat_indices / n123 - pos)=}")
    
    return r, ξ, vr_avg


def write_diffvatom(filename: str, 
                    bulk: AtomConfig | VATOM, supcl: AtomConfig | VATOM, 
                    order, epsilon: NDArray, charge, 
                    r: NDArray, ξ: NDArray, 
                    dv_bulk: NDArray, dv_supcl0: NDArray, dv_supcl: NDArray):
    """
    natoms
    ...
    # bulk.ac.itypes, bulk.ac.positions, bulk.vatom, supcl.ac.itypes, supcl.ac.positions, supcl.vatom, diff_vatom,
    # metadata
    """
    if exists(filename):
        pass  # TODO overwrite or rename

    ξ[ξ < 1e-10] = 1e-10  # avoid divide by zero
    try:
        eps_ii = eigvalsh(epsilon)
    except np.linalg.LinAlgError:
        eps_ii = np.nan
    pointv      = charge / ξ / sqrt(prod(eps_ii)) 
    pointv[ξ <= 1e-10] = charge * np.inf

    r           = r * BOHR
    ξ           = ξ * BOHR
    pointv      = pointv * HA        
    dv_bulk     = dv_bulk * HA
    dv_supcl0   = dv_supcl0 * HA
    dv_supcl    = dv_supcl * HA
    diff_vatom0 = (dv_supcl0 - dv_bulk)
    diff_vatom  = (dv_supcl - dv_bulk)

    supcl.revise_atomsposition()
    bulk.revise_atomsposition()
    with open(filename, "w") as f:
        f.write(f" {supcl.natoms} atoms\n")
        f.write(" EPSILON\n")
        for i in range(3):
            f.write("%15.8f %15.8f %15.8f\n" % (
                epsilon[i, 0], epsilon[i, 1], epsilon[i, 2]))
        f.write(" diff_VATOM  # "
                "%29s %49s %15s %49s %31s %35s\n" % (
                    f"r, ξ, {charge}/√(ε₁ε₂ε₃)/ξ",
                    "bulk.ac.(itypes, positions[0:3])", 
                    "bulk.vatom",
                    "supcl.ac.(itypes, positions[0:3])", 
                    "supcl.[vatom, vatom_corrected]",
                    "diff_vatom, diff_vatom_corrected"))
        for i in range(supcl.natoms):
            k = order[i]
            f.write(f"{r[i]:14.8f} {ξ[i]:14.8f} {pointv[i]:14.8f} ")
            f.write("%10d %12.8f %12.8f %12.8f %15.8f " % (
                bulk.itypes[k], *bulk.positions[k], 
                dv_bulk[i]))
            f.write("%10d %12.8f %12.8f %12.8f %15.8f %15.8f " % (
                supcl.itypes[i], *supcl.positions[i], 
                dv_supcl0[i], dv_supcl[i]))
            f.write(f"{diff_vatom0[i]:20.8f} {diff_vatom[i]:14.8f}\n")
