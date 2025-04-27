from itertools import product
from typing import overload

import numpy as np
from numpy import pi
from numpy import prod, abs, sqrt
from numpy.linalg import norm, inv, det, eigvalsh
from numba import jit, guvectorize, prange

from potpatch.utils import simpson, timing, NameTuple, revise_epsilon
from potpatch.constant import HA, EPSILON0, BOHR
from potpatch.datatype import REAL_8
from potpatch.objects import Lattice, VR, AtomConfig, MaterialSystemInfo


def gen_charge_correct(supclInfo: MaterialSystemInfo, 
                       correction: NameTuple):
    """
    supclInfo: provide AL & n123 & charge
        
    `charge` at `charge_density` & `plus_V`
    `charge_pos` in (-0.5,0.5] (because of -n//2+1...n//2+1); zeros(3) if is None; at decive `Rmin`
    `epsilon` at `minus/plus_V`
    don't modify them
    
    in this function, all quantity is in a.u., 1 electron has 1 a.u. charge.

    generate functions
    """
    @jit(nopython=True)
    def charge_density(charge, r, R_cutoff):
        return np.sinc(r/R_cutoff) * (np.pi / (4*R_cutoff**3)) * charge \
            if r < R_cutoff else 0

    AL      = supclInfo.lattice.AL_AU
    n123    = supclInfo.vr.n123
    charge  = correction.charge
    charge_pos = supclInfo.charge_pos if supclInfo.charge_pos is not None else np.zeros(3)
    epsilon = correction.epsilon
    epsilon = revise_epsilon(epsilon)

    for i in range(3):
        pos_i = charge_pos[i] % 1
        charge_pos[i] = pos_i if pos_i <= 0.5 else pos_i-1  # n=6  //2->  -3...2  +1->  -2...0 ∪ 1...3( `=` matter)
    for i, j in product(range(3), range(3)):
        assert epsilon[i, j] == epsilon[j, i], "epsilon is not symmetric"
    eps_ii = eigvalsh(epsilon)
    assert len(eps_ii) == 3, "epsilon is not a 3x3 matrix"
    inv_eps = inv(epsilon)

    vol = abs(det(AL))
    rlattice = (2 * pi * inv(AL)).T

    R = chargemodel_R(AL, epsilon, charge_pos)
    _n = (4/3*np.pi * R**3) / (vol / prod(n123)) * sqrt(prod(eps_ii))
    _n = int(round(_n))
    print(
        f"gen_charge_correct: charge_model.radius = "
        f"{R*BOHR*max(sqrt(eps_ii)):.3f} angstrom\n"
        f"    about {_n}/{prod(n123)} grid points in it. "
    )

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # minus_V_periodic
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    delta_AL = np.array([AL[i]/n123[i] for i in range(3)])
    mesh_rho = np.zeros(n123, dtype=REAL_8)

    @timing()
    @jit(nopython=True)
    def _make_mesh_rho(mesh_rho, n123):
        """a subtle implementation

        mesh_rho is a *zero* array for skipping; 
        mesh_rho can be indexed with negative numbers 
        """
        # `AL`, `charge_pos` are in supclInfo
        n1, n2, n3 = n123
        for i in prange(-n1//2+1, n1//2+1):
            for j in prange(-n2//2+1, n2//2+1):
                for k in prange(-n3//2+1, n3//2+1):
                    r = delta_AL[0] * i - AL[0] * charge_pos[0] + \
                        delta_AL[1] * j - AL[1] * charge_pos[1] + \
                        delta_AL[2] * k - AL[2] * charge_pos[2]
                    r = r[0] * inv_eps[0, 0] * r[0] + \
                        r[0] * inv_eps[0, 1] * r[1] + \
                        r[0] * inv_eps[0, 2] * r[2] + \
                        r[1] * inv_eps[1, 0] * r[0] + \
                        r[1] * inv_eps[1, 1] * r[1] + \
                        r[1] * inv_eps[1, 2] * r[2] + \
                        r[2] * inv_eps[2, 0] * r[0] + \
                        r[2] * inv_eps[2, 1] * r[1] + \
                        r[2] * inv_eps[2, 2] * r[2]
                    r = np.sqrt(r)
                    mesh_rho[i, j, k] = charge_density(charge, r, R)
    _make_mesh_rho(mesh_rho, mesh_rho.shape)
    
    if (debug := False):
        VR(lattice=supclInfo.lattice, mesh=mesh_rho).write_vr("model_charge")
    
    fourier_coeff = np.fft.rfftn(mesh_rho, axes=(0, 1, 2))

    @timing()
    @jit(nopython=True)
    def _rho2V_fourier_coeff(fourier_coeff, rlattice):
        n1, n2, n3_r = fourier_coeff.shape

        eps_prod = np.prod(eps_ii)
        eps_sqrt = np.sqrt(eps_prod)

        r0 = rlattice[0]
        r1 = rlattice[1]
        r2 = rlattice[2]
        e0 = epsilon[0]
        e1 = epsilon[1]
        e2 = epsilon[2]
        
        # n=5  //2->  -3...1  +1->  -2...2 (matter)
        # n=6  //2->  -3...2  +1->  -2...0 ∪ 1...3
        for i in prange(-n1//2+1, n1//2+1):
            for j in prange(-n2//2+1, n2//2+1):
                # n3=5  n3_r=3  (0..2 ~ -2...2)
                # n3=6  n3_r=4  (0..3 ~ -2...3)
                for k in prange(0, n3_r):
                    Gx = i*r0[0] + j*r1[0] + k*r2[0]
                    Gy = i*r0[1] + j*r1[1] + k*r2[1]
                    Gz = i*r0[2] + j*r1[2] + k*r2[2]

                    eGx = e0[0]*Gx + e0[1]*Gy + e0[2]*Gz
                    eGy = e1[0]*Gx + e1[1]*Gy + e1[2]*Gz
                    eGz = e2[0]*Gx + e2[1]*Gy + e2[2]*Gz
                    GeG = Gx*eGx + Gy*eGy + Gz*eGz

                    if GeG == 0:
                        fourier_coeff[0, 0, 0] = 0
                    else:
                        fourier_coeff[i, j, k] = (
                            -1                         # ∇⋅ε∇V(r)=-ρ(r）
                            * -fourier_coeff[i, j, k]  # exp(iG⋅r)
                            / GeG / (EPSILON0)) / eps_sqrt
    _rho2V_fourier_coeff(fourier_coeff, rlattice)

    mesh_V = np.fft.irfftn(fourier_coeff, axes=(0, 1, 2))
    mesh_rho = None
    fourier_coeff = None
    
    def minus_V_periodic(supclInfo: MaterialSystemInfo, inverse=False):
        """
        input supercell MSInfo
        modify its 3D-array vr mesh, 
        TODO input VR but not MatInfo
        """
        if inverse:
            supclInfo.vr.mesh += mesh_V
        else:
            supclInfo.vr.mesh -= mesh_V

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # plus_V_single
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    R_cutoff = R 

    if correction.plus_V == "single":
        @timing()
        @jit(nopython=True)
        def _plus_V_single(suuuupcl_mesh, inverse=1.0):
            # `AL`, `charge_pos` are in supclInfo
            n1, n2, n3 = suuuupcl_mesh.shape
            for i in prange(-n1//2+1, n1//2+1):
                for j in prange(-n2//2+1, n2//2+1):
                    for k in prange(-n3//2+1, n3//2+1):
                        r = delta_AL[0] * i - AL[0] * charge_pos[0] + \
                            delta_AL[1] * j - AL[1] * charge_pos[1] + \
                            delta_AL[2] * k - AL[2] * charge_pos[2]
                        r = r[0] * inv_eps[0, 0] * r[0] + \
                            r[0] * inv_eps[0, 1] * r[1] + \
                            r[0] * inv_eps[0, 2] * r[2] + \
                            r[1] * inv_eps[1, 0] * r[0] + \
                            r[1] * inv_eps[1, 1] * r[1] + \
                            r[1] * inv_eps[1, 2] * r[2] + \
                            r[2] * inv_eps[2, 0] * r[0] + \
                            r[2] * inv_eps[2, 1] * r[1] + \
                            r[2] * inv_eps[2, 2] * r[2]
                        r = np.sqrt(r)
                        if r > R_cutoff:
                            v = 1/r * charge / sqrt(prod(eps_ii)) 
                        else:
                            v = (np.sinc(r/R_cutoff)/R_cutoff + 1/R_cutoff) \
                                * charge / sqrt(prod(eps_ii)) 
                        suuuupcl_mesh[i, j, k] += v * inverse

        def plus_V_single(suuuupclInfo: MaterialSystemInfo, inverse=False):
            """
            input suuuupercell MSInfo
            modify its 3D-array vr mesh, 
            TODO input VR but not MatInfo
            """
            if inverse:
                _plus_V_single(suuuupclInfo.vr.mesh, inverse=-1.0)
            else:
                _plus_V_single(suuuupclInfo.vr.mesh)

        plus_V = plus_V_single
    
    elif correction.plus_V == "periodic":
        def plus_V_periodic(suuuupclInfo: MaterialSystemInfo, inverse=False):
            # _make_mesh_rho, AL, delta_AL, R, charge 
            mesh_rho = np.zeros_like(suuuupclInfo.vr.mesh)
            _make_mesh_rho(mesh_rho, n123=supclInfo.vr.n123)

            rlattice = (2 * np.pi * inv(suuuupclInfo.lattice.AL_AU)).T
            fourier_coeff = np.fft.rfftn(mesh_rho, axes=(0, 1, 2))
            mesh_rho = None
            _rho2V_fourier_coeff(fourier_coeff, rlattice)
            fourier_coeff = None
            mesh2_V = np.fft.irfftn(fourier_coeff, axes=(0, 1, 2))

            if inverse:
                suuuupclInfo.vr.mesh -= mesh2_V
            else:
                suuuupclInfo.vr.mesh += mesh2_V

        plus_V = plus_V_periodic

    else:
        assert False, 'correction.plus_V is not "single" or "periodic"'

    return minus_V_periodic, plus_V


def gen_charge_correct_gaussian(supclInfo: MaterialSystemInfo):
    raise NotImplementedError()

    def minus_V_periodic(supclInfo: MaterialSystemInfo):
        def _make_V_inG(epsilon):
            pass

        pass

    def plus_V_single(suuuupclInfo: MaterialSystemInfo):
        """
        ε₀ * ∇⋅ε∇ V(ξ) = - n(ξ), where
            V(ξ) = 1/√(ε₁ε₂ε₃) 1/ξ erf[√p ξ] 
            n(ξ) = 1/√(ε₁ε₂ε₃) √(p/π)³ exp(-p ξ²) 
            ξ    = √(r₁²/ε₁ + r₂²/ε₂ + r₃²/ε₃)
        """
        pass

    return minus_V_periodic, plus_V_single


def chargemodel_R(AL, epsilon, charge_pos=np.zeros(3)):
    rlattice = (2 * pi * inv(AL)).T
    perp_nvecs = [rlattice[i]/norm(rlattice[i]) for i in range(3)]
    heights = [perp_nvecs[i].T @ AL[i] for i in range(3)]
    charge_pos = (np.array(charge_pos) + 0.5) % 1
    inv_eps = inv(epsilon)

    R = np.zeros(3)
    for i in range(3):
        p, L = perp_nvecs[i], heights[i]
        alpha = min(abs(0-charge_pos[i]), abs(1-charge_pos[i]))
        cos_theta = (p.T @ epsilon @ p) / (norm(p) * norm(epsilon @ p))
        R[i] = alpha * L / cos_theta * sqrt(p.T @ inv_eps @ p)

    return min(R)


# TODO 如果奇数超胞, 奇数扩胞, 会出现什么问题吗
# TODO 如果是奇数网格怎么办
def edge_match_correct(supclInfo: MaterialSystemInfo, bulkInfo: MaterialSystemInfo):
    """
        TODO input VR but not MatInfo
    """
    supclmesh     = supclInfo.vr.mesh
    bulkmesh      = bulkInfo.vr.mesh

    diff, num, num_left = edge_diff(supclmesh, bulkmesh, (0, 1, 2))
    edge_shift, sigma = mean_std(diff)

    print(f"edge shift (V_align): {edge_shift*HA*1000:6.2f} meV, "
          f"{num-num_left}/{num} data detached")
    print(f"standard deviation of diffs at the boundary: "
          f"{sigma*HA*1000:6.2f} meV")

    if (verbose := True):
        diff_yz, m1, ml1 = edge_diff(supclmesh, bulkmesh, (0,))
        diff_xz, m2, ml2 = edge_diff(supclmesh, bulkmesh, (1,))
        diff_xy, m3, ml3 = edge_diff(supclmesh, bulkmesh, (2,))
        mean_yz, std_yz = mean_std(diff_yz)
        mean_xz, std_xz = mean_std(diff_xz)
        mean_xy, std_xy = mean_std(diff_xy)
        print(f"    mean_yz: {mean_yz*HA*1000:6.2f} meV, "
              f"std_yz: {std_yz*HA*1000:6.2f} meV, "
              f"({m1-ml1}/{m1}) data detached")
        print(f"    mean_xz: {mean_xz*HA*1000:6.2f} meV, "
              f"std_xz: {std_xz*HA*1000:6.2f} meV, "
              f"({m2-ml2}/{m2}) data detached")
        print(f"    mean_xy: {mean_xy*HA*1000:6.2f} meV, "
              f"std_xy: {std_xy*HA*1000:6.2f} meV, "
              f"({m3-ml3}/{m3}) data detached")

    supclInfo.vr.mesh -= edge_shift


@timing()
@jit(nopython=True)
def edge_diff(supclmesh, bulkmesh, axes: tuple):
    """
    if 0 in `axes`, then include the yz plane in `diff`
    """
    n1,  n2,  n3  = supclmesh.shape
    n1b, n2b, n3b = bulkmesh.shape

    num = sum(
        [(n2 * n3, n1 * n3, n1 * n2)[ax] for ax in axes]
        ) * 2

    n, diff = 0, np.zeros(num, dtype=REAL_8)
    if 0 in axes:
        for i in range(-n2//2+1, n2//2+1):
            for j in range(-n3//2+1, n3//2+1):
                diff[n] = supclmesh[n1//2, i, j] - \
                    bulkmesh[(n1//2) % n1b, i % n2b, j % n3b]
                diff[n+1] = supclmesh[-(n1//2), i, j] - \
                    bulkmesh[(-(n1//2)) % n1b, i % n2b, j % n3b]
                n += 2
    if 1 in axes:
        for i in range(-n1//2+1, n1//2+1):
            for j in range(-n3//2+1, n3//2+1):
                diff[n] = supclmesh[i, n2//2, j] - \
                    bulkmesh[i % n1b, (n2//2) % n2b, j % n3b]
                diff[n+1] = supclmesh[i, -(n2//2), j] - \
                    bulkmesh[i % n1b, (-(n2//2)) % n2b, j % n3b]
                n += 2
    if 2 in axes:
        for i in range(-n1//2+1, n1//2+1):
            for j in range(-n2//2+1, n2//2+1):
                diff[n] = supclmesh[i, j, n3//2] - \
                    bulkmesh[i % n1b, j % n2b, (n3//2) % n3b]
                diff[n+1] = supclmesh[i, j, -(n3//2)] - \
                    bulkmesh[i % n1b, j % n2b, (-(n3//2)) % n3b]
                n += 2
    assert n == num

    while True:
        sigma = np.std(diff)
        diffmean = np.mean(diff)
        mask = np.array([np.abs(i - diffmean) < 3*sigma for i in diff])
        if diff.size <= 0 or np.all(mask): 
            break
        diff = diff[mask]
    num_left = len(diff)

    return diff, num, num_left


def mean_std(diff):
    diff = np.array(diff)
    return np.mean(diff), np.std(diff)
