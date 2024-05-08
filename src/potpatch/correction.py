from time import perf_counter

import numpy as np
from numba import jit, guvectorize

from potpatch.utils import simpson
from potpatch.constant import HA, EPSILON0
from potpatch.objects import Lattice, VR, AtomConfig, MaterialSystemInfo


def gen_charge_correct(supclInfo:MaterialSystemInfo):
    """
    supclInfo: provide AL & n123 & charge
        
    `charge` at `charge_density` & `plus_V`
    `epsilon` at `minus/plus_V`
    don't modify them
    
    in this function, all quantity is in a.u., 1 electron has 1 a.u. charge.

    generate functions
    """
    @jit(nopython=True)
    def charge_density(charge, r, R_cutoff):
        return np.sinc(r/R_cutoff) * (np.pi / (4*R_cutoff**3)) * charge if r<R_cutoff else 0

    AL      = supclInfo.lattice.AL_AU
    n123    = supclInfo.vr.n123
    charge  = supclInfo.charge
    epsilon = supclInfo.epsilon

    vol = np.abs(np.linalg.det(AL))
    rlattice = (2 * np.pi * np.linalg.inv(AL)).T
    Rmin = np.min( 
        [ np.dot(rlattice[i]/np.linalg.norm(rlattice[i]), AL[i]) for i in range(0,3) ]
    ) / 2

    # minus_V_periodic
    little_cube = np.array( [AL[i]/n123[i] for i in range(3)] )
    mesh_rho = np.ndarray(n123, dtype=np.float64)

    @jit(nopython=True)
    def _make_mesh_rho(mesh_rho):
        n1, n2, n3 = mesh_rho.shape
        for i in range(-n1//2,n1//2):
            for j in range(-n2//2,n2//2):
                for k in range(-n3//2,n3//2):
                    r = little_cube[0]*i + little_cube[1]*j + little_cube[2]*k
                    r = np.sqrt(r[0]**2 + r[1]**2 + r[2]**2)
                    mesh_rho[i,j,k] = charge_density(charge, r, Rmin)
    t0 = perf_counter()
    _make_mesh_rho(mesh_rho)
    # print(f"_make_mesh_rho time cost: {perf_counter() - t0} s")
    

    def minus_V_periodic(supclInfo:MaterialSystemInfo):
        """
        input supercell MSInfo
        modify its 3D-array vr mesh, 
        """
        fourier_coeff = np.fft.fftn(mesh_rho, axes=(0,1,2))

        @jit(nopython=True)
        def _rho2V_fourier_coeff(fourier_coeff):
            n1, n2, n3 = fourier_coeff.shape
            # n=5 -> -3...1 -> -2...2
            # n=6 -> -3...2 -> -2...0 ∪ 1...3
            for i in range(-n1//2+1,n1//2+1):
                for j in range(-n2//2+1,n2//2+1):
                    for k in range(-n3//2+1,n3//2+1):
                        G = rlattice[0]*i + rlattice[1]*j + rlattice[2]*k
                        G2 = G[0]*G[0] + G[1]*G[1] + G[2]*G[2]
                        if G2 == 0:
                            fourier_coeff[0,0,0] = 0
                        else:
                            fourier_coeff[i,j,k] = -1 * -fourier_coeff[i,j,k] / G2 / (EPSILON0*epsilon)
        t0 = perf_counter()
        _rho2V_fourier_coeff(fourier_coeff)
        # print(f"_rho2V_fourier_coeff time cost: {perf_counter() - t0} s")

        mesh_V = np.fft.ifftn(fourier_coeff, axes=(0,1,2))
        mesh_V = np.real(mesh_V)
        assert np.sum(np.abs(np.imag(mesh_V))) < 1e-10, "V_mesh from iFFT is not real"
        supclInfo.vr.mesh -= mesh_V

    # plus_V_single
    R_cutoff = Rmin # TODO R_cutoff 为了那些有尾巴的函数准备; Rmin 切于 unit cell 的球半径, 作为一个参考
    ## n = 1000 # TODO 把这个变成一个可控参数
    ## # charge solid sphere
    ## charge_shll = np.ndarray(4*n+1)
    ## for i in range(len(charge_shll)):
    ##     r = i*R_cutoff/(4*n)
    ##     charge_shll[i] = charge_density(charge, r, Rmin) * 4 * np.pi * r**2
    ## charge_sph = np.cumsum(simpson(R_cutoff/(4*n), charge_shll)) 
    ## # electric field intensity
    ## r = np.linspace(0,R_cutoff,num=len(charge_sph))
    ## ele_fld = np.ndarray(len(charge_sph))
    ## ele_fld[0] = 0
    ## ele_fld[1:] = charge_sph[1:] / np.square(r[1:]) / epsilon # in a.u., 1/4πε_0==1
    ## # potential
    ## V_sph = np.cumsum(simpson(R_cutoff/(2*n), ele_fld[::-1]))
    ## V_sph = V_sph[::-1] 
    ## V_sph = V_sph + 1/R_cutoff * charge / epsilon

    def plus_V_single(suuuupclInfo:MaterialSystemInfo):
        """
        input suuuupercell MSInfo
        modify its 3D-array vr mesh, 
        """
        @jit(nopython=True)
        def _plus_V_single(suuuupcl_mesh):
            n1, n2, n3 = suuuupcl_mesh.shape
            for i in range(-n1//2,n1//2):
                for j in range(-n2//2,n2//2):
                    for k in range(-n3//2,n3//2):
                        r = little_cube[0]*i + little_cube[1]*j + little_cube[2]*k
                        r = np.sqrt(r[0]**2 + r[1]**2 + r[2]**2)
                        if r > R_cutoff:
                            suuuupcl_mesh[i,j,k] += 1/r * charge / epsilon
                        else:
                            ## dr = R_cutoff / n
                            ## i = int(r/dr)
                            ## (V_sph[i]*(r-i*dr)+V_sph[i+1]*((i+1)*dr-r)) / dr
                            suuuupcl_mesh[i,j,k] += \
                                (np.sinc(r/R_cutoff)/R_cutoff + 1/R_cutoff) * charge  / epsilon
        t0 = perf_counter()
        _plus_V_single(suuuupclInfo.vr.mesh)
        # print(f"_plus_V_single time cost: {perf_counter() - t0} s")

    return minus_V_periodic, plus_V_single

# TODO 如果奇数超胞, 奇数扩胞, 会出现什么问题吗
# TODO 如果是奇数网格怎么办
def edge_match_correct(supclInfo:MaterialSystemInfo, bulkInfo:MaterialSystemInfo):
    supclmesh = supclInfo.vr.mesh
    bulkmesh = bulkInfo.vr.mesh
    n1, n2, n3 = supclmesh.shape
    n1b, n2b, n3b = bulkmesh.shape
    num = n1*n2+n2*n3+n3*n1-(n1+n2+n3)+1
    diff = np.ndarray(num, dtype=np.float64)
    n = 0
    for i in range(n1):
        for j in range(n2):
            for k in range(n3):
                if i==n1//2 or j==n2//2 or k==n3//2 :
                    diff[n] = supclmesh[i,j,k] - bulkmesh[i%n1b, j%n2b, k%n3b]
                    n += 1
    while True:
        sigma = np.std(diff)
        diffmean = np.mean(diff)
        mask = np.array( [ np.abs( i - diffmean ) < 3*sigma for i in diff ] )
        if diff.size <= 0 or np.all(mask) : break
        diff = diff[mask]
    edge_shift, sigma = np.mean(diff), np.std(diff)
    print(f"edge shift (V_align): {edge_shift*HA*1000} meV, {num-len(diff)}/{num} data detached")
    print(f"standard deviation of diffs at the boundary: {sigma*HA*1000} meV")
    supclInfo.vr.mesh -= edge_shift

