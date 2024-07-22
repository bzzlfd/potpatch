import warnings
from textwrap import dedent
from time import perf_counter

import numpy as np
from numba import jit, guvectorize

from potpatch.objects import Lattice, VR, AtomConfig, MaterialSystemInfo
from potpatch.supercell import make_supercell, modify_supercell, closed_to_edge

def inspect_ingredient(supclInfo:MaterialSystemInfo, bulkInfo:MaterialSystemInfo, size_confirm=None, frozen_confirm=None):
    """
    size_confirm: optional parameter `potpatch.supercell.size` in "potpatch.input"
                the `supercell_size` information is infered from `atom.config` and `VR` files, so
                it is not necessary to pass the `supercell_size` to `patch`
                this parameter was designed to confirm `supercell_size` as user expect
    frozen_confirm: optional parameter `potpatch.supercell.frozen_range` in "potpatch.input"
                when open it, this function will check whether the atom positions in frozen range 
                are coincide of bulk and supercell
    """
    
    # supcl size inference
    supcl_vrsize = supclInfo.vr.n123 / bulkInfo.vr.n123 
    # supcl size ?integer mag
    if not all(np.abs(supcl_vrsize - np.int32(supcl_vrsize)) < 1e-6):
        raise ValueError("magnifacation between the two VR is not integer")
    # Lattice and VR.n123: ?same
    lattice_mulmag = bulkInfo.lattice * supcl_vrsize
    if not supclInfo.lattice == lattice_mulmag:
        raise ValueError(f"""
                         magnifacation between Lattice and VR_mesh is not equal
                         supclInfo.lattice({supclInfo.lattice.in_unit("angstrom")})
                         bulkInfo.lattice * mag({lattice_mulmag.in_unit("angstrom")})
                         """)
    
    supcl_size = np.int32(supcl_vrsize)

    if not all( bulkInfo.vr.n123 % 2 == 0 ) or \
       not all( supclInfo.vr.n123 % 2 == 0 ) :
        warnings.warn(dedent(f"""
                    `N123` in bulk.vr is {bulkInfo.vr.n123}
                    `N123` in supercll.vr is {supclInfo.vr.n123}
                    I am assuming that I'm dealing with even-numbered mesh for VR 
                    when I write this code. if your grid is odd, it could lead to 
                    potential issues. 
                    It's recommended that forcibly set `N123` to the even numver.
                    """))

    if supclInfo.charge is None:
        raise ValueError("information of supercell `charge` is not given")
    if supclInfo.epsilon is None:
        raise ValueError("information of supercell `epsilon` is not given")

    if size_confirm is not None:
        print(f"{size_confirm = } is specified...", end="")
        if not all(supcl_size == size_confirm):
            warnings.warn(f"supcl_size = {list(supcl_size)} and size_confirm = {list(size_confirm)}")
        else:
            print("no error")
    if frozen_confirm is not None:
        print(f"{frozen_confirm = } is specified...", end="")
        supclAtom_1 = make_supercell(bulkInfo.atomconfig, supcl_size)
        supclAtom_1 = modify_supercell(supclAtom_1, frozen_confirm)
        supclAtom_2 = supclInfo.atomconfig
        err = False
        for pos1 in supclAtom_1.positions:
            if closed_to_edge(supclAtom_1.lattice, pos1, frozen_confirm):
                if min(sum(np.abs(pos1-pos2)) for pos2 in supclAtom_2.positions) != 0:
                    err = True
                    l2 = [ sum(np.abs(pos1-pos2)) for pos2 in supclAtom_2.positions ] 
                    pos2_index = l2.index(sorted(l2)[0])
                    pos2 = supclAtom_2.positions[pos2_index]
                    warnings.warn(f"{pos1}(from bulk make_supercell) and {pos2}(from supcl) doesn't coincide ")
        if not err:
            print("no error")
    return supcl_size


def patch(supclInfo:MaterialSystemInfo, bulkInfo:MaterialSystemInfo, supcl_size, target_size) -> MaterialSystemInfo:
    suuuupclInfo = MaterialSystemInfo()
    suuuupclInfo.charge=supclInfo.charge

    t0 = perf_counter()
    suuuupclInfo.atomconfig = patch_atom_v2(supclInfo.atomconfig, bulkInfo.atomconfig, supcl_size, target_size)
    # print(f"patch atom time cost: {perf_counter() - t0} s")
    t0 = perf_counter()
    suuuupclInfo.vr = patch_vr(supclInfo.vr, bulkInfo.vr, supcl_size, target_size)
    # print(f"patch vr time cost: {perf_counter() - t0} s")
    return suuuupclInfo


# 这个代码只适合偶数 target_size
def patch_atom(supclAtom:AtomConfig, bulkAtom:AtomConfig, supcl_size, target_size) -> AtomConfig:
    """
    resolution is bulk, supcl_size should be even numbers
    """
    supclAtom.revise_atomsposition()
    bulkAtom.revise_atomsposition()
    
    indx = supclAtom.positions >= 0.5
    supclAtom.positions[indx] -= 1.0
    # code prototype from AtomConfig.__mul__
    nrepeat         = np.prod(target_size)
    natoms          = bulkAtom.natoms * nrepeat
    lattice         = bulkAtom.lattice * target_size
    atoms_itype     = np.tile(bulkAtom.itypes, nrepeat)
    atoms_position  = np.tile(bulkAtom.positions, (nrepeat,1))
    atoms_move      = np.tile(bulkAtom.moves, (nrepeat,1))
    
    # supcl atoms
    atoms_itype[0:supclAtom.natoms]     = supclAtom.itypes
    atoms_position[0:supclAtom.natoms]  = supclAtom.positions * supcl_size
    atoms_move[0:supclAtom.natoms]      = supclAtom.moves

    # bulk atoms
    cnt = 0
    for i in range(-target_size[0]//2, target_size[0]//2): # 5 -3...1; 4 -2...1
        for j in range(-target_size[1]//2,target_size[1]//2):
            for k in range(-target_size[2], target_size[2]//2):
                if not ((-supcl_size[0]//2 <= i < supcl_size[0]//2) or\
                        (-supcl_size[1]//2 <= j < supcl_size[1]//2) or\
                        (-supcl_size[2]//2 <= k < supcl_size[2]//2)):
                    indx = bulkAtom.natoms * (cnt + np.prod(supcl_size))
                    atoms_position[indx:indx+bulkAtom.natoms] += np.array([i,j,k])
                    atoms_itype[indx:indx+bulkAtom.natoms] = bulkAtom.itypes
                    atoms_move[indx:indx+bulkAtom.natoms] = bulkAtom.moves
                    cnt += 1
    
    atoms_position /= target_size
    rtac = AtomConfig(natoms=natoms, lattice=lattice, \
                    itypes=atoms_itype, positions=atoms_position, moves=atoms_move)
    rtac.revise_atomsposition()
    return rtac

# TODO 如果奇数超胞, 奇数扩胞, 会出现什么问题吗
# TODO 如果是奇数网格怎么办
def patch_atom_v2(supclAtom:AtomConfig, bulkAtom:AtomConfig, supcl_size, target_size) -> AtomConfig:
    """
    another implement of patch_atom, break the bulkcell-resolution. 
    but hidden issues may arise. 
    """
    supclAtom.revise_atomsposition()
    bulkAtom.revise_atomsposition()
    
    indx = supclAtom.positions >= 0.5
    supclAtom.positions[indx] -= 1.0
    # code prototype from AtomConfig.__mul__
    nrepeat         = np.prod(target_size)
    natoms          = bulkAtom.natoms * nrepeat
    lattice         = bulkAtom.lattice * target_size
    atoms_itype     = np.tile(bulkAtom.itypes, nrepeat)
    atoms_position  = np.tile(bulkAtom.positions, (nrepeat,1))
    atoms_move      = np.tile(bulkAtom.moves, (nrepeat,1))

    # supcl atoms
    atoms_itype[0:supclAtom.natoms]     = supclAtom.itypes
    atoms_position[0:supclAtom.natoms]  = supclAtom.positions * supcl_size
    atoms_move[0:supclAtom.natoms]      = supclAtom.moves

    # bulk atoms
    cnt, offset = 0, supclAtom.natoms
    def in_imbox(pos, planl, planr):
        if all( planl <= pos ) and all( pos < planr ): return True
        return False

    target_size = np.array(target_size)
    planl, planr = -supcl_size/2, +supcl_size/2
    for i in range(target_size[0]):
        for j in range(target_size[1]):
            for k in range(target_size[2]):
                for l in range(bulkAtom.natoms):
                    ijk = np.array((i,j,k))
                    pos = bulkAtom.positions[l] + ijk
                    pos[pos >= target_size/2] -= target_size[pos >= target_size/2]
                    if not in_imbox(pos, planl, planr):
                        indx = cnt + offset
                        atoms_position[indx]    = pos
                        atoms_itype[indx]       = bulkAtom.itypes[l]
                        atoms_move[indx]        = bulkAtom.moves[l]
                        cnt += 1
    assert cnt + offset == natoms, "some thing wrong when setting atoms at patch_atom_v2"

    atoms_position /= target_size
    rtac = AtomConfig(natoms=natoms, lattice=lattice, \
                    itypes=atoms_itype, positions=atoms_position, moves=atoms_move)
    rtac.revise_atomsposition()
    return rtac


def patch_vr(supclVR:VR, bulkVR:VR, supcl_size, target_size) -> VR:
    suuuupclVR = bulkVR * target_size

    # @jit(nopython=True) 
    @guvectorize('void(float64[:,:,:], float64[:,:,:])','(n1,n2,n3),(m1,m2,m3)')
    def _overwrite_supcl_mesh(supcl_mesh, suuuupcl_mesh):
        n1, n2, n3 = supcl_mesh.shape
        for i in range(-n1//2,n1//2):
            for j in range(-n2//2,n2//2):
                for k in range(-n3//2,n3//2):
                    suuuupcl_mesh[i, j, k] = supcl_mesh[i,j,k]
    t0 = perf_counter()
    _overwrite_supcl_mesh(supclVR.mesh, suuuupclVR.mesh)
    # print(f"_overwrite_supcl_mesh time cost: {perf_counter() - t0} s")
    # print(_overwrite_supcl_mesh.inspect_types())
    return suuuupclVR
