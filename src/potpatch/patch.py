import warnings
from textwrap import indent, dedent

import numpy as np
from numba import jit, guvectorize

from potpatch.objects import Lattice, VR, AtomConfig, MaterialSystemInfo
from potpatch.supercell import (make_supercell, modify_supercell, closed_to_edge,
                                infer_supercell_size)
from potpatch.atompos_coin import check_atompos_consistency
from potpatch.utils import timing, log
from potpatch.datatype import REAL_8, INTEGER
from potpatch.validation import PATCH_REQUIRED_FIELDS


def _resample_periodic_mesh(mesh: np.ndarray,
                            shape: tuple[int, int, int]) -> np.ndarray:
    """Fourier-resample a real periodic mesh onto ``shape`` grid points.

    Samples remain registered at fractional coordinate zero on every axis.
    The zero-frequency coefficient, and hence the mean potential, is kept.
    """
    result = mesh
    for axis, target in enumerate(shape):
        source = result.shape[axis]
        if source == target:
            continue

        spectrum = np.fft.rfft(result, axis=axis)
        new_shape = list(spectrum.shape)
        new_shape[axis] = target // 2 + 1
        resized = np.zeros(new_shape, dtype=spectrum.dtype)
        common = min(source, target) // 2 + 1
        sl = [slice(None)] * mesh.ndim
        sl[axis] = slice(0, common)
        resized[tuple(sl)] = spectrum[tuple(sl)]

        # An even grid stores its Nyquist pair in one real coefficient.
        # Split that pair when enlarging and combine it when shrinking.
        nyquist = [slice(None)] * mesh.ndim
        if source % 2 == 0 and target > source:
            nyquist[axis] = source // 2
            resized[tuple(nyquist)] *= 0.5
        elif target % 2 == 0 and source > target:
            nyquist[axis] = target // 2
            resized[tuple(nyquist)] = 2 * resized[tuple(nyquist)].real

        result = np.fft.irfft(resized, n=target, axis=axis) * target / source
    return result


def inspect_ingredient(supclInfo: MaterialSystemInfo,
                       bulkInfo: MaterialSystemInfo,
                       size_confirm=None, frozen_confirm=None):
    log(f"{inspect_ingredient.__name__}")
    """
    size_confirm: optional parameter `potpatch.supercell.size` in "potpatch.input"
                the `supercell_size` information is inferred from the bulk and
                supercell lattices, so
                it is not necessary to pass the `supercell_size` to `patch`
                this parameter was designed to confirm `supercell_size` as user expect
    frozen_confirm: optional parameter `potpatch.supercell.frozen_range` in "potpatch.input"
                when open it, this function will check whether the atom positions in frozen range 
                are coincide of bulk and supercell
    """

    required = ("lattice", "atomconfig.lattice", "vr.lattice", "vr.mesh")
    if frozen_confirm is not None:
        required += ("atomconfig.natoms", "atomconfig.positions")
    for info in (bulkInfo, supclInfo):
        info.validate(required_paths=required).raise_for_errors()

    # Infer the physical supercell magnification from the lattices, then
    # check the VR mesh ratio required by patch_vr().
    log("infer supercell size from lattice")  # >log
    supcl_size, supcl_lattice_size = infer_supercell_size(
        bulkInfo.lattice, supclInfo.lattice)
    if not all(np.abs(supcl_lattice_size - supcl_size) < 1e-6):
        raise ValueError(
            "magnification between the bulk and supercell lattices is not integer\n"
            "bulkInfo.lattice:\n"
            f"{np.array2string(bulkInfo.lattice.in_unit('angstrom'), prefix='    ')}\n"
            "supclInfo.lattice:\n"
            f"{np.array2string(supclInfo.lattice.in_unit('angstrom'), prefix='    ')}\n"
            f"magnification: {supcl_lattice_size}"
            )

    # Keep the bulk grid as the output resolution.  Align a slightly different
    # supercell grid before charge correction, boundary matching, and patching.
    expected_n123 = bulkInfo.vr.n123 * supcl_size
    actual_n123 = supclInfo.vr.n123
    spacing_error = np.abs(expected_n123 / actual_n123 - 1)
    if np.any(spacing_error > 0.05):
        message = (
            "VR grid spacing mismatch exceeds 5% on an axis; cannot safely patch: "
            f"bulk N123={bulkInfo.vr.n123}, lattice size={supcl_size}, "
            f"expected supercell N123={expected_n123}, "
            f"actual supercell N123={actual_n123}, "
            f"spacing difference={np.round(spacing_error * 100, 2)}%"
        )
        warnings.warn(message, stacklevel=2)
        raise ValueError(message)
    if np.any(actual_n123 != expected_n123):
        log(f"resample supercell VR from N123={actual_n123} "
            f"to N123={expected_n123}")
        supclInfo.vr.mesh = _resample_periodic_mesh(
            supclInfo.vr.mesh, tuple(int(n) for n in expected_n123))

    if not all(bulkInfo.vr.n123 % 2 == 0) or \
       not all(supclInfo.vr.n123 % 2 == 0):
        warnings.warn(dedent(f"""
            `N123` in bulk.vr is {bulkInfo.vr.n123}
            `N123` in supercll.vr is {supclInfo.vr.n123}
            I am assuming that I'm dealing with even-numbered mesh for VR
            when I write this code. if your grid is odd, it could lead to
            potential issues.
            It's recommended that forcibly set `N123` to the even numver.
            """))

    if supclInfo.charge_pos is None:
        supclInfo.charge_pos = np.array([0., 0., 0.])
    else:
        if any(supclInfo.charge_pos != np.zeros(3)):
            for i in range(3):
                pos_i = supclInfo.charge_pos[i] % 1
                supclInfo.charge_pos[i] = pos_i if pos_i < 0.5 else pos_i - 1
            print(dedent(f"""
                Warning: supercell.charge_pos is {supclInfo.charge_pos}
                    In principle, it should be [0, 0, 0].
                    Parameter `charge_pos` provides more flexiblity, but with that comes limited/dangerous.
                    In previous supercell SCF calculation, patch boundary atoms were fixed. If
                    `charge_pos` is too close to the boundary, full relaxation is questionable.
                    Do your own check.
                """))

    if size_confirm is not None:       
        log(f"{size_confirm = } is specified ... ", end="")  # >log
        if not all(supcl_size == size_confirm):
            warnings.warn(f"supcl_size = {list(supcl_size)} and size_confirm = {list(size_confirm)}")
            log("", unadorned=True)
        else:
            log("no error", unadorned=True)
    if (frozen_confirm is not None):
        log(f"{frozen_confirm = } is specified ... ", end="")  # >log
        nwarn, max_dist = check_atompos_consistency(
                                             bulkInfo.atomconfig, 
                                             supclInfo.atomconfig, 
                                             frozen_range=frozen_confirm, 
                                             supcl_size=supcl_size)
        if nwarn == 0:
            log("no error", unadorned=True)
        elif nwarn > 0:
            log(f"{nwarn} warnings:  max_distance = {max_dist} angstrom", unadorned=True)
    return supcl_size


def patch(supclInfo: MaterialSystemInfo, bulkInfo: MaterialSystemInfo,
          supcl_size, target_size) -> MaterialSystemInfo:
    log(f"{patch.__name__}")
    for info in (bulkInfo, supclInfo):
        info.validate(required_paths=PATCH_REQUIRED_FIELDS).raise_for_errors()
    suuuupclInfo = MaterialSystemInfo()
    suuuupclInfo.charge     = supclInfo.charge
    suuuupclInfo.charge_pos = supclInfo.charge_pos

    suuuupclInfo.atomconfig = patch_atom_v2(  # >log
        supclInfo.atomconfig, bulkInfo.atomconfig, supcl_size, target_size)
    suuuupclInfo.vr = patch_vr(  # >log
        supclInfo.vr, bulkInfo.vr, supcl_size, target_size)
    return suuuupclInfo


# 这个代码只适合偶数 target_size

def patch_atom(supclAtom: AtomConfig, bulkAtom: AtomConfig,
               supcl_size, target_size) -> AtomConfig:
    """
    resolution is bulk, supcl_size should be even numbers
    """
    log(f"{patch_atom.__name__}")
    supclAtom.revise_atomsposition()
    bulkAtom.revise_atomsposition()

    indx = supclAtom.positions >= 0.5
    supclAtom.positions[indx] -= 1.0
    # code prototype from AtomConfig.__mul__
    nrepeat         = np.prod(target_size)
    natoms          = bulkAtom.natoms * nrepeat
    lattice         = bulkAtom.lattice * target_size
    atoms_itype     = np.tile(bulkAtom.itypes, nrepeat)
    atoms_position  = np.tile(bulkAtom.positions, (nrepeat, 1))
    atoms_move      = np.tile(bulkAtom.moves, (nrepeat, 1))

    # supcl atoms
    atoms_itype[0:supclAtom.natoms]     = supclAtom.itypes
    atoms_position[0:supclAtom.natoms]  = supclAtom.positions * supcl_size
    atoms_move[0:supclAtom.natoms]      = supclAtom.moves

    # bulk atoms
    cnt = 0
    for i in range(-target_size[0]//2, target_size[0]//2): # 5 -3...1; 4 -2...1
        for j in range(-target_size[1]//2, target_size[1]//2):
            for k in range(-target_size[2], target_size[2]//2):
                if not ((-supcl_size[0]//2 <= i < supcl_size[0]//2) or
                        (-supcl_size[1]//2 <= j < supcl_size[1]//2) or
                        (-supcl_size[2]//2 <= k < supcl_size[2]//2)):
                    indx = bulkAtom.natoms * (cnt + np.prod(supcl_size))
                    atoms_position[indx:indx+bulkAtom.natoms] += np.array([i, j, k])
                    atoms_itype[indx:indx+bulkAtom.natoms] = bulkAtom.itypes
                    atoms_move[indx:indx+bulkAtom.natoms] = bulkAtom.moves
                    cnt += 1

    atoms_position /= target_size
    rtac = AtomConfig(natoms=natoms, lattice=lattice,
                      itypes=atoms_itype, positions=atoms_position,
                      moves=atoms_move)
    rtac.revise_atomsposition()
    return rtac


# TODO 如果奇数超胞, 奇数扩胞, 会出现什么问题吗
# TODO 如果是奇数网格怎么办
@timing()
def patch_atom_v2(supclAtom: AtomConfig, bulkAtom: AtomConfig,
                  supcl_size, target_size) -> AtomConfig:
    """
    another implement of patch_atom, break the bulkcell-resolution.
    but hidden issues may arise.
    """
    log(f"{patch_atom_v2.__name__}")
    supclAtom.revise_atomsposition()
    bulkAtom.revise_atomsposition()

    assert len(supcl_size) == 3
    assert len(target_size) == 3

    supclAtom.positions[supclAtom.positions >= 0.5] -= 1.0
    # code prototype from AtomConfig.__mul__
    nrepeat         = np.prod(target_size)
    natoms          = bulkAtom.natoms * nrepeat
    lattice         = bulkAtom.lattice * target_size
    atoms_itype     = np.tile(bulkAtom.itypes, nrepeat)
    atoms_position  = np.tile(bulkAtom.positions, (nrepeat, 1))
    atoms_move      = np.tile(bulkAtom.moves, (nrepeat, 1))

    # supcl atoms
    atoms_itype[0:supclAtom.natoms]     = supclAtom.itypes
    atoms_position[0:supclAtom.natoms]  = supclAtom.positions * supcl_size
    atoms_move[0:supclAtom.natoms]      = supclAtom.moves

    # bulk atoms
    cnt, offset = 0, supclAtom.natoms

    def in_imbox(pos, planl, planr):
        if all(planl <= pos) and all(pos < planr):
            return True
        else:
            return False

    target_size = np.array(target_size)
    planl, planr = -supcl_size/2, +supcl_size/2
    for i in range(target_size[0]):
        for j in range(target_size[1]):
            for k in range(target_size[2]):
                ijk = np.array((i, j, k))
                for ia in range(bulkAtom.natoms):
                    pos = bulkAtom.positions[ia] + ijk
                    pos[pos >= target_size / 2] -= \
                        target_size[pos >= target_size / 2]
                    if not in_imbox(pos, planl, planr):
                        indx = cnt + offset
                        atoms_position[indx]  = pos
                        atoms_itype[indx]     = bulkAtom.itypes[ia]
                        atoms_move[indx]      = bulkAtom.moves[ia]
                        cnt += 1
    assert cnt + offset == natoms, \
        "some thing wrong when setting atoms at patch_atom_v2"

    atoms_position /= target_size
    rtac = AtomConfig(natoms=natoms, lattice=lattice,
                      itypes=atoms_itype, positions=atoms_position,
                      moves=atoms_move)
    rtac.revise_atomsposition()
    return rtac


@timing()
def patch_vr(supclVR: VR, bulkVR: VR, supcl_size, target_size) -> VR:
    log(f"{patch_vr.__name__}")
    # The alloy example uses fractional size and intentionally overwrites a
    # smaller physical cell; the grid-spacing requirement applies to ordinary
    # integer supercells only.
    if np.all(np.asarray(supcl_size) == np.round(supcl_size)):
        expected_n123 = bulkVR.n123 * supcl_size
        if not np.array_equal(supclVR.n123, expected_n123):
            raise ValueError(
                f"supercell VR N123={supclVR.n123} must equal "
                f"bulk N123 * supercell size={expected_n123}; "
                "call inspect_ingredient first"
            )
    suuuupclVR = bulkVR * target_size

    # @jit(nopython=True)
    @timing()
    @guvectorize('void(float64[:,:,:], float64[:,:,:])',
                 '(n1,n2,n3),(m1,m2,m3)')
    def _overwrite_supcl_mesh(supcl_mesh, suuuupcl_mesh):
        n1, n2, n3 = supcl_mesh.shape
        for i in range(-n1//2+1, n1//2+1):
            for j in range(-n2//2+1, n2//2+1):
                for k in range(-n3//2+1, n3//2+1):
                    suuuupcl_mesh[i, j, k] = supcl_mesh[i, j, k]
    _overwrite_supcl_mesh(supclVR.mesh, suuuupclVR.mesh)
    return suuuupclVR
