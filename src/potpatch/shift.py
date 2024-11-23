import numpy as np
from numpy import abs, exp, pi
from numba import jit, guvectorize

from potpatch.objects import MaterialSystemInfo, AtomConfig, VR, Lattice
from potpatch.supercell import infer_supercell_size
from potpatch.utils import timing, gen_counter
from potpatch.datatype import INTEGER

shift_operation_counter = gen_counter()


def shift_oneAtomConfig(ac: AtomConfig, shift) -> None:
    """
    shift *fractional* point `shift` to `(0, 0, 0)` 
    """
    if shift_operation_counter() == 1: 
        print("shift *fractional* point `shift` to `(0, 0, 0)`")
    assert np.shape(shift) == np.zeros(3).shape

    shift_atomconfig(ac, shift)


def shift_twoAtomConfig(bulk_ac: AtomConfig, 
                        supcl_ac: AtomConfig, supcl_shift) -> None:
    """
    `bulk_shift` refer from `supcl_shift`, `supcl_ac` and `bulk_ac`
    
    shift `supcl_ac` *fractional* point `supcl_shift` to `(0, 0, 0)` 
    shift `bulk_ac`  *fractional* point `bulk_shift`  to `(0, 0, 0)` 
    """
    if shift_operation_counter() == 1: 
        print("... hint: shift supercell *fractional* position"
              " `shift` to `[0,0,0]`")
    assert np.shape(supcl_shift) == np.zeros(3).shape

    # infer the supercell size 
    mag, mag_f = infer_supercell_size(bulk_ac.lattice, supcl_ac.lattice)
    print("I infer that the size of supercell is "
          f"{','.join([str(i) for i in mag_f])}")
    
    shift_oneAtomConfig(supcl_ac, supcl_shift)
    shift_oneAtomConfig(bulk_ac, supcl_shift * mag_f)


def shift_materialsystem(ms: MaterialSystemInfo, shift) -> None:
    """
    shift *fractional* point `shift` to `(0, 0, 0)` 
    """
    if shift_operation_counter() == 1: 
        print("shift *fractional* point `shift` to `(0, 0, 0)`")
    assert np.shape(shift) == np.zeros(3).shape

    shift_atomconfig(ms.atomconfig, shift)
    shift_vr(ms.vr, shift)
    if ms.charge_pos:
        ms.charge_pos = ms.charge_pos - np.array(shift)


def shift_atomconfig(ac: AtomConfig, shift, revise_pos=True) -> None:
    """
    shift *fractional* point `shift` to `(0, 0, 0)` 
    """
    if shift_operation_counter() == 1: 
        print("shift *fractional* point `shift` to `(0, 0, 0)`")
    assert np.shape(shift) == np.zeros(3).shape

    shift = np.array(shift)
    ac.positions -= shift
    if revise_pos:
        ac.revise_atomsposition()


def shift_vr(vr: VR, shift) -> None:
    """
    shift $vr(r)$ to $vr(r+r_shift)$, i.e. 
    shift *fractional* point `shift` to `(0, 0, 0)` 
    """
    if shift_operation_counter() == 1: 
        print("shift *fractional* point `shift` to `(0, 0, 0)`")
    assert np.shape(shift) == np.zeros(3).shape

    shift = np.array(shift)
    meshG = np.fft.rfftn(vr.mesh, axes=(0, 1, 2))
    _shift_Rmesh_atG(meshG, shift)
    vr.mesh = np.fft.irfftn(meshG, axes=(0, 1, 2))


@timing()
@guvectorize('void(complex128[:,:,:], float64[:]) ', 
             '(n1,n2,n3),(n)', nopython=True)
def _shift_Rmesh_atG(meshG: np.ndarray, shift):
    "V(r) -> V(r + r_shift) in G space"
    n1, n2, n3 = meshG.shape
    for i in range(0, n1):
        for j in range(0, n2):
            for k in range(0, n3):
                ip = shift[0]*i + shift[1]*j + shift[2]*k
                meshG[i, j, k] *= exp(1j*2*pi*ip)
