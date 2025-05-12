from os.path import join
from textwrap import indent, dedent
from numbers import Rational
import warnings
from pathlib import Path

from numpy import ndarray, array, zeros, diag, prod, abs, round

from potpatch import (AtomConfig, VR, MaterialSystemInfo,
                      gen_charge_correct, check_atompos_consistency, 
                      edge_match_correct, patch, 
                      diff_vatom, write_diffvatom,
                      BOHR, INTEGER)


# a copy from __main__.potpatch
# main difference is the float element type of `supcl_size``
def potpatch():
    # potpatch
    root_dir = join(Path(__file__).parent, "")
    charge             = 1                    # (setting NUM_ELECTRON) - (default NUM_ELECTRON)
    charge_pos         = array([0., 0., 0.])  # deprecated, just set to 0.
    epsilon            = 12.34   # can be a scalar or a 3x3 matrix
    bulk_atomconfig    = ""
    bulk_vr            = ""
    supcl_atomconfig   = ""
    supcl_vr           = ""
    output_atomconfig  = ""
    output_vr          = ""

    # optional, for inspect_ingredient
    frozen_confirm = 1.0  # angstrom
    supcl_size_confirm = [0.2, 0.2, 0.2]

    # optional, for diff_vatom
    check_diff_vatom = True
    diff_vatom_outfile = "OUT.diff_vatom"
    sigma = 1.0  # angstrom

    # ============================================================ 
    target_size        = [1, 1, 1]  # in this approach, alloy `bulk` is larger than `supercell`
    # internal epsilon must be a 3x3 np.ndarray
    if isinstance(epsilon, (list, tuple, ndarray)):
        epsilon = array(epsilon)
    else:
        epsilon = float(epsilon)
        epsilon = diag([epsilon, epsilon, epsilon])

    # 
    bulk_atomconfig    = join(root_dir, bulk_atomconfig )
    bulk_vr            = join(root_dir, bulk_vr         )
    supcl_atomconfig   = join(root_dir, supcl_atomconfig)
    supcl_vr           = join(root_dir, supcl_vr        )
    bulkInfo           = MaterialSystemInfo(atoms_filename=bulk_atomconfig,  vr_filename=bulk_vr)
    supclInfo          = MaterialSystemInfo(atoms_filename=supcl_atomconfig, vr_filename=supcl_vr, charge=charge, charge_pos=charge_pos, epsilon=epsilon)
    output_atomconfig  = join(root_dir, output_atomconfig)
    output_vr          = join(root_dir, output_vr        )
    
    supcl_size = inspect_ingredient(supclInfo, bulkInfo, 
                                    frozen_confirm, supcl_size_confirm)
    if print_ingredient := True:
        print(
            r"bulk:",
            r"    lattice: (in angstrom)",
            indent(f"{bulkInfo.lattice.in_unit('angstrom').__str__()}", " "*8),
            f"    VR({bulkInfo.vr.filename}):",
            f"        n123: {bulkInfo.vr.n123}",
            f"        nnodes: {bulkInfo.vr.nnodes}",
            f"    atom.config({bulkInfo.atomconfig.filename}):",
            r"supercell:",
            f"    size: {supcl_size}",
            r"    lattice: (in angstrom)",
            indent(f"{supclInfo.lattice.in_unit('angstrom').__str__()}", " "*8),
            f"    VR({supclInfo.vr.filename}):",
            f"        n123: {supclInfo.vr.n123}",
            f"        nnodes: {supclInfo.vr.nnodes}",
            f"    atom.config({supclInfo.atomconfig.filename}):",
            f"    charge: {supclInfo.charge}",
            f"    charge_pos: {supclInfo.charge_pos}",
            f"epsilon: {indent(supclInfo.epsilon.__str__(), ' '*(4+9))}",
            sep="\n"
        )
    
    minus_V_periodic, plus_V_single = gen_charge_correct(supclInfo)
    minus_V_periodic(supclInfo)
    edge_match_correct(supclInfo, bulkInfo)
    if (debg := False):
        supclInfo.vr.write_vr(filename="supcl.VR.debug")
    suuuupclInfo = patch(supclInfo, bulkInfo, supcl_size, target_size)
    plus_V_single(suuuupclInfo)

    suuuupclInfo.atomconfig.write(filename=output_atomconfig)
    suuuupclInfo.vr.write(filename=output_vr)

    if check_diff_vatom:
        outfile = diff_vatom_outfile
        sigma  = float(sigma) / BOHR

        plus_V_single(supclInfo)
        r, ξ, dv, ac_bulk, ac_supcl, order = diff_vatom(bulkInfo, supclInfo, sigma)
        write_diffvatom(outfile, ac_bulk, ac_supcl, order, epsilon, r, ξ, dv)


"""
patch is compatible for alloy bulk case:

    patch_atom_v2 is a overwrite verison, and supcl_size is a "float needed"
      just set target_size to [1, 1, 1]

    patch_vr is a overwrite verison, and supcl_size is not needed
      just set target_size to [1, 1, 1]
"""


def inspect_ingredient(supclInfo, bulkInfo, frozen_confirm=None, size_confirm=None):
    # supcl size inference
    supcl_vrsize = supclInfo.vr.n123 / bulkInfo.vr.n123 
    # Lattice and VR.n123: ?same
    lattice_mulmag = bulkInfo.lattice * supcl_vrsize
    if not supclInfo.lattice == lattice_mulmag:
        raise warnings.warn(dedent(f"""
                magnifacation between Lattice and VR_mesh is not equal
                supclInfo.lattice({supclInfo.lattice.in_unit("angstrom")})
                bulkInfo.lattice * mag({lattice_mulmag.in_unit("angstrom")})
                """))
    supcl_size = array(supcl_vrsize, dtype=float)

    if supclInfo.charge is None:
        warnings.warn("information of supercell `charge` is not given")
    if supclInfo.epsilon is None:
        raise ValueError("information of supercell `epsilon` is not given")
    
    if size_confirm is not None:
        print(f"{size_confirm = } is specified ... ", end="")
        if not all(supcl_size == size_confirm):
            warnings.warn(f"supcl_size = {list(supcl_size)} and size_confirm = {list(size_confirm)}")
            print("")
        else:
            print("no error")
    if (frozen_confirm is not None):
        print(f"{frozen_confirm = } is specified ... ", end="")
        nwarn, _ = check_atompos_consistency(bulkInfo.atomconfig, 
                                             supclInfo.atomconfig, 
                                             frozen_range=frozen_confirm, 
                                             )
        if nwarn == 0:
            print("no error")
        elif nwarn > 0:
            print(f"{nwarn} warnings")

    return supcl_size


if __name__ == "__main__":
    potpatch()
