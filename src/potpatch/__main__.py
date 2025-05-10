from os.path import join, dirname, basename, isabs
from os import getcwd
from pathlib import Path
from textwrap import indent, dedent
from time import perf_counter

import numpy as np
from numpy import prod, array, diag, zeros

from potpatch.constant import BOHR, HA, EPSILON0
from potpatch.objects import (Lattice, 
                              VR, AtomConfig, VATOM, EIGEN, 
                              MaterialSystemInfo)
from potpatch.correction import gen_charge_correct, edge_match_correct
from potpatch.patch import (patch, patch_vr, patch_atom, patch_atom_v2,
                            inspect_ingredient)

from potpatch.supercell import (make_supercell, modify_supercell, 
                                which_lattice_is_bulk)
from potpatch.shift import shift_oneAtomConfig, shift_twoAtomConfig
from potpatch.atompos_coin import check_atompos_consistency
from potpatch.diff_vatom import diff_vatom, write_diffvatom
from potpatch.parse import cli_arg_parse, file_input_parse
from potpatch.utils import revise_epsilon


def main():
    PROG, args = cli_arg_parse()
    PROG, args = file_input_parse(PROG, args)
    if PROG == "potpatch":
        potpatch(args)
    elif PROG == "mksupcl":
        mksupcl(args)
    elif PROG == "shift":
        shift(args)
    elif PROG == "check_atompos":
        check_atompos(args)
    else:
        assert False, f"Invalid PROG {PROG}"


def potpatch(args):
    target_size        = args.output.size
    # [
    correction = args.correction
    charge     = correction.charge
    epsilon    = correction.epsilon
    epsilon    = revise_epsilon(epsilon)

    supcl_size         = args.supcl.size
    frozen_range       = args.supcl.frozen_range
    
    basedir = args.bulk.basedir if isabs(args.bulk.basedir) \
        else join(args.inputfile_dir, args.bulk.basedir)
    bulk_atomconfig    = join(basedir, args.bulk.atomconfig)
    bulk_vr            = join(basedir, args.bulk.vr        )
    bulkInfo    = MaterialSystemInfo(
        atoms_filename=bulk_atomconfig,  vr_filename=bulk_vr, 
        charge=0, epsilon=epsilon)
    # [
    basedir = args.supcl.basedir if isabs(args.supcl.basedir) \
        else join(args.inputfile_dir, args.supcl.basedir)
    supcl_atomconfig   = join(basedir, args.supcl.atomconfig)
    supcl_vr           = join(basedir, args.supcl.vr        )
    supclInfo   = MaterialSystemInfo(
        atoms_filename=supcl_atomconfig, vr_filename=supcl_vr, 
        charge=charge, epsilon=epsilon)
    # [
    output_atomconfig  = args.output.atomconfig \
        if args.output.atomconfig is not None \
        else f"atom.config_{bulkInfo.atomconfig.natoms*prod(target_size)}"
    output_vr          = args.output.vr         \
        if args.output.vr is not None \
        else f"IN.VR_{bulkInfo.atomconfig.natoms*prod(target_size)}"
    basedir = args.output.basedir if isabs(args.output.basedir) \
        else join(args.inputfile_dir, args.output.basedir)
    output_atomconfig  = join(basedir, output_atomconfig)
    output_vr          = join(basedir, output_vr        )
    
    supcl_size = inspect_ingredient(
        supclInfo, bulkInfo, 
        size_confirm=supcl_size, frozen_confirm=frozen_range)
    if args.onlyinspect:
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
            r"corretion:",
            f"    charge: {supclInfo.charge}",
            # f"    charge_pos: {supclInfo.charge_pos}",
            f"    epsilon:",
            f"{indent(supclInfo.epsilon.__str__(), ' '*(4+9))}",
            f"    plus_V: {correction.plus_V}"
            r"output:",
            f"    target_size: {target_size}",
            f"    VR({output_vr})",
            f"    atom.config({output_atomconfig})",
            sep="\n"
        )
        return
    
    minus_V_periodic, plus_V_sop = gen_charge_correct(supclInfo, correction)
    minus_V_periodic(supclInfo)
    edge_match_correct(supclInfo, bulkInfo)
    if (debg := True):
        supclInfo.vr.write_vr(filename="supcl.VR.debug")
        
    suuuupclInfo = patch(supclInfo, bulkInfo, supcl_size, target_size)
    suuuupclInfo.charge, suuuupclInfo.epsilon = charge, epsilon
    plus_V_sop(suuuupclInfo)

    suuuupclInfo.atomconfig.write_atoms(filename=output_atomconfig)
    suuuupclInfo.vr.write_vr(filename=output_vr)

    if args.check.diff_vatom is not None:
        outfile = args.check.diff_vatom.output \
            if args.check.diff_vatom.output is not None \
            else join(args.inputfile_dir, "OUT.diff_vatom")
        sigma  = args.check.diff_vatom.sigma  \
            if args.check.diff_vatom.sigma  is not None else 1 * BOHR
        sigma  = float(sigma) / BOHR

        plus_V_sop(supclInfo)
        r, ξ, dv_bulk, dv_supcl, ac_bulk, ac_supcl, order = \
            diff_vatom(bulkInfo, supclInfo, sigma)
        
        # plus_V_single(supclInfo,    inverse=True)
        # minus_V_periodic(supclInfo, inverse=True)
        supclInfo   = MaterialSystemInfo(
            atoms_filename=supcl_atomconfig, vr_filename=supcl_vr, 
            charge=charge, epsilon=epsilon)
        _, _, _, dv_supcl_0, _, _, _ = \
            diff_vatom(bulkInfo, supclInfo, sigma)

        write_diffvatom(outfile, ac_bulk, ac_supcl, 
                        order, epsilon, charge, r, ξ,
                        dv_bulk, dv_supcl_0, dv_supcl)


def mksupcl(args):
    input_          = args.input_
    output          = args.output
    size            = args.size 
    frozen_range    = args.frozen_range 
    outsider        = args.outsider

    bulk_ac = AtomConfig(filename=input_)
    if output is None:
        output = f"atom.config_mksupcl_{bulk_ac.natoms*prod(size)}"
    
    supcl_ac = make_supercell(bulk_ac, size)
    supcl_ac = modify_supercell(supcl_ac, frozen_range, mv2outsider=outsider)  # `supcl_ac` has new attribute _natoms_frozen
    
    comment = f"/{supcl_ac._natoms_frozen} atoms frozen/ " + \
        f"potpatch mksupcl -i {input_} -o {output} -r {frozen_range}" + \
        f" -s {' '.join(str(f) for f in size)} "
    comment += f"-d {outsider} " if outsider is not None else ""
    comment += f"# pwd={getcwd()}"

    supcl_ac.write_atoms(output, comment=comment)


def shift(args):
    bulk    = args.bulk
    supcl   = args.supcl
    shift   = [float(s) for s in args.shift]
    comment = f"potpatch shift -B {bulk} -S {supcl} " + \
        f"-s {' '.join(str(f) for f in shift)} # pwd={getcwd()}"

    assert args.count == 1 or args.count == 2
    if args.count == 1:
        if bulk is not None:
            ac    = AtomConfig(filename=bulk)
        if supcl is not None:
            ac    = AtomConfig(filename=supcl)
        shift_oneAtomConfig(ac, shift)
        ac.write_atoms(ac.filename + "_shift", comment=comment)
    elif args.count == 2:
        bulk_ac   = AtomConfig(filename=bulk)
        supcl_ac  = AtomConfig(filename=supcl)
        shift_twoAtomConfig(bulk_ac, supcl_ac, shift)
        bulk_ac.write_atoms(bulk_ac.filename + "_shift", comment=comment)
        supcl_ac.write_atoms(supcl_ac.filename + "_shift", comment=comment)


def check_atompos(args):
    ac_1 = AtomConfig(filename=args.ac_1)
    ac_2 = AtomConfig(filename=args.ac_2)
    whichbulk = which_lattice_is_bulk(ac_1.lattice, ac_2.lattice)
    if whichbulk == 1:
        bulk, supcl = ac_1, ac_2
    elif whichbulk == 2:
        bulk, supcl = ac_2, ac_1
    nwarns, max_dist = check_atompos_consistency(bulk, supcl,
                                                 tol=args.tol)
    if nwarns == 0:
        print(f"No atom position inconsistency found.")
    else:
        print(
            f"{nwarns} atom position inconsistency found "
            f"(max_dist={max_dist} angstrom).")


if __name__ == "__main__":
    main()
