from os.path import join, dirname, basename
from os import getcwd
from pathlib import Path
from textwrap import indent, dedent
from time import perf_counter

import numpy as np
from numpy import prod, array, zeros

from potpatch.constant import BOHR, HA, EPSILON0
from potpatch.objects import (Lattice, 
                              VR, AtomConfig, VATOM, EIGEN, 
                              MaterialSystemInfo)
from potpatch.correction import gen_charge_correct, edge_match_correct
from potpatch.patch import (patch, patch_vr, patch_atom, patch_atom_v2,
                            inspect_ingredient)

from potpatch.supercell import make_supercell, modify_supercell
from potpatch.shift import shift_oneAtomConfig, shift_twoAtomConfig
from potpatch.parse import cli_arg_parse, file_input_parse


def main():
    PROG, args = cli_arg_parse()
    PROG, args = file_input_parse(PROG, args)
    if PROG == "potpatch":
        potpatch(args)
    elif PROG == "mksupcl":
        mksupcl(args)
    elif PROG == "shift":
        shift(args)
    else:
        assert False, f"Invalid PROG {PROG}"


def potpatch(args):
    target_size        = args.output.size
    charge             = args.supcl.charge
    charge_pos         = np.array(args.supcl.charge_pos) if args.supcl.charge_pos is not None else np.array([0.,0.,0.])
    epsilon            = args.supcl.epsilon

    supcl_size         = args.supcl.size
    frozen_range       = args.supcl.frozen_range
    
    bulk_atomconfig    = join(args.inputfile_dir, args.bulk.atomconfig )
    bulk_vr            = join(args.inputfile_dir, args.bulk.vr         )
    bulkInfo    = MaterialSystemInfo(atoms_filename=bulk_atomconfig,  vr_filename=bulk_vr)
    supcl_atomconfig   = join(args.inputfile_dir, args.supcl.atomconfig)
    supcl_vr           = join(args.inputfile_dir, args.supcl.vr        )
    supclInfo   = MaterialSystemInfo(atoms_filename=supcl_atomconfig, vr_filename=supcl_vr, charge=charge, charge_pos=charge_pos, epsilon=epsilon)
    output_atomconfig  = args.output.atomconfig if args.output.atomconfig is not None else f"atom.config_{bulkInfo.atomconfig.natoms*prod(target_size)}"
    output_vr          = args.output.vr         if args.output.vr         is not None else f"IN.VR_{bulkInfo.atomconfig.natoms*prod(target_size)}"
    output_atomconfig  = join(args.inputfile_dir, output_atomconfig)
    output_vr          = join(args.inputfile_dir, output_vr        )
    
    supcl_size = inspect_ingredient(supclInfo, bulkInfo, size_confirm=supcl_size, frozen_confirm=frozen_range)
    if args.onlyinspect:
        space4 = " "*4
        print(
        # TODO 这么写太丑了
f"""
bulk:
    lattice: (in angstrom)
{indent(f"{bulkInfo.lattice.in_unit('angstrom').__str__()}", space4*2)}
    VR({bulkInfo.vr.filename}):
{indent(f"n123: {bulkInfo.vr.n123}", space4*2)}
{indent(f"nnodes: {bulkInfo.vr.nnodes}", space4*2)}
    atom.config({bulkInfo.atomconfig.filename}):
supercell:
    size: {supcl_size}
    lattice: (in angstrom)
{indent(f"{supclInfo.lattice.in_unit('angstrom').__str__()}", space4*2)}
    VR({supclInfo.vr.filename}):
{indent(f"n123: {supclInfo.vr.n123}", space4*2)}
{indent(f"nnodes: {supclInfo.vr.nnodes}", space4*2)}
    atom.config({supclInfo.atomconfig.filename}):
    charge: {supclInfo.charge}
    charge_pos: {supclInfo.charge_pos}
    epsilon: {supclInfo.epsilon}
"""
        )
        return
    
    minus_V_periodic, plus_V_single = gen_charge_correct(supclInfo)
    minus_V_periodic(supclInfo)
    edge_match_correct(supclInfo, bulkInfo)
    if (debg := False):
        supclInfo.vr.write_vr(filename="supcl.VR.debug")
    suuuupclInfo = patch(supclInfo, bulkInfo, supcl_size, target_size)
    plus_V_single(suuuupclInfo)

    suuuupclInfo.atomconfig.write_atoms(filename=output_atomconfig)
    suuuupclInfo.vr.write_vr(filename=output_vr)


def mksupcl(args):
    input_          = args.input_
    output          = args.output
    size            = args.size 
    frozen_range    = args.frozen_range 
    outsider        = args.outsider

    bulk_ac = AtomConfig(filename=input_)
    if output is None:
        output = f"atom.config_mksupcl_{bulk_ac.natoms*prod(size)}"
    comment = f"potpatch mksupcl -i {input_} -o {output} " + \
        f"-s {size} -r {frozen_range} # pwd={getcwd()}"
    
    supcl_ac = make_supercell(bulk_ac, size)
    supcl_ac = modify_supercell(supcl_ac, frozen_range, mv2outsider=outsider)
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


if __name__ == "__main__":
    main()
