# 边写注释边执行
from os.path import join
from textwrap import indent, dedent
from time import perf_counter

from numpy import prod

from potpatch.constant import BOHR, HA, EPSILON0
from potpatch.objects import (Lattice, 
                              VR, AtomConfig, VATOM, EIGEN, 
                              MaterialSystemInfo)
from potpatch.correction import gen_charge_correct, edge_match_correct
from potpatch.patch import (patch, patch_vr, patch_atom, patch_atom_v2,
                            inspect_ingredient)

from potpatch.supercell import make_supercell, modify_supercell
from potpatch.parse import cli_arg_parse, file_input_parse

def main():
    PROG, args = cli_arg_parse()
    PROG, args = file_input_parse(PROG, args)
    if PROG == "potpatch":
        potpatch(args)
    elif PROG == "mksupcl":
        mksupcl(args)
    else:
        assert False, f"Invalid PROG {PROG}"


def potpatch(args):
    target_size        = args.output.size
    charge             = args.supcl.charge
    epsilon            = args.supcl.epsilon

    supcl_size         = args.supcl.size
    frozen_range       = args.supcl.frozen_range
    
    bulk_atomconfig    = join(args.inputfile_path, args.bulk.atomconfig  )
    bulk_vr            = join(args.inputfile_path, args.bulk.vr          )
    bulkInfo    = MaterialSystemInfo(atoms_filename=bulk_atomconfig,  vr_filename=bulk_vr)
    supcl_atomconfig   = join(args.inputfile_path, args.supcl.atomconfig )
    supcl_vr           = join(args.inputfile_path, args.supcl.vr         )
    supclInfo   = MaterialSystemInfo(atoms_filename=supcl_atomconfig, vr_filename=supcl_vr, charge=charge, epsilon=epsilon)
    output_atomconfig  = args.output.atomconfig if args.output.atomconfig is not None else f"atom.config_{bulkInfo.atomconfig.natoms*prod(target_size)}"
    output_vr          = args.output.vr         if args.output.vr         is not None else f"IN.VR_{bulkInfo.atomconfig.natoms*prod(target_size)}"
    output_atomconfig  = join(args.inputfile_path, output_atomconfig)
    output_vr          = join(args.inputfile_path, output_vr        )
    
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
    epsilon: {supclInfo.epsilon}
"""
        )
        return
    
    minus_V_periodic, plus_V_single = gen_charge_correct(supclInfo)
    minus_V_periodic(supclInfo)
    edge_match_correct(supclInfo, bulkInfo)
    suuuupclInfo = patch(supclInfo, bulkInfo, supcl_size, target_size)
    plus_V_single(suuuupclInfo)

    suuuupclInfo.atomconfig.write_atoms(filename=output_atomconfig)
    suuuupclInfo.vr.write_vr(filename=output_vr)


def mksupcl(args):
    supcl_size      = args.supcl_size 
    frozen_range    = args.frozen_range if args.frozen_range is not None else 1.0

    bulk_atomconfig    = join(args.inputfile_path, args.bulk_atomconfig)
    bulkAtom        = AtomConfig(filename=bulk_atomconfig)
    output          = args.output       if args.output       is not None else f"atom.config_{bulkAtom.natoms*prod(supcl_size)}"
    output             = join(args.inputfile_path, output              )
    
    supclAtom = make_supercell(bulkAtom, supcl_size)
    supclAtom = modify_supercell(supclAtom, frozen_range)
    supclAtom.write_atoms(output)


if __name__ == "__main__":
    main()
