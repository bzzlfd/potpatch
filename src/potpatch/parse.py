import argparse
import tomllib
from os.path import dirname, abspath
from os import getcwd

from potpatch.utils import NameTuple
from potpatch.version import __version__


# TODO replace `if-elif-else`` into `match-case`` when everyone in the world use Python >= 3.10
def cli_arg_parse():
    parser_potpatch = argparse.ArgumentParser()
    parser_potpatch.add_argument(
         "-i", "--input", 
         help='specify the input file name. default is "potpatch.input"', 
         default="potpatch.input")
    parser_potpatch.add_argument(
         "-I", "--only-inspect", 
         help='just inspecting ingredients without performing calculations.', 
         action="store_true", 
         dest="onlyinspect", 
         default=False)
    parser_potpatch.add_argument(
         "-V", "--version", 
         action="store_true")
    
    subparsers = parser_potpatch.add_subparsers(
         help='sub-command help', 
         dest="func")

    parser_mksupcl = subparsers.add_parser(
         "mksupcl", 
         help="a script to make supercell")
    parser_mksupcl.add_argument(
         "-i", "--input", 
         help='specify the input file name. default is "potpatch.input"', 
         default="potpatch.input")
    # TODO mksupcl is a script. should support passing param in cli, params outside will overwrite the param in input file. 

    args = parser_potpatch.parse_args()
    if args.func is None:
        PROG = "potpatch"
        if args.version:
            print(f"potpatch version {__version__}")
            exit()
        else:
            inputfile = args.input
            onlyinspect = args.onlyinspect
            ret_nt = NameTuple(inputfile=inputfile, 
                               onlyinspect=onlyinspect)
    elif args.func == "mksupcl":
        PROG = "mksupcl"
        inputfile = args.input
        ret_nt = NameTuple(inputfile=inputfile)
    else:
        raise ValueError(f'Invalid command {args.func}')

    return PROG, ret_nt


def file_input_parse(PROG, args):
    inputfile = args.inputfile
    with open(inputfile, "rb") as f:
        inputfile_path = dirname(abspath(f.name))

        settings = tomllib.load(f)

        if PROG == "potpatch":
            onlyinspect = args.onlyinspect
            sub_settings = settings["potpatch"]

            bulk = NameTuple(
                atomconfig = sub_settings["bulk"]["atomconfig"], 
                vr         = sub_settings["bulk"]["vr"])
            supcl = NameTuple(
                atomconfig   = sub_settings["supercell"]["atomconfig"], 
                vr           = sub_settings["supercell"]["vr"], 
                charge       = sub_settings["supercell"]["charge"], 
                charge_pos   = sub_settings["supercell"].get("charge_pos", None), 
                epsilon      = sub_settings["supercell"]["epsilon"], 
                size         = sub_settings["supercell"].get("size", None),
                frozen_range = sub_settings["supercell"].get("frozen_range", None))
            
            output = NameTuple(
                size       = sub_settings["output"]["size"], 
                atomconfig = sub_settings["output"].get("atomconfig", None), 
                vr         = sub_settings["output"].get("vr", None))
            
            ret_nt = NameTuple(inputfile_path = inputfile_path, 
                               bulk           = bulk, 
                               supcl          = supcl, 
                               output         = output, 
                               onlyinspect    = onlyinspect, 
                               )
        
        elif PROG == "mksupcl":
                sub_settings = settings["mksupcl"]
                bulk_atomconfig = sub_settings["bulk_atomconfig"]
                supcl_size      = sub_settings["supercell_size"]
                frozen_range    = sub_settings.get("frozen_range", None)
                output          = sub_settings.get("output", None)

                ret_nt = NameTuple(inputfile_path  = inputfile_path, 
                                   bulk_atomconfig = bulk_atomconfig, 
                                   supcl_size      = supcl_size, 
                                   frozen_range    = frozen_range,
                                   output          = output
                                   )
        else:
            assert False, f"Invalid PROG {PROG}"
    
    return PROG, ret_nt
