import argparse
import tomllib
from os.path import dirname, abspath
from os import getcwd

from potpatch.utils import NameTuple
from potpatch.version import __version__


# TODO replace `if-elif-else`` into `match-case`` when everyone in the world use Python >= 3.10
def cli_arg_parse():
    parser_potpatch = argparse.ArgumentParser()
    # cli potpatch
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
    
    # cli sub-commands
    subparsers = parser_potpatch.add_subparsers(
         help='sub-command help', 
         dest="func")

    # cli mksupcl 
    parser_mksupcl = subparsers.add_parser(
         "mksupcl", 
         help="a script to make supercell and freeze atoms on the patching edge")
    parser_mksupcl.add_argument(
         "-i", "--input", 
         help='specify the input bulk file name.', 
         dest="input_",
         required=True)
    parser_mksupcl.add_argument(
         "-o", "--output", 
         help='specify the output supercell file name. '
              'default is `atom.config_mksupcl_{natoms}`', 
         default=None)
    parser_mksupcl.add_argument(
         "-r", "--range", 
         help='specify the frozen range (angstrom) on the edge. '
              'default is `1.0`', 
         dest="frozen_range",
         type=float,
         default=1.0)
    parser_mksupcl.add_argument(
         "-s", "--size", 
         help='specify the size of supercell (e.g. `--size 4 4 4`). ' 
              'default is `[1, 1, 1]`', 
         type=int,
         default=[1, 1, 1],
         nargs=3)

    # cli shift atom.config 
    parser_shift = subparsers.add_parser(
         "shift", 
         help="a script to shift one atom.config / "
              "both bulk&supercell two atom.config(s)")
    parser_shift.add_argument(
         "-B", "--bulk", 
         help='specify the bulk file name. '
              'output is the same filename suffixed by "_shift"', 
         default=None)
    parser_shift.add_argument(
         "-S", "--supcl", 
         help='specify the supercell file name. '
              'output is the same filename suffixed by "_shift"', 
         default=None)
    parser_shift.add_argument(
         "-s", "--shift",
         help="specify the fractional position `shift` "
              "(e.g. `--shift 0.25 0.25 0.75`.) "
              "shift `[shift, shift, shift]` to [0, 0, 0]. "
              "if both -B and -S are specified, it will be regarded as "
              "supercell's `shift`, bulk's `shift` can be inferred from "
              "the existing information", 
         type=float,
         nargs=3,
         required=True)

    # parse_args
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
        input_       = args.input_  # `input` is a function, and `input_` is aligned to `output``
        output       = args.output
        frozen_range = args.frozen_range
        size         = args.size
        ret_nt = NameTuple(input_=input_, output=output,
                           frozen_range=frozen_range, size=size)
    elif args.func == "shift":
        PROG = "shift"
        bulk  = args.bulk
        supcl = args.supcl
        shift = args.shift
        if bulk is None and supcl is None: 
            print("one of {--bulk, --supcl} must be specified")
            exit(255)
        elif bulk is None or supcl is None: 
            count = 1
        else:
            count = 2
        ret_nt = NameTuple(count=count, bulk=bulk, supcl=supcl, shift=shift)
    else:
        raise ValueError(f'Invalid command {args.func}')

    return PROG, ret_nt


def file_input_parse(PROG, args):
    if hasattr(args, "inputfile"):  # case `potpatch`
        inputfile = args.inputfile
        with open(inputfile, "rb") as f:
            inputfile_dir = dirname(abspath(f.name))
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
        
        ret_nt = NameTuple(inputfile_dir  = inputfile_dir, 
                           bulk           = bulk, 
                           supcl          = supcl, 
                           output         = output, 
                           onlyinspect    = onlyinspect, 
                           )
    
    # scripts
    elif PROG == "mksupcl":
        ret_nt = args
    elif PROG == "shift":
        ret_nt = args
    else:
        assert False, f"Invalid PROG {PROG}"
    
    return PROG, ret_nt
