from potpatch import VR, Atom


def main():
    input_vr           = ""
    output_vr          = ""
    input_atom         = ""
    output_atom        = ""
    

    bulkVR = VR()
    bulkVR.read(filename=input_vr, fmt="PWmat")
    bulkVR.write(filename=output_vr, fmt="Escan", nnodes=None)
    
    bulkVR = Atom()
    bulkVR.read(filename=input_atom, fmt="PWmat")
    bulkVR.write(filename=output_atom, fmt="Escan")
    

if __name__ == "__main__":
    main()
