from potpatch import VR, Atom


def main():
    bulk_ac            = ""
    output_ac          = ""
    bulk_vr            = ""
    output_vr          = ""

    bulkAC = Atom()
    bulkAC.read(filename=bulk_ac, fmt="PWmat")
    bulkAC.write(filename=output_ac, fmt="Escan")

    bulkVR = VR()
    bulkVR.read(filename=bulk_vr, fmt="PWmat")
    bulkVR.write(filename=output_vr, fmt="Escan", nnodes=None)
    

if __name__ == "__main__":
    main()
