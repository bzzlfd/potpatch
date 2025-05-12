from potpatch import VR


def main():
    bulk_vr            = ""
    output_vr          = ""

    bulkVR = VR()
    bulkVR.read(filename=bulk_vr, vr_fmt="PWmat")
    bulkVR.write(filename=output_vr, fmt="Escan", nnodes=None)
    

if __name__ == "__main__":
    main()
