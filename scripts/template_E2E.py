from potpatch import MaterialSystemInfo, inspect_ingredient, gen_charge_correct, edge_match_correct, patch


def potpatch():
    target_size        = [8, 8, 8]
    charge             = 1
    epsilon            = 12.34

    supcl_size         = [4, 4, 4]
    frozen_range       = 1.0
    
    bulk_atomconfig    = "./1.lda/1.bulk/2.scf/atom.config"
    bulk_vr            = "./1.lda/1.bulk/2.scf/OUT.VR"
    bulkInfo    = MaterialSystemInfo()
    bulkInfo.atomconfig.read(filename=bulk_atomconfig)
    bulkInfo.vr.read(filename=bulk_vr, vr_fmt="Escan")

    supcl_atomconfig    = "./1.lda/2.supercell/2.scf/atom.config"
    supcl_vr            = "./1.lda/2.supercell/2.scf/OUT.VR"
    supclInfo   = MaterialSystemInfo(charge=charge, epsilon=epsilon)
    supclInfo.atomconfig.read(filename=supcl_atomconfig)
    supclInfo.vr.read(filename=supcl_vr, vr_fmt="Escan")

    output_atomconfig  = "./3.escan/2.supercell_08a/atom.config2"
    output_vr          = "./3.escan/2.supercell_08a/IN.VR2"
    
    supcl_size = inspect_ingredient(supclInfo, bulkInfo, size_confirm=supcl_size, frozen_confirm=frozen_range)

    minus_V_periodic, plus_V_single = gen_charge_correct(supclInfo)
    minus_V_periodic(supclInfo)
    edge_match_correct(supclInfo, bulkInfo)
    suuuupclInfo = patch(supclInfo, bulkInfo, supcl_size, target_size)
    plus_V_single(suuuupclInfo)

    suuuupclInfo.atomconfig.write(filename=output_atomconfig, fmt="Escan")
    suuuupclInfo.vr.write(filename=output_vr, fmt="Escan", nnodes=None)


if __name__ == "__main__":
    potpatch()
