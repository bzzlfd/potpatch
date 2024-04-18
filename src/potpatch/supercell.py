import numpy as np

from potpatch.objects import AtomConfig, Lattice

def make_supercell(bulkAtom:AtomConfig, m123) -> AtomConfig :
    "``m123`` should be a tuple contains 3 integer numbers, e.g. ``(4,4,4)``"
    supclAtom = bulkAtom * m123
    supclAtom.revise_atomsposition()
    return supclAtom

def modify_supercell(supclAtom:AtomConfig, frozen_edge_width=1.0) -> AtomConfig:
    """
    for supercell RELAX job: freeze atoms whose positions are closed to the pathcing edge(interface between supercell and bulk)
    frozen_edge_width: in angstrom
    """
    rtac = supclAtom * (1,1,1) # make a copy
    
    for (i, pos) in enumerate(rtac.positions):
        if closed_to_edge(rtac.lattice, pos, frozen_edge_width):
            rtac.moves[i] = np.array([0,0,0])
    
    rtac.revise_atomsposition()
    rtac.sort_atomposition()
    return rtac

def closed_to_edge(lattice:Lattice, position, frozen_edge_width):
    position = np.array(position)
    position %= 1
    indx = position >= 0.5
    position[indx] -= 1.0

    AL = lattice.in_unit("angstrom")
    iAL = (np.linalg.inv(AL))
    iuAL = (iAL / np.linalg.norm(iAL, axis=0)).T
    L = np.array([np.abs(np.dot(AL[i],iuAL[i])) for i in range(3)])
    # print(L*position)
    if any(position*L <= -0.5*L+frozen_edge_width) or any(position*L >= 0.5*L-frozen_edge_width): 
        return True
    else:
        return False
    

