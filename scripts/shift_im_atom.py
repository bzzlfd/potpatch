"""
在 im_pos 推出前的临时办法
"""

from potpatch import AtomConfig
import numpy as np


def main():
    bulk_ac     = ""
    output_ac   = ""
    im_pos      = np.array([0.0, 0.0, 0.0])

    bulkAC = AtomConfig()
    bulkAC.read(filename=bulk_ac, fmt="PWmat")
    bulkAC.positions -= im_pos
    bulkAC.revise_atomsposition()
    bulkAC.write(filename=output_ac, fmt="PWmat")
    

if __name__ == "__main__":
    main()
