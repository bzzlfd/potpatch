from potpatch import VR, Lattice, Atom

from numpy import zeros, ones

ref = ""
out = ""
n123 = (32, 32, 32)

ac = Atom(filename=ref)
mesh = ones(n123)
gen_uniform_vr = VR(lattice=ac.lattice, mesh=mesh)
gen_uniform_vr.write(out)
