from potpatch.shift import shift_vr
from potpatch import VR

file = ""
shift = (0.5, 0.5, 0.5)

vr = VR(filename=file)
shift_vr(vr, shift)
vr.write(file + "_shift")
