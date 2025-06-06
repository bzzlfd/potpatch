import os
from copy import deepcopy, copy
import warnings
from math import ceil

import numpy as np

from potpatch.utils import (read_fortran_binary_block, 
                            read_fortran_binary_block_varioustype, 
                            write_fortran_binary_block)
from potpatch.constant import BOHR, HA, EPSILON0
from potpatch.datatype import INTEGER, INTEGER_IN, INTEGER_OUT, REAL_8


"""
TODO type hint, typing and docs
TODO assert, raise 
"""


class Lattice():
    """
    unit: "angstrom" | "atomic unit"
    fromwhere: a `list` recording set-lattice chain
    """
    registered_units = ["angstrom", "atomic unit"]
    au_transformer = dict(zip(registered_units, [BOHR, 1]))  # 1/au_transformer (a.u.) = 1 (registered_unit)
    
    def __init__(self, AL: np.ndarray, unit: str = None, 
                 fromwhere: str = ["unkown"]) -> None:
        self.unit_check(unit)
        self.AL         = AL
        self.unit       = unit
        self.fromwhere  = fromwhere if type(fromwhere) is list else [fromwhere]

    @property
    def AL_AU(self):
        return self.in_unit("atomic unit")

    def modify_unit(self, target_unit: str):
        """
        convert `(self.AL, self.unit)` to `target_unit`
        """
        self.unit_check(target_unit)
        self.AL = self.AL \
            / self.au_transformer[self.unit] * self.au_transformer[target_unit]
        self.unit = target_unit

    def in_unit(self, target_unit: str):
        """
        return `self.AL(self.unit)` to `(target_unit)`
        """
        self.unit_check(target_unit)
        return self.AL \
            / self.au_transformer[self.unit] * self.au_transformer[target_unit]
        
    def unit_check(self, unit):
        assert unit is not None, f"""
            AL doesn't know its unit, 
            you need to provide the unit format along with the input for AL.
            {self.registered_units} are the optional units
            """
        assert unit in self.registered_units, f"""
            unit format of AL is not in {self.registered_units},
            """

    def copy(self, appendwhere: str | None = None) -> 'Lattice':
        """
        return a new lattice which `fromwhere` append `appendwhere` if `appendwhere` is not None
        """
        newlattice = deepcopy(self)

        newlattice.fromwhere = copy(self.fromwhere)
        if appendwhere is not None:
            newlattice.fromwhere.append(appendwhere)

        return newlattice

    def __eq__(self, lattice: 'Lattice') -> bool:
        # TODO 判断类型是不是 Lattice
        return np.all(np.abs(self.AL_AU-lattice.AL_AU) < 1e-5) 

    def __mul__(self, mag):
        """
        magification can be not a integer

        self * magnification
        magnification example: ``(4,4,4)``
        return a new Lattice
        """
        for i in list(mag):
            if not i == round(i, 0):
                print(f"WARNING: the magnification({mag}) is not integer")
        assert len(mag) == 3, f"length of magnification must be 3, but not {len(mag)}"
        AL = np.zeros((3, 3))
        for i in range(3):
            AL[i] = self.AL[i] * mag[i]
        fromwhere = copy(self.fromwhere)
        fromwhere.append(f"(    magnification = {mag})")
        return Lattice(AL, self.unit, fromwhere=fromwhere)

    def __rmul__(self, magnification):
        return self.__mul__(magnification)

    def __imul__(self, magnification):
        lattice = self.__mul__(magnification)
        self.AL = lattice.AL
        self.AL_AU = lattice.AL_AU
        return self

    def __truediv__(self, lattice: 'Lattice') -> np.ndarray[int]:
        # TODO 或许写成矩阵形式更好, 还能照顾到并不是
        # TODO 或许这里不应该限定整数
        raise NotImplementedError()
        assert isinstance(lattice, Lattice), \
            "__truediv__ is designed to get supclLattice/bulkLattice magnification. if you want Lattice/4, use Lattice * (1/4)"

        AL1 = self.in_unit("angstrom")
        AL2 = lattice.in_unit("angstrom")

        div3 = np.zeros(3, dtype=INTEGER)
        for i in range(3):
            _div = [1145, 1419, 19810]
            for j in range(3):
                if (AL1[i, j] == 0 and AL2[i, j] != 0) or (AL1[i, j] != 0 and AL2[i, j] == 0):
                    assert False, f"a{i+1}({AL1[i], AL2[i]})Å in the two lattice is not parallel"  # TODO 一个是0向量, 它平行于任何向量, 但是被放在这个情况里了
                elif (AL1[i, j] == 0 and AL2[i, j] == 0):
                    _div[j] = np.nan
                else:
                    _div[j] = AL1[i, j] / AL2[i, j]
            assert not all([np.isnan(i) for i in _div]), f"a{i+1}({AL1[i], AL2[i]})Å in the two lattice is zero vector"
            
            _div = [i for i in _div if not np.isnan(i)]
            for j in range(len(_div)):
                assert abs(_div[0] - _div[j]) < 1e-6, f"a{i+1}({AL1[i], AL2[i]})Å in the two lattice is not parallel"
            if not _div[0] == round(_div[0], 0):
                warnings.warn(f"the magnification of a{i+1}({AL1[i], AL2[i]})Å between the two lattices is not integer ")
            div3[i] = _div[0]
        return div3

    def __floordiv__(self, lattice: 'Lattice'):
        raise NotImplementedError()
        """
        return a Rational number Vector. 
        Julia user: very in river
        """
        assert isinstance(lattice, Lattice), \
            "__truediv__ is designed to get supclLattice/bulkLattice magnification. if you want Lattice/4, use Lattice * (1/4)"
        pass


class VR():
    """
    OUT.VR
    """
    # TODO 没写关于使用这两个变量时的各种检查代码
    registered_fmt = ["PWmat", "Escan"]
    fmt2unit = dict(zip(registered_fmt, ["angstrom", "atomic unit"]))
    
    def __init__(self, filename=None, 
                 lattice: Lattice = None, mesh=None, 
                 vr_fmt="PWmat", lattice_check_trigger=None,
                 comment: str = "vr"):
        # TODO 如果给传了 filename, ①检查(/提示) vr_fmt ②read
        self._latticeflag   = None
        self.lattice        = None
        self.mesh           = None
        self.comment = comment
        self.vr_fmt = vr_fmt
        self.lattice_check_trigger = lattice_check_trigger

        if filename is not None:
            self.read(filename=filename, vr_fmt=vr_fmt)

        if lattice  is not None:
            self._latticeflag = f"({self.comment}) (init parameter)"
            self.lattice: Lattice = lattice
        if mesh     is not None: 
            self.mesh: np.ndarray = mesh

    def __setattr__(self, name: str, value) -> None:
        if (name == "lattice") and (value is not None):
            assert type(value) is Lattice, "type of self.lattice is not Lattice"
            lattice: Lattice = value
            if self.lattice_check_trigger is not None:
                self.lattice_check_trigger(lattice)
            if self._latticeflag is not None:
                self.__dict__[name] = lattice.copy(self._latticeflag)
            else:
                self.__dict__[name] = lattice.copy(f"({self.comment}) (via setattr)")
        else:
            self.__dict__[name] = value

    @property
    def n123(self):
        assert self.mesh is not None, "the current `mesh` attribute in this instance is `None`"
        return np.array(self.mesh.shape, dtype=INTEGER)

    def revise_keyattribute_type(self):
        # TODO 看哪些定义哪些没定义
        # TODO 检查类型, 修改类型
        pass

    def read(self, filename, vr_fmt="PWmat"):
        self.filename = os.path.abspath(filename)
        self.vr_fmt = vr_fmt
        with open(filename, "br") as io:
            n1, n2, n3, nnodes = read_fortran_binary_block(io, INTEGER_IN)
            assert (n1*n2*n3) % nnodes == 0, "`n1*n2*n3` is not divisible by `nnodes`"
            self.nnodes = nnodes

            AL = read_fortran_binary_block(io, REAL_8)
            AL = np.reshape(AL, (3, 3))
            self._latticeflag = f"({self.comment}) read fromfile({self.filename})"
            self.lattice = Lattice(AL, self.fmt2unit[vr_fmt], 
                                   fromwhere=self._latticeflag)

            self.mesh = np.zeros(n1*n2*n3, dtype=REAL_8)
            nr = n1*n2*n3//nnodes
            for i in range(0, nnodes):
                self.mesh[i*nr: (i+1)*nr] = read_fortran_binary_block(io, REAL_8)[0:nr]
            self.mesh = np.reshape(self.mesh, (n1, n2, n3))
        
        return self

    def write(self, filename: str, fmt="PWmat", 
            nnodes: int | None = None, nnodes_base: int = 1):
        """
        nnodes will be converted into `INTEGER` when writing into file
        """
        if nnodes is None:
            nnodes = ceil(self.mesh.size / (128*1024*1024))
            if nnodes < 1:
                nnodes = 1
            while self.mesh.size % nnodes != 0 or nnodes % nnodes_base != 0:
                nnodes += 1
        # {mnodes} processors, {nnodes} vr slices
        # in Escan pot_input, 
        #     mod(mnodes, nnodes) == 0, if mnodes > nnodes => nnodes is 2^n
        #     mod(nnodes, mnodes) == 0, if nnodes > mnodes
        # but mnodes is unknown in potpatcch
        # 1. when nnodes is     specified, check mnodes 
        # 2. when nnodes is not specified, another while loop
        # but mnodes is unknown in potpatcch
        if fmt == "Escan":
            print("be careful of the relation between Escan parallel nnodes and vr slice nnodes")
        #     while self.mesh.size % nnodes != 0 and (nnodes & (nnodes - 1)): # nnodes is 2^n
        #         nnodes += 1
        assert self.mesh.size % nnodes == 0, \
            f"vr.mesh.size({self.mesh.size}) can't be devided evenly by nnodes({nnodes})"
        assert nnodes % nnodes_base == 0, \
            f"nnodes({nnodes}) is not divisible by nnodes_base({nnodes_base})"
        assert self.mesh.size / nnodes <= (128*1024*1024), \
            f"vr.mesh.size({self.mesh.size}) / nnodes({nnodes}) is larger than 128M"

        with open(filename, "bw") as io:
            # meta data
            write_fortran_binary_block(
                io, np.array(self.mesh.shape, dtype=INTEGER_OUT), INTEGER_OUT(nnodes))
            # AL data
            AL = self.lattice.in_unit(self.fmt2unit[fmt])
            AL1d = np.reshape(AL, AL.size)
            write_fortran_binary_block(io, AL1d)
            
            nr_n = self.mesh.size // nnodes
            mesh1d = np.reshape(self.mesh, self.mesh.size)
            for i in range(nnodes):
                write_fortran_binary_block(io, mesh1d[i*nr_n:(i+1)*nr_n])

    def __mul__(self, magnification) -> 'VR':
        """
        is used to make supercell
        self * magnification
        magnification example: ``(4,4,4)``
        """
        for i in list(magnification):
            assert isinstance(i, (int, np.integer)), "the magnification is not (all) integer "
        assert len(magnification) == 3, f"length of magnification must be 3, but not {len(magnification)}"
        
        lattice = self.lattice * magnification
        mesh = np.tile(self.mesh, magnification)
        return VR(lattice=lattice, mesh=mesh, vr_fmt=self.vr_fmt)
        
    def __rmul__(self, magnification) -> 'VR':
        return self.__mul__(magnification)

    def __imul__(self, magnification) -> 'VR':
        vr = self.__mul__(magnification)
        self.lattice = vr.lattice
        self.mesh    = vr.mesh
        return self


class AtomConfig():
    """
    atom.config
    """
    registered_fmt = ["PWmat", "Escan"]
    fmt2unit = dict(zip(registered_fmt, ["angstrom", "atomic unit"]))
    
    def __init__(self, filename=None, natoms=None, lattice: Lattice | None = None,
                 itypes=None, positions=None, moves=None, 
                 atoms_fmt="PWmat", lattice_check_trigger=None, 
                 comment: str = "atomconfig") -> None:
        self._latticeflag   = None
        self.lattice        = None
        self.natoms         = None
        self.itypes         = None
        self.positions      = None
        self.moves          = None

        self.comment = comment
        self.lattice_check_trigger = lattice_check_trigger
        self.atoms_fmt = atoms_fmt

        if filename is not None:
            self.read(filename=filename, fmt=atoms_fmt)

        # TODO 像 VR 一样注释类型
        if lattice        is not None: 
            self._latticeflag = f"({self.comment}) (init parameter)"
            self.lattice: Lattice = lattice 
        if natoms         is not None: self.natoms         = natoms        
        if itypes         is not None: self.itypes         = itypes   
        if positions      is not None: self.positions      = positions
        if moves          is not None: self.moves          = moves    

    def __setattr__(self, name: str, value) -> None:
        if (name == "lattice") and (value is not None):
            assert type(value) is Lattice, "type of self.lattice is not Lattice"
            lattice: Lattice = value
            if self.lattice_check_trigger is not None:
                self.lattice_check_trigger(lattice)
            if self._latticeflag is not None:
                self.__dict__[name] = lattice.copy(self._latticeflag)
            else:
                self.__dict__[name] = lattice.copy(f"({self.comment}) (via setattr)")
        else:
            self.__dict__[name] = value

    def revise_atomsposition(self):
        """
        revise fractional atom positions so that they are in the range [0,1)
        who cares: [X]make supercell([o]modify), [X]read, [○]write, [○]manually construct, [√]patch, [√]revise keyattribute
        """
        self.positions %= 1

    def sort_atomposition(self):
        sorted_indx = np.lexsort(self.positions.T)
        self.itypes     = self.itypes[sorted_indx]
        self.positions  = self.positions[sorted_indx]
        self.moves      = self.moves[sorted_indx]

    def revise_keyattribute_type(self):
        # TODO 看哪些定义哪些没定义
        # TODO 检查类型, 修改类型, 这对 write 很重要
        pass

    def read(self, filename: str, fmt: str = "PWmat") -> np.ndarray:
        self.filename = os.path.abspath(filename)
        with open(filename, "r") as io:
            natoms = int(io.readline().split()[0])
            self.natoms = natoms

            if fmt == "PWmat": 
                io.readline()
            AL = np.zeros((3, 3))
            for i in range(3):
                AL[i] = np.array([REAL_8(i) for i in io.readline().split()[0:3]])
            self._latticeflag = f"({self.comment}) read fromfile({self.filename})"
            self.lattice = Lattice(AL, self.fmt2unit[fmt], 
                                   fromwhere=self._latticeflag)

            if fmt == "PWmat": 
                io.readline()
            itype_list = np.zeros(natoms, dtype=INTEGER_IN)
            position_list = np.zeros((natoms, 3), dtype=REAL_8)
            move_list = np.zeros((natoms, 3), dtype=INTEGER_IN)
            for i in range(natoms):
                s = io.readline().split()[0:7]
                itype_list[i] = INTEGER_IN(s[0])
                position_list[i] = np.array([REAL_8(i) for i in s[1:4]])
                move_list[i] = np.array([INTEGER_IN(i) for i in s[4:7]])
            self.itypes     = itype_list
            self.positions  = position_list
            self.moves      = move_list
        
        return self

    def write(self, filename: str, comment: str | None = None, fmt: str = "PWmat"):
        with open(filename, "w") as io:
            if comment is not None:
                comment = comment.replace("\n", " ")
                io.write(f"      {self.natoms} atoms !! {comment}\n")
            else:
                io.write(f"      {self.natoms} atoms\n")

            if fmt == "PWmat": 
                io.write(" Lattice vector (Angstrom)\n")
            AL = self.lattice.in_unit(self.fmt2unit[fmt])
            for i in range(3):
                io.write("%20.14f%20.14f%20.14f\n" % tuple(AL[i]))

            if fmt == "PWmat": 
                io.write(" Position, move_x, move_y, move_z\n")
            for i in range(self.natoms):
                if fmt == "PWmat":
                    io.write("%6d %20.17f %20.17f %20.17f    %1d %1d %1d\n" % 
                             (self.itypes[i], *self.positions[i], *self.moves[i]))
                elif fmt == "Escan":
                    io.write("%6d %20.17f %20.17f %20.17f    0 1\n" % 
                             (self.itypes[i], *self.positions[i]))

    def __mul__(self, magnification) -> 'AtomConfig' :
        """
        is used to make supercell
        self * magnification
        magnification example: ``(4,4,4)``
        """
        for i in list(magnification):
            assert isinstance(i, (int, np.integer)), "the magnification is not (all) integer "
        assert len(magnification) == 3, f"length of magnification must be 3, but not {len(magnification)}"
        
        nrepeat = np.prod(magnification)
        natoms = self.natoms * nrepeat
        lattice = self.lattice * magnification
        atoms_itype = np.tile(self.itypes, nrepeat)
        atoms_position = np.tile(self.positions, (nrepeat, 1))
        atoms_move = np.tile(self.moves, (nrepeat, 1))
        for i in range(magnification[0]):
            for j in range(magnification[1]):
                for k in range(magnification[2]):
                    indx = self.natoms * (i*magnification[1]*magnification[2] + j*magnification[2] + k)
                    atoms_position[indx:indx+self.natoms] += np.array([i, j, k])
        atoms_position /= magnification
        return AtomConfig(natoms=natoms, lattice=lattice,
                          itypes=atoms_itype, positions=atoms_position, moves=atoms_move,
                          atoms_fmt=self.atoms_fmt)  # TODO 这个格式重要吗 (VR也是)
        
    def __rmul__(self, magnification) -> 'AtomConfig' :
        return self.__mul__(magnification)

    def __imul__(self, magnification) -> 'AtomConfig' :
        atomconfig = self.__mul__(magnification)
        self.natoms         = atomconfig.natoms
        self.lattice        = atomconfig.lattice
        self.itypes    = atomconfig.itypes
        self.positions = atomconfig.positions
        self.moves     = atomconfig.moves
        return self


class VATOM():
    """OUT.VATOM
    The unit of potential is eV .
    """
    def __init__(self, filename=None, comment: str = "unkown vatom") -> None:
        if filename is not None:
            self.read_vatom(filename)

        self.comment    = comment

    def read_vatom(self, filename) -> None:
        self.filename = os.path.abspath(filename)
        with open(filename, "r") as io:
            natoms = int(io.readline().split()[0])
            self.natoms = natoms

            itype_list      = np.zeros(natoms, dtype=INTEGER)
            position_list   = np.zeros((natoms, 3), dtype=REAL_8)
            vatom_list      = np.zeros((natoms), dtype=REAL_8)
            for i in range(natoms):
                s = io.readline().split()[0:5]
                itype_list[i]       = INTEGER(s[0])
                position_list[i]    = np.array([REAL_8(i) for i in s[1:4]])
                vatom_list[i]       = REAL_8(s[4])
            self.itypes     = itype_list
            self.positions  = position_list
            self._vatoms     = vatom_list / HA  # convert to a.u.

        return self
    
    @property
    def vatoms(self):  # default unit is eV (Consistent with the OUT.VATOM)
        return self._vatoms * HA 

    @property
    def vatoms_ev(self): 
        return self._vatoms * HA 

    @property
    def vatoms_au(self):
        return self._vatoms 

    def revise_atomsposition(self):
        """
        revise fractional atom positions so that they are in the range [0,1)
        who cares: [√]write_diffvatom
        """
        self.positions %= 1

class EIGEN():
    """
    OUT.EIGEN
    this class is only for reading propose, but not for constructing from 
    `eig=EIGEN();eig.eigenvals=...` or `EIGEN(eigenvals=...)`

    Attributes:
        eigenvals[islda, ikpt, iband]
        kpoints[ikpt, 3]
        weighkpt[ikpt]

    TODO too many dimensions, override __getitem__, use like key word
    """
    def __init__(self, filename=None, comment: str = "unkown eigen") -> None:
        if filename is not None:
            self.read_eigen(filename=filename)
        
        self.comment    = comment

    def read_eigen(self, filename: str) -> None:
        self.filename = os.path.abspath(filename)
        with open(filename, "br") as io:
            self.islda, self.nkpt, self.nband, self.nref_tot_8, self.natom, self.nnodes = read_fortran_binary_block(io, INTEGER_IN)
            self.eigenvals   = np.zeros((self.islda, self.nkpt, self.nband), dtype=REAL_8)
            self.kpoints     = np.zeros((self.nkpt, 3),                      dtype=REAL_8)
            self.weighkpt    = np.zeros((self.nkpt),                         dtype=REAL_8)
            
            for iislda in range(self.islda):
                for ikpt in range(self.nkpt):
                    ak = self.kpoints[ikpt, :]
                    iislda_, ikpt_, self.weighkpt[ikpt], ak[0], ak[1], ak[2] = \
                        read_fortran_binary_block_varioustype(
                            io, INTEGER_IN, INTEGER_IN, REAL_8, REAL_8, REAL_8, REAL_8)
                    self.eigenvals[iislda, ikpt, :] = \
                        read_fortran_binary_block(io, REAL_8)

        return self


class MaterialSystemInfo():
    """
    TODO __repr__, verbose
    """
    def __init__(self, 
                 lattice: Lattice = None, 
                 atomconfig: AtomConfig | None = None,  atoms_filename: str | None = None, atoms_fmt: str = "PWmat",
                 vr: VR | None = None,                  vr_filename:    str | None = None, vr_fmt: str = "PWmat", 
                 vatom: VATOM | None = None,            vatom_filename: str | None = None, 
                 eigen: EIGEN | None = None,            eigen_filename: str | None = None, 
                 charge=None, charge_pos=None, epsilon=None, 
                 comment: str = "unkown MaterialSystemInfoSummary") -> None:
        """
        self.lattice can take the data from vr/atomconfig modified by the `lattice_check_trigger=check_modify_lattice`
        self.charge is the periodic charge in unit cell. `-1` means 1 electron
        """
        if np.all(charge_pos):
            assert np.shape(charge_pos) == np.zeros(3).shape

        self._latticeflag = None
        self.lattice    = None
        self.comment    = comment

        self.charge     = charge
        self.charge_pos = charge_pos
        self.epsilon    = epsilon

        self.vr         = VR(lattice_check_trigger=self.check_modify_lattice)
        self.atomconfig = AtomConfig(lattice_check_trigger=self.check_modify_lattice)
        self.vatom      = VATOM()
        self.eigen      = EIGEN()
        
        if atoms_filename is not None:
            self.atomconfig.read(atoms_filename, atoms_fmt)
        if vr_filename    is not None:
            self.vr.read(vr_filename, vr_fmt)
        if vatom_filename is not None:
            self.vatom.read_vatom(vatom_filename)
        if eigen_filename is not None:
            self.eigen.read_eigen(eigen_filename)

        if lattice      is not None:
            self._latticeflag = f"({self.comment}) (init parameter)"
            self.lattice: Lattice = lattice 
        if atomconfig   is not None: self.atomconfig    = atomconfig
        if vr           is not None: self.vr            = vr
        if vatom        is not None: self.vatom         = vatom
        if eigen        is not None: self.eigen         = eigen

    def __setattr__(self, name: str, value) -> None:
        if name in ["vr", "atomconfig"]:
            self.__dict__[name] = value
            getattr(self, name).lattice_check_trigger = self.check_modify_lattice
            getattr(self, name).lattice_check_trigger(getattr(self, name).lattice)
        elif (name == "lattice") and (value is not None):
            assert type(value) is Lattice, "type of self.lattice is not Lattice"
            if (verbose := False):
                print(f"setting lattice({value.fromwhere}) in MaterialSystemInfo")
            lattice: Lattice = value
            if self._latticeflag is not None:
                _latticeflag = self._latticeflag
            else:
                _latticeflag = f"({self.comment}) (via setattr)"
            self.__dict__["lattice"]            = lattice.copy(_latticeflag)
            self.atomconfig.__dict__["lattice"] = self.lattice.copy(f"({self.atomconfig.comment}) (via setattr)")
            self.vr.__dict__["lattice"]         = self.lattice.copy(f"({self.vr.comment}) (via setattr)")
        else:
            self.__dict__[name] = value

    def check_modify_lattice(self, lattice: Lattice):
        if lattice is None:
            return 
        
        if self.lattice is None:
            self._latticeflag = f"({self.comment}) (via lattice_trigger)"
            self.lattice = lattice
        else:
            if not self.lattice == lattice: 
                s1 = " -> ".join(str(i) for i in self.lattice.fromwhere)
                s2 = " -> ".join(str(i) for i in lattice.fromwhere)
                print(
                    "lattice vector conflict with existing AL",
                    "old AL (angstrom): ",
                    f"{s1}",
                    f"{self.lattice.in_unit('angstrom')}",
                    "new AL (angstrom):",
                    f"{s2}",
                    f"{lattice.in_unit('angstrom')}",
                    sep="\n"
                )
