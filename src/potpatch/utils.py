from textwrap import dedent
from time import perf_counter

import numpy as np


# =============================================
# fortran binary io
# =============================================
from potpatch.datatype import REAL_8, INTEGER, INTEGER_4


def read_fortran_binary_block_arraytype(io, dtype: np.dtype):
    # TODO 用 memoryview 实现 multi data type
    len_Byte = np.fromfile(io, INTEGER_4, 1)[0]
    sizeofd = np.dtype(dtype).itemsize
    assert len_Byte % sizeofd == 0, dedent(f"\
        the length({len_Byte}) of fortran data block is not divisible \
        by the length({sizeofd}) of given dtype")

    # "begin read"
    count = len_Byte // sizeofd
    data = np.fromfile(io, dtype, count=count)
    # "end read"

    assert len_Byte == np.fromfile(io, INTEGER_4, 1)[0], "\
        there is something wrong in Fortran binary file format"
    return data


def read_fortran_binary_block_varioustype(io, *dtypes):
    # TODO 用 memoryview 实现 multi data type
    len_Byte = np.fromfile(io, INTEGER_4, 1)[0]
    sum_sizeofd = sum(
        (np.dtype(dtype).itemsize)
        for dtype in dtypes
    )
    assert sum_sizeofd == len_Byte, dedent(f"\
        the length({len_Byte}) of fortran data block does not match \
        the sum of length({sum_sizeofd}) of given dtypes")

    # "begin read "
    data = list(range(len(dtypes)))
    for (i, dtype) in enumerate(dtypes):
        data[i] = np.fromfile(io, dtype, count=1)
    data = tuple(data)
    # "end read"

    assert len_Byte == np.fromfile(io, INTEGER_4, 1)[0], "\
        there is something wrong in Fortran binary file format"
    return data


def read_fortran_binary_block(io, dtype: np.dtype):
    return read_fortran_binary_block_arraytype(io, dtype)


def write_fortran_binary_block(io, *datas):
    """
    datas: tuple of numpy.ndarray|numpy.generic
    """
    len_Byte = 0
    for data in datas:
        # TODO 或许这里可以换成 warnings.warn
        assert isinstance(data, np.ndarray) or isinstance(data, np.generic), \
            "stroooongly recommend the data in datas is the instance of "\
            "numpy.ndarray|numpy.generic, so you can control the type of "\
            "data, which is important in fortran binary file read/write."
        sizeofd = data.dtype.itemsize
        len_Byte += sizeofd * data.size
    len_Byte = INTEGER_4(len_Byte)

    len_Byte.tofile(io)
    for data in datas:
        if not (isinstance(data, np.ndarray) or isinstance(data, np.generic)):
            data = np.array(data)
        data.tofile(io)
    len_Byte.tofile(io)


# =============================================
# numerical method
# =============================================
def simpson(dx: REAL_8, y: np.ndarray) -> np.ndarray:
    """
    input list length must be odd, The corresponding index::

        input:  0 1 2 3 4 5 6
        output: 0   1   2   3

    Example:
    >>> simpson(0.5, np.log([1,1.5,2,2.5,3]))
    array([0,0.3858...,0.9094...])

    i.e.::

        input:  ln(1)     ln(1.5) ln(2)     ln(2.5) ln(3)
        output: 0                 0.3858...         0.9094...

        ( ∫₁² lnx dx ≈ 0.38629436...
          ∫₂³ lnx dx ≈ 0.90954250... )

    """
    assert len(y) % 2 == 1, "length of y must be odd"
    intglength = (len(y)+1)//2
    integrate = np.ndarray(intglength)
    integrate[0] = 0
    for i in range(1, intglength):
        integrate[i] = (y[2*i-2] + 4*y[2*i-1] + y[2*i]) * dx / 3
    return integrate


# =============================================
# basic classes
# =============================================
class NameTuple():
    def __init__(self, **kwargs) -> None:
        # print(kwargs)
        self.__slots__ = tuple(kwargs.keys())
        for attr in self.__slots__:
            setattr(self, attr, kwargs[attr])

    def __repr__(self) -> str:
        s = [f"{attr} = {getattr(self, attr)}" for attr in self.__slots__]
        s = ", ".join(s)
        return f"({s})"

    def __str__(self) -> str:
        return self.__repr__()


# =============================================
# functions 
# =============================================
def gen_counter():
    count = 0

    def plusplus():
        nonlocal count
        count += 1
        return count
    return plusplus


# =============================================
# decorators 
# =============================================
PRINT_TIMING_DEFAULT   = False
PRINT_TIMING_ALL_CLOSE = None
PRINT_TIMING_ALL_OPEN  = None


def timing(tprint=PRINT_TIMING_DEFAULT):
    assert PRINT_TIMING_ALL_CLOSE is None or PRINT_TIMING_ALL_OPEN is None
    
    def timing_decorator(func):
        def _timing(*args, **kwargs):
            nonlocal tprint
            t0 = perf_counter()
            results = func(*args, **kwargs)
            t1 = perf_counter()
            tprint = True  if PRINT_TIMING_ALL_OPEN  is True else tprint 
            tprint = False if PRINT_TIMING_ALL_CLOSE is True else tprint 
            if tprint:
                print("...", f"{func.__name__:20s} time cost: {t1 - t0} s")
            return results
        return _timing

    return timing_decorator
