import numpy as np


# =============================================
# fortran binary io
# =============================================
from textwrap import dedent

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
    assert len_Byte < 2**31, \
        "the length of fortran binary data record block is too large"
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


def revise_epsilon(epsilon):
    """
    input epsilon is float or list[float]
    output epsilon is 2DArray[float]
    """
    if isinstance(epsilon, list):
        epsilon = np.array(epsilon)
    else:
        epsilon = float(epsilon)
        epsilon = np.diag([epsilon, epsilon, epsilon])
    return epsilon


# =============================================
# decorators 
# =============================================
from time import perf_counter
from functools import wraps

PRINT_TIMING_DEFAULT   = False
PRINT_TIMING_ALL_CLOSE = None
PRINT_TIMING_ALL_OPEN  = None


def timing(tprint=PRINT_TIMING_DEFAULT):
    assert PRINT_TIMING_ALL_CLOSE is None or PRINT_TIMING_ALL_OPEN is None
    
    def timing_decorator(func):
        @wraps(func)
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


# =============================================
# log component
# =============================================
from textwrap import indent
from datetime import datetime
import inspect
import re


def log(*message, end="\n",
        unadorned=False, log_level=1, pre_idt = 0, timestamp=False):
    """
    Notes
    -----
    - indent
    """
    adorn = ""
    if timestamp:
        ts = datetime.now().strftime("%Y-%m-%d %H:%M")
        ts = f"({ts})"

    if codecontext_has_n(r"#.*>silence") >= 1:
        return
    
    n_idt = max(0, codecontext_has_n(r"#.*?>log") - 1)
    idt = ' ' * (4 * (pre_idt + n_idt) + len(adorn))
    if unadorned:
        adorn = ''
        idt = ''
    for msg in message:
        print(indent(msg, idt), end=end)


def codecontext_has_n(mark: str) -> int:   
    """Count occurrences of a marker string in the call stack's code contexts.

    Parameters
    ----------
    mark : str
        String pattern to search for (will be regex-escaped). Special regex
        characters are automatically escaped.

    Returns
    -------
    : int
        Number of frames where the marker appears in code context lines.
        Returns 0 if no matches found or code contexts unavailable.

    Notes
    -----
    - Only checks the active line (frame.index) of each frame
    - `findall` but not `search`
    - Safe against None code_context (e.g. for builtin functions)
    """
    n = 0

    stack = inspect.stack()
    for frame in stack[2:]:
        code_context = frame.code_context
        if code_context:
            n += len(re.findall(mark, code_context[frame.index]))
    
    return n


def extract_logmark(mark: str, mark_f):
    """Extract and process log context from caller's code line using regex 
    patterns.
    
    Parameters
    ----------
    mark : str
        Regex pattern to search in log content
    mark_f : callable
        Function to process matched groups from mark pattern
        
    Returns
    -------
    any or None
        Processed result from mark_f if both patterns match, otherwise None
        
    Examples
    --------
    >>> def f():
        codecontext_get(r'(\d+)', lambda x: int(x)*2)
    >>> f()  # >log test123
    246
    """
    stack = inspect.stack()
    "stack[0]: inspect.stack()"
    "stack[1]: extract_logmark(...)"
    caller_frame = stack[2] 
    code_context = caller_frame.code_context
    code_line = code_context[caller_frame.index]

    # Match log pattern format: #...>log...
    log_expr_match = re.search(r'#.*>log(.*)$', code_line)
    if not log_expr_match:
        return None
        
    log_content = log_expr_match.group(1)
    
    mark_match = re.search(mark, log_content)
    if not mark_match:
        return None
        
    return mark_f(*mark_match.groups())
