"""Validation rules and reports for :class:`MaterialSystemInfo`.

This module implements the checks exposed by
``potpatch.objects.MaterialSystemInfo.validate()``. Checks inspect the current
state without changing the material system; validate again after editing it.
"""

from dataclasses import dataclass
from collections.abc import Iterable
from numbers import Integral
from typing import TYPE_CHECKING, Literal

import numpy as np

from potpatch.objects import Lattice

if TYPE_CHECKING:
    from potpatch.objects import MaterialSystemInfo


CheckStatus = Literal["pass", "fail", "insufficient"]

PATCH_REQUIRED_FIELDS = (
    "lattice", "atomconfig.lattice", "atomconfig.natoms",
    "atomconfig.itypes", "atomconfig.positions", "atomconfig.moves",
    "vr.lattice", "vr.mesh",
)


@dataclass(frozen=True)
class CheckResult:
    """Result of one check on a ``MaterialSystemInfo`` instance."""

    name: str
    status: CheckStatus
    detail: str
    sources: tuple[str, ...] = ()


@dataclass(frozen=True)
class ValidationReport:
    """Snapshot of checks from one ``MaterialSystemInfo.validate()`` call."""

    checks: tuple[CheckResult, ...]

    @property
    def failures(self) -> tuple[CheckResult, ...]:
        return tuple(check for check in self.checks if check.status == "fail")

    @property
    def passed(self) -> bool:
        return bool(self.checks) and all(
            check.status == "pass" for check in self.checks
        )

    def raise_for_errors(self) -> None:
        if self.failures:
            details = "\n".join(
                f"{check.name}: {check.detail}" for check in self.failures
            )
            raise ValueError(f"Material system validation failed:\n{details}")


def _check_array(checks, name, value, ndim, tail_shape=None, integer=False):
    if value is None:
        return None
    try:
        array = np.asarray(value)
        shape_ok = array.ndim == ndim and (
            tail_shape is None or array.shape[-len(tail_shape):] == tail_shape
        )
        type_ok = isinstance(value, np.ndarray) and np.issubdtype(
            array.dtype, np.integer if integer else np.floating
        )
    except (TypeError, ValueError):
        shape_ok = type_ok = False
        array = None
    if not shape_ok or not type_ok:
        checks.append(CheckResult(name, "fail", "invalid shape or numeric type"))
        return None
    checks.append(CheckResult(name, "pass", f"shape {array.shape}"))
    return array


def _check_lattice(checks, name, lattice):
    if lattice is None:
        return None
    if not isinstance(lattice, Lattice):
        checks.append(CheckResult(name, "fail", "expected a Lattice"))
        return None
    if lattice.unit not in Lattice.registered_units:
        checks.append(CheckResult(name, "fail", f"unknown unit {lattice.unit!r}"))
        return None
    try:
        matrix = np.asarray(lattice.AL_AU, dtype=float)
        valid = (matrix.shape == (3, 3) and np.isfinite(matrix).all()
                 and np.linalg.matrix_rank(matrix) == 3)
    except (TypeError, ValueError, np.linalg.LinAlgError):
        valid = False
        matrix = None
    if not valid:
        checks.append(CheckResult(name, "fail", "expected a finite, nonsingular 3x3 lattice"))
        return None
    checks.append(CheckResult(name, "pass", "valid lattice"))
    return matrix


def _get_path(info, path):
    value = info
    for part in path.split("."):
        value = getattr(value, part, None)
        if value is None:
            break
    return value


def validate_material_system(
    info: "MaterialSystemInfo", required_paths: Iterable[str] = ()
) -> ValidationReport:
    """Check one ``MaterialSystemInfo`` without changing its data.

    ``required_paths`` names fields that must be present for the caller's
    operation. Other missing fields are allowed, while present fields and
    relationships between them are checked.
    """
    checks: list[CheckResult] = []

    for path in required_paths:
        present = _get_path(info, path) is not None
        checks.append(CheckResult(
            f"required.{path}", "pass" if present else "fail",
            "present" if present else "required for this operation but missing",
        ))

    atomconfig = info.atomconfig
    vr = info.vr
    lattice_sources = (
        ("lattice", info._lattice),
        ("atomconfig.lattice", getattr(atomconfig, "lattice", None)),
        ("vr.lattice", getattr(vr, "lattice", None)),
    )
    lattice_values = {
        name: _check_lattice(checks, name, value)
        for name, value in lattice_sources
    }
    for index, (left_name, left) in enumerate(lattice_sources):
        for right_name, right in lattice_sources[index + 1:]:
            if left_name == "lattice" and left is None:
                continue
            if left is None and right is None:
                continue
            name = f"{left_name} == {right_name}"
            left_matrix = lattice_values[left_name]
            right_matrix = lattice_values[right_name]
            sources = (
                " -> ".join(map(str, getattr(left, "fromwhere", []))),
                " -> ".join(map(str, getattr(right, "fromwhere", []))),
            )
            if left_matrix is None or right_matrix is None:
                checks.append(CheckResult(
                    name, "insufficient", "both valid lattices are needed", sources
                ))
                continue
            delta = float(np.max(np.abs(left_matrix - right_matrix)))
            checks.append(CheckResult(
                name, "pass" if delta < 1e-5 else "fail",
                f"maximum difference {delta:.6g} atomic units", sources,
            ))

    if vr is not None:
        mesh = _check_array(checks, "vr.mesh", vr.mesh, 3)
        if mesh is not None and not all(mesh.shape):
            checks.append(CheckResult(
                "vr.mesh.dimensions", "fail", "expected nonempty grid axes"
            ))

    if atomconfig is not None:
        natoms = atomconfig.natoms
        natoms_valid = False
        if natoms is not None:
            natoms_valid = (
                isinstance(natoms, Integral) and not isinstance(natoms, bool)
                and natoms >= 0
            )
            checks.append(CheckResult(
                "atomconfig.natoms", "pass" if natoms_valid else "fail",
                "nonnegative integer" if natoms_valid
                else "expected a nonnegative integer",
            ))
        for field, ndim, tail, integer in (
            ("itypes", 1, None, True),
            ("positions", 2, (3,), False),
            ("moves", 2, (3,), True),
        ):
            array = _check_array(
                checks, f"atomconfig.{field}", getattr(atomconfig, field),
                ndim, tail, integer,
            )
            if array is not None and field == "positions" and not np.isfinite(array).all():
                checks.append(CheckResult(
                    "atomconfig.positions.finite", "fail",
                    "positions contain non-finite values",
                ))
            if array is not None and natoms_valid:
                checks.append(CheckResult(
                    f"atomconfig.natoms == len({field})",
                    "pass" if len(array) == natoms else "fail",
                    f"natoms={natoms}, len({field})={len(array)}",
                ))

    vatom_natoms = getattr(info.vatom, "natoms", None)
    atom_natoms = getattr(atomconfig, "natoms", None)
    if vatom_natoms is not None:
        if atom_natoms is None:
            checks.append(CheckResult(
                "atomconfig.natoms == vatom.natoms", "insufficient",
                "atom configuration atom count is missing",
            ))
        else:
            checks.append(CheckResult(
                "atomconfig.natoms == vatom.natoms",
                "pass" if atom_natoms == vatom_natoms else "fail",
                f"atomconfig={atom_natoms}, vatom={vatom_natoms}",
            ))

    eigen_natoms = getattr(info.eigen, "natom", None)
    if eigen_natoms is not None:
        if atom_natoms is None:
            checks.append(CheckResult(
                "atomconfig.natoms == eigen.natom", "insufficient",
                "atom configuration atom count is missing",
            ))
        else:
            checks.append(CheckResult(
                "atomconfig.natoms == eigen.natom",
                "pass" if atom_natoms == eigen_natoms else "fail",
                f"atomconfig={atom_natoms}, eigen={eigen_natoms}",
            ))

    if info.charge is not None:
        try:
            charge = np.asarray(info.charge, dtype=float)
            valid = charge.shape == () and np.isfinite(charge).all()
        except (TypeError, ValueError):
            valid = False
        checks.append(CheckResult(
            "charge", "pass" if valid else "fail",
            "finite scalar" if valid else "expected a finite scalar",
        ))

    if info.charge_pos is not None:
        try:
            charge_pos = np.asarray(info.charge_pos, dtype=float)
            valid = charge_pos.shape == (3,) and np.isfinite(charge_pos).all()
        except (TypeError, ValueError):
            valid = False
        checks.append(CheckResult(
            "charge_pos", "pass" if valid else "fail",
            "finite fractional 3-vector" if valid else "expected a finite 3-vector",
        ))

    if info.epsilon is not None:
        try:
            epsilon = np.asarray(info.epsilon, dtype=float)
            valid = epsilon.shape in ((), (3, 3)) and np.isfinite(epsilon).all()
        except (TypeError, ValueError):
            valid = False
        checks.append(CheckResult(
            "epsilon", "pass" if valid else "fail",
            "finite scalar or 3x3 tensor" if valid
            else "expected a finite scalar or 3x3 tensor",
        ))

    return ValidationReport(tuple(checks))
