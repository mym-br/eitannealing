from __future__ import annotations

from ._core import __doc__, __version__
from .complex_wrapper import EitComplexSolver
from .core_wrapper import EitSolver

__all__ = ["__doc__", "__version__", "EitSolver", "EitComplexSolver"]
