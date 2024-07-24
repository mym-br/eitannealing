import pathlib

from ._core import EitSolver as EitCppSolver
from .validation import validate_currents

DEFAULT_CURRENTS_COUNT = 32


class EitSolver(EitCppSolver):
    def __init__(
        self,
        mesh_filename: str,
        currents_filename: str,
        currents_count: int = DEFAULT_CURRENTS_COUNT,
    ):
        """
        Initialize the EitSolver with mesh and current filenames.

        Parameters:
        mesh_filename (str): The path to the mesh file.
        currents_filename (str): The path to the current file.
        """
        # Perform some checks before calling the C++ constructor

        # Validate currents. TODO: retrieve electrodes count from mesh validation
        with open(currents_filename, "r", encoding="utf-8") as currents_file:
            validate_currents(currents_file, currents_count)

        # Call the C++ constructor
        super().__init__(mesh_filename, currents_filename)
