import pathlib

from ._core import EitSolver as EitSolverCpp


class EitSolver(EitSolverCpp):
    def __init__(self, mesh_filename, currents_filename):
        """
        Initialize the EitSolver with mesh and current filenames.

        Parameters:
        mesh_filename (str): The path to the mesh file.
        currents_filename (str): The path to the current file.
        """
        # Perform some checks before calling the C++ constructor
        if pathlib.Path(mesh_filename).exists() == False:
            raise FileNotFoundError(
                "Input mesh file not found. The value was: {}".format(mesh_filename)
            )
        if pathlib.Path(currents_filename).exists() == False:
            raise FileNotFoundError(
                "Input currents file not found. The value was: {}".format(
                    currents_filename
                )
            )

        # Call the C++ constructor
        super().__init__(mesh_filename, currents_filename)
