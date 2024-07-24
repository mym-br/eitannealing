import pathlib

from ._complex import EitComplexSolver as EitComplexCppSolver
from .validation import validate_currents

if not hasattr(EitComplexCppSolver, "solve_forward_problem"):
    import numpy as np
    from scipy.sparse import coo_array, diags, linalg

    try:
        import petsc4py

        use_petsc = True
    except ModuleNotFoundError:
        from scipy.sparse import linalg

        use_petsc = False

DEFAULT_CURRENTS_COUNT = 32


class EitComplexSolver(EitComplexCppSolver):
    def __init__(
        self,
        mesh_filename: str,
        currents_filename: str,
        currents_count: int = DEFAULT_CURRENTS_COUNT,
    ):
        """
        Initialize the EitComplexSolver with mesh and current filenames.

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

    @property
    def type(self):
        if hasattr(EitComplexCppSolver, "solve_forward_problem"):
            return "PETSc KSP"
        elif use_petsc:
            return "petsc4py KSP"
        else:
            return "scipy.sparse.linalg"

    def solve_forward_problem(self, conds, mesh_potentials=False):

        if self.type == "PETSc KSP":
            return super().solve_forward_problem(conds, mesh_potentials=mesh_potentials)
        else:
            b_np = self.get_currents_vectors()
            row, col, data = self.getCOO_formatted_stiffness(conds)

            A_triangular = coo_array(
                (data, (row, col)),
                shape=(self.info["nodes_count"], self.info["nodes_count"]),
            )
            A_coo = A_triangular + A_triangular.T - diags(A_triangular.diagonal())

            potentials = np.zeros(
                (self.info["currents_count"], self.info["nodes_count"])
            ).astype(np.complex128)

            if self.type == "scipy.sparse.linalg":
                # solve using scipy linear algebra functions
                for i in range(self.info["currents_count"]):
                    x = linalg.spsolve(A_coo, b_np[i])
                    potentials[i] = x * self.current
            else:
                # solve using petsc
                petsc4py.init()

                # create PETSc Mat from CSR
                A_csr = A_coo.tocsr()
                A = petsc4py.PETSc.Mat().createAIJWithArrays(
                    size=(self.info["nodes_count"], self.info["nodes_count"]),
                    csr=(A_csr.indptr, A_csr.indices, A_csr.data),
                )
                A.assemble()

                # create linear solver context
                ksp = petsc4py.PETSc.KSP().create(petsc4py.PETSc.COMM_WORLD)
                ksp.setOperators(A)
                ksp.setFromOptions()

                for i in range(self.info["currents_count"]):
                    b = petsc4py.PETSc.Vec().createWithArray(b_np[i])
                    x = b.duplicate()
                    ksp.solve(b, x)
                    potentials[i] = x.getArray() * self.current

            return potentials
            # raise NotImplementedError("Eit complex forward solving requires scipy or PETSc package")
