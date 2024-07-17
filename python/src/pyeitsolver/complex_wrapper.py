import pathlib

from ._complex import EitComplexSolver as EitComplexCppSolver
if not hasattr(EitComplexCppSolver, "solve_forward_problem"):
    from scipy.sparse import coo_array, diags, linalg
    import numpy as np
    try:
        import petsc4py
        use_petsc = True
    except ModuleNotFoundError:
        from scipy.sparse import linalg

class EitComplexSolver(EitComplexCppSolver):
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
            
            A_triangular = coo_array((data, (row, col)), shape=(self.info["nodes_count"], self.info["nodes_count"]),)
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
                A = petsc4py.PETSc.Mat().createAIJWithArrays(size=(self.info["nodes_count"], self.info["nodes_count"]), csr=(A_csr.indptr,  A_csr.indices, A_csr.data) )
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
