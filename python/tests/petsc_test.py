import argparse

import petsc4py
from petsc4py import PETSc
import numpy as np
import pyeitsolver
from icecream import ic
from mesh_io import save_complex_potentials
from scipy.sparse import coo_array, diags

def main(mesh_file, currents_file, output_file):
    solver = pyeitsolver.EitComplexSolver(mesh_file, currents_file)
    b_np = ic(solver.get_currents_vectors())
    row, col, data = ic(
        solver.getCOO_formatted_stiffness(np.ones(solver.info["nodes_count"]) * 0.3810)
    )
    A_triangular = ic(
        coo_array(
            (data, (row, col)),
            shape=(solver.info["nodes_count"], solver.info["nodes_count"]),
        )
    )
    A_coo = ic(A_triangular + A_triangular.T - diags(A_triangular.diagonal()))

    # Initialize PETSc
    petsc4py.init()

    # Convert COO matrix to CSR format
    A_csr = ic(A_coo.tocsr())

    # Convert CSR to PETSc Mat
    A = PETSc.Mat().createAIJWithArrays(size=(solver.info["nodes_count"], solver.info["nodes_count"]), csr=(A_csr.indptr,  A_csr.indices, A_csr.data) )

    # Set symmetric option and assemble the matrix
    # A.setOption(PETSc.Mat.Option.SYMMETRIC, True)
    A.assemble()

    # A = ic(A_triangular + A_triangular.T - diags(A_triangular.diagonal()))


    # Create linear solver context
    ksp = PETSc.KSP().create(PETSc.COMM_WORLD)
    ksp.setOperators(A)
    ksp.setFromOptions()

    potentials = np.zeros(
        (solver.info["currents_count"], solver.info["nodes_count"])
    ).astype(np.complex128)

    for i in range(solver.info["currents_count"]):
        # Convert numpy array to PETSc Vec
        b = PETSc.Vec().createWithArray(b_np[i])

        # Create solution vector
        x = b.duplicate()

        # Solve linear system
        ksp.solve(b, x)

        potentials[i] = x.getArray() * solver.current

    ic(potentials)

    save_complex_potentials(
        np.real(potentials),
        np.imag(potentials),
        output_file,
        mesh_file,
        solver.info["nodes_count"],
    )

    print(f"Saved potentials to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run EIT solver with given mesh and currents files."
    )
    parser.add_argument("mesh_file", type=str, help="Path to the mesh file")
    parser.add_argument("currents_file", type=str, help="Path to the currents file")
    parser.add_argument("output_file", type=str, help="Path to the output file")

    args = parser.parse_args()
    main(args.mesh_file, args.currents_file, args.output_file)
