import argparse

import numpy as np
import pyeitsolver
from icecream import ic
from mesh_io import save_complex_potentials
from scipy.sparse import coo_array, diags, linalg


def add_matrix_capacitances(
    stiffness, first_electrode_idx, electrodes_count, freq, capacitance
):
    jwc = 1j * 2 * np.pi * freq * capacitance
    for i in range(first_electrode_idx, first_electrode_idx + electrodes_count):
        stiffness[i, i] += jwc
    return stiffness


def main(mesh_file, currents_file, output_file):
    solver = pyeitsolver.EitSolver(mesh_file, currents_file)
    b = ic(solver.get_currents_vectors())
    row, col, data = ic(
        solver.getCOO_formatted_stiffness(np.ones(solver.info["nodes_count"]) * 0.3810)
    )

    # add 1 to the last diagonal element
    row = np.append(row, solver.info["nodes_count"] - 1)
    col = np.append(col, solver.info["nodes_count"] - 1)
    data = np.append(data, 1.0)
    A_triangular = ic(
        coo_array(
            (data, (row, col)),
            shape=(solver.info["nodes_count"], solver.info["nodes_count"]),
        )
    )

    A = ic(A_triangular + A_triangular.T - diags(A_triangular.diagonal()))

    A_complex = ic(
        add_matrix_capacitances(
            A.astype(np.complex128),
            solver.info["first_electrode_idx"],
            solver.info["electrodes_count"],
            275000,
            80e-12,
        )
    )

    potentials = np.zeros(
        (solver.info["currents_count"], solver.info["nodes_count"])
    ).astype(np.complex128)

    for i in range(solver.info["currents_count"]):
        x = ic(linalg.spsolve(A_complex, np.append(b[i], 0)))
        potentials[i] = x * solver.current

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
