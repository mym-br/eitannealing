import argparse

import numpy as np
import pyeitsolver
from icecream import ic
from mesh_io import save_complex_potentials


def main(mesh_file, currents_file, output_file):
    solver = ic(pyeitsolver.EitComplexSolver(mesh_file, currents_file))
    ic(type(solver))
    potentials = ic(solver.solve_forward_problem(np.ones(solver.info["nodes_count"]) * 0.3810))
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
