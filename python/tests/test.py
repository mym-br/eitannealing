import argparse

import numpy as np
import pyeitsolver


def main(mesh_file, currents_file):
    solver = pyeitsolver.EitSolver(mesh_file, currents_file)
    print(solver)
    print(solver.info)

    conductivities = np.ones(solver.info["coeff_count"]) * 0.3810
    print(conductivities.shape)

    electrode_potentials = solver.solve_forward_problem(conductivities)

    print(electrode_potentials[0])
    print(electrode_potentials[1][0])
    print(electrode_potentials[1].shape)

    mesh_potentials = solver.solve_forward_problem(
        np.ones(solver.info["coeff_count"]) * 0.3810, mesh_potentials=True
    )
    print(mesh_potentials[0])
    print(mesh_potentials[1])
    print(mesh_potentials[1].shape)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run EIT solver with given mesh and currents files."
    )
    parser.add_argument("mesh_file", type=str, help="Path to the mesh file")
    parser.add_argument("currents_file", type=str, help="Path to the currents file")

    args = parser.parse_args()
    main(args.mesh_file, args.currents_file)
