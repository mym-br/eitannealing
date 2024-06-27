import argparse

import numpy as np
import pyeitsolver


def main(mesh_file, currents_file):
    problem_data = pyeitsolver.init(mesh_file, currents_file)
    print(problem_data)

    electrode_potentials = pyeitsolver.solve_forward_problem(
        np.ones(problem_data["coeffCount"]) * 0.3810
    )
    print(electrode_potentials[1])
    print(electrode_potentials[1][:32])
    print(electrode_potentials[1].shape)

    mesh_potentials = pyeitsolver.solve_full_forward_problem(
        np.ones(problem_data["coeffCount"]) * 0.3810
    )
    print(mesh_potentials)
    print(mesh_potentials.shape)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run EIT solver with given mesh and currents files."
    )
    parser.add_argument("mesh_file", type=str, help="Path to the mesh file")
    parser.add_argument("currents_file", type=str, help="Path to the currents file")

    args = parser.parse_args()
    main(args.mesh_file, args.currents_file)
