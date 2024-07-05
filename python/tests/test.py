import argparse

import numpy as np
import pyeitsolver
from icecream import ic


def save_potentials(mesh_potentials, output_filename, mesh_filename, nodes_count):
    # Read and copy the mesh file
    with open(mesh_filename, "r") as mesh_file:
        mesh_content = mesh_file.readlines()

    # Open the output file
    with open(output_filename, "w") as output_file:
        # Write the mesh content to the output file
        for line in mesh_content:
            output_file.write(line)

        # Append the solutions in the specified format
        for pattern_number, potentials in enumerate(mesh_potentials):
            output_file.write('$NodeData\n1\n"Electric Potential"\n1\n0.0\n3\n')
            output_file.write(f"{pattern_number}\n1\n{nodes_count}\n")
            enumerated_data = [f"{j + 1}\t{val}" for j, val in enumerate(potentials)]
            output_file.write("\n".join(enumerated_data))
            output_file.write("\n$EndNodeData\n")


def main(mesh_file, currents_file, output_file):
    solver = pyeitsolver.EitSolver(mesh_file, currents_file)
    ic(solver)
    ic(solver.info)

    conductivities = np.ones(solver.info["nodes_count"]) * 0.3810
    ic(conductivities.shape)

    electrode_solve_info, electrode_potentials = solver.solve_forward_problem(
        conductivities
    )

    ic(electrode_solve_info)
    ic(electrode_potentials[0])
    ic(electrode_potentials.shape)

    mesh_solve_info, mesh_potentials = solver.solve_forward_problem(
        np.ones(solver.info["nodes_count"]) * 0.3810, mesh_potentials=True
    )
    ic(mesh_solve_info)
    ic(mesh_potentials)
    ic(mesh_potentials.shape)

    save_potentials(mesh_potentials, output_file, mesh_file, solver.info["nodes_count"])

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
