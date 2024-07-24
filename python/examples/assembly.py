import argparse
from pathlib import Path

import h5py
import numpy as np
import pyeitsolver
from icecream import ic
from mesh_io import save_potentials
from scipy.sparse import coo_array, diags, linalg


def main(mesh_file, currents_file, output_file):
    solver = pyeitsolver.EitSolver(mesh_file, currents_file)
    b = ic(solver.get_currents_vectors())
    row, col, data = ic(
        solver.getCOO_formatted_stiffness(np.ones(solver.info["nodes_count"]) * 0.3810)
    )
    A_triangular = ic(
        coo_array(
            (data, (row, col)),
            shape=(solver.info["nodes_count"] - 1, solver.info["nodes_count"] - 1),
        )
    )

    A = ic(A_triangular + A_triangular.T - diags(A_triangular.diagonal()))

    # Save matrix A to numpy format
    output_path = Path(output_file)
    A_output_file = output_path.with_stem(output_path.stem + "_A").with_suffix(".npy")
    np.save(A_output_file, A.toarray())

    print(f"Saved matrix A to {A_output_file}")

    # why does it not work?
    # ic(linalg.spsolve_triangular(A_triangular, b[0], lower=False))

    potentials = np.zeros((solver.info["currents_count"], solver.info["nodes_count"]))

    for i in range(solver.info["currents_count"]):
        x = ic(linalg.spsolve(A, b[i]))
        potentials[i] = np.insert(x, solver.info["ground_node"], 0.0) * solver.current

    ic(potentials)

    save_potentials(potentials, output_file, mesh_file, solver.info["nodes_count"])

    # Save data to HDF5 dataset
    hdf5_output_file = output_path.with_suffix(".hdf5")
    with h5py.File(hdf5_output_file, "w") as f:
        f.create_dataset("A", data=A.toarray())
        f.create_dataset("b", data=b)
        f.create_dataset("x", data=potentials)

    print(f"Saved potentials to {output_file} and dataset to {hdf5_output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run EIT solver with given mesh and currents files."
    )
    parser.add_argument("mesh_file", type=str, help="Path to the mesh file")
    parser.add_argument("currents_file", type=str, help="Path to the currents file")
    parser.add_argument("output_file", type=str, help="Path to the output file")

    args = parser.parse_args()
    main(args.mesh_file, args.currents_file, args.output_file)
