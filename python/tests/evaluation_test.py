import argparse
import math
import os
from pathlib import Path

import h5py
import numpy as np
import pyeitsolver
from icecream import ic


def read_mesh_values(file_path):
    with open(file_path, "r") as file:
        lines = file.readlines()

    # Find the index of the $NodeData line
    node_data_index = next(i for i, line in enumerate(lines) if "$NodeData" in line)

    # Drop all lines before $NodeData + 8
    relevant_lines = lines[node_data_index + 8 :]

    # The first line in relevant_lines should be the number of values
    values_count = int(relevant_lines[0].strip())

    # Extract the values, ignoring the indices
    values = [float(line.split()[1]) for line in relevant_lines[1 : 1 + values_count]]

    return np.array(values)


def read_measurements_values(file_path):
    measurements_array = np.genfromtxt(file_path, delimiter="   ")
    size_sqrt = int(math.sqrt(measurements_array.shape[0]))
    measurements_reshaped = np.reshape(measurements_array, (size_sqrt, size_sqrt))
    # Subtract the last value of each row from each element in that row using broadcasting
    measurements_adjusted = (
        measurements_reshaped - measurements_reshaped[:, -1][:, np.newaxis]
    )

    return measurements_adjusted


def mean_squared_error(array1, array2):
    # Ensure both arrays have the same shape
    if array1.shape != array2.shape:
        raise ValueError("Input arrays must have the same shape")

    # Compute the mean squared error
    mse = np.mean((array1 - array2) ** 2)
    return mse


def main(mesh_file, currents_file, results_folder, data_folder, output_file):
    solver = ic(pyeitsolver.EitSolver(mesh_file, currents_file))
    ic(solver.info)

    # Find all .msh files in results_folder
    mesh_files = list(Path(results_folder).glob("*.msh"))

    # Create corresponding .txt paths in data_folder
    measurements_files = [
        Path(os.path.join(data_folder, mesh_file.with_suffix(".txt").name))
        for mesh_file in mesh_files
    ]

    mse_matrix = np.zeros((len(mesh_files), len(measurements_files)))
    for mesh_idx, mesh_file in enumerate(mesh_files):
        conductivities = read_mesh_values(mesh_file)
        _, calculated_potentials = solver.solve_forward_problem(conductivities)
        for measurement_idx, measurement_file in enumerate(measurements_files):
            measured_potentials = read_measurements_values(measurement_file)
            mse_matrix[mesh_idx, measurement_idx] = ic(
                mean_squared_error(calculated_potentials, measured_potentials)
            )
    ic(mse_matrix)

    with h5py.File(output_file, "w") as f:
        dset = f.create_dataset("mses", data=mse_matrix)
        dset.attrs["experiments"] = [mesh.name for mesh in mesh_files]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run EIT solver with given mesh and currents files."
    )
    parser.add_argument("mesh_file", type=str, help="Path to the mesh file")
    parser.add_argument("currents_file", type=str, help="Path to the currents file")
    parser.add_argument("results_folder", type=str, help="Path to the results folder")
    parser.add_argument("data_folder", type=str, help="Path to the data folder")
    parser.add_argument("output_file", type=str, help="Path to the output file")

    args = parser.parse_args()
    main(
        args.mesh_file,
        args.currents_file,
        args.results_folder,
        args.data_folder,
        args.output_file,
    )
