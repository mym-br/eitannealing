import numpy as np


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


def save_complex_potentials(
    mesh_potentials_real,
    mesh_potentials_imag,
    output_filename,
    mesh_filename,
    nodes_count,
):
    # Read and copy the mesh file
    with open(mesh_filename, "r") as mesh_file:
        mesh_content = mesh_file.readlines()

    # Open the output file
    with open(output_filename, "w") as output_file:
        # Write the mesh content to the output file
        for line in mesh_content:
            output_file.write(line)

        for pattern_number, (potentials_real, potentials_imag) in enumerate(
            zip(mesh_potentials_real, mesh_potentials_imag)
        ):
            # Append the real part of the solutions in the specified format
            output_file.write('$NodeData\n1\n"Real Electric Potential"\n1\n0.0\n3\n')
            output_file.write(f"{pattern_number}\n1\n{nodes_count}\n")
            enumerated_real_data = [
                f"{j + 1}\t{val}" for j, val in enumerate(potentials_real)
            ]
            output_file.write("\n".join(enumerated_real_data))
            output_file.write("\n$EndNodeData\n")

            # Append the imaginary part of the solutions in the specified format
            output_file.write(
                '$NodeData\n1\n"Imaginary Electric Potential"\n1\n0.0\n3\n'
            )
            output_file.write(f"{pattern_number}\n1\n{nodes_count}\n")
            enumerated_imag_data = [
                f"{j + 1}\t{val}" for j, val in enumerate(potentials_imag)
            ]
            output_file.write("\n".join(enumerated_imag_data))
            output_file.write("\n$EndNodeData\n")

            # Compute magnitude and phase
            potentials_complex = np.array(potentials_real) + 1j * np.array(
                potentials_imag
            )
            potentials_magnitude = np.abs(potentials_complex)
            potentials_phase = np.angle(potentials_complex)

            # Append the magnitude part of the solutions in the specified format
            output_file.write(
                '$NodeData\n1\n"Magnitude Electric Potential"\n1\n0.0\n3\n'
            )
            output_file.write(f"{pattern_number}\n1\n{nodes_count}\n")
            enumerated_magnitude_data = [
                f"{j + 1}\t{val}" for j, val in enumerate(potentials_magnitude)
            ]
            output_file.write("\n".join(enumerated_magnitude_data))
            output_file.write("\n$EndNodeData\n")

            # Append the phase part of the solutions in the specified format
            output_file.write('$NodeData\n1\n"Phase Electric Potential"\n1\n0.0\n3\n')
            output_file.write(f"{pattern_number}\n1\n{nodes_count}\n")
            enumerated_phase_data = [
                f"{j + 1}\t{val}" for j, val in enumerate(potentials_phase)
            ]
            output_file.write("\n".join(enumerated_phase_data))
            output_file.write("\n$EndNodeData\n")
