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
