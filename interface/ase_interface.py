import numpy as np
from ase.io import read, write
from ase.visualize import view

def atom_set_finder(atoms, isymbol):
    cell_vectors = atoms.get_cell()

    a_vector = cell_vectors[0]
    b_vector = cell_vectors[1]
    xy_plane_vector_1 = a_vector[:2]
    xy_plane_vector_2 = b_vector[:2]
    xy_plane_midpoint = 0.5 * (xy_plane_vector_1 + xy_plane_vector_2)

    closest_distance = float(9999)
    closest_index = None

    for iatom in range(0, len(atoms)):
        if isymbol in atoms.symbols[iatom]:
            distance_to_centre = np.linalg.norm(atoms.positions[iatom][:2] - xy_plane_midpoint)
            if distance_to_centre < closest_distance:
                closest_index = iatom
                closest_distance = distance_to_centre

    hip_hip_array = [atoms.positions[closest_index]]
    atom_array = [closest_index]

    for secatom in range(0, len(atoms)):
        if isymbol in atoms.symbols[secatom] and secatom != closest_index:  # Exclude the central Pd atom
            dist = atoms.get_distance(closest_index, secatom, mic=True)
            if 2.5 < dist < 3.5:
                hip_hip_array.append(atoms.positions[secatom])
                atom_array.append(secatom)

    return hip_hip_array, atom_array

def find_parallel_atoms(atoms, seed_indices, direction_vector, isymbol, perpendicular_range, tolerance=0.1):
    parallel_atoms = []
    seed_positions = [atoms.positions[i] for i in seed_indices]

    norm_direction_vector = direction_vector / np.linalg.norm(direction_vector)

    for iatom in range(len(atoms)):
        if iatom in seed_indices:
            continue

        if isymbol in atoms.symbols[iatom]:
            position = atoms.positions[iatom]
            vector = position - seed_positions[0]
            #projection might be useful: corresponds to opposite side of right angled triangle formed by bond and direction vector
            projection = np.dot(vector, norm_direction_vector)
            if -perpendicular_range <= projection <= perpendicular_range:
                parallel_atoms.append(iatom)

    return parallel_atoms

def find_hexagonal_strip(atoms, seed_indices, seed_index, perpendicular_range, tolerance):
    direction_vector = (atoms.positions[19] - atoms.positions[16])
    # Tidy this up later and generalize if it works!
    # Need to set it as the vector perpendicular to the desired strip for some reason
    parallel_atoms = find_parallel_atoms(atoms, seed_indices, direction_vector, "Pd", perpendicular_range, tolerance)

    return parallel_atoms

atoms = read("Pd_O_isolated.cif")
view(atoms)
Pd_hip_hip, Pd_indices = atom_set_finder(atoms, "Pd")
print("Pd_indices:", Pd_indices)
O_hip_hip, O_indices = atom_set_finder(atoms, "O")

seed_Pd_index = Pd_indices[0]  #This is the centre
print("seed_Pd_index", seed_Pd_index)

strip_atoms_indices = find_hexagonal_strip(atoms, Pd_indices, seed_Pd_index, perpendicular_range=3.0, tolerance=0.1)
combined_indices = Pd_indices + strip_atoms_indices
strip_atoms = atoms[combined_indices]
write("strip_atoms.cif", strip_atoms)
view_strip = read("strip_atoms.cif")
view(view_strip)

print("The strip of overlapping Pd hexagons contains %d atoms." % len(strip_atoms))
