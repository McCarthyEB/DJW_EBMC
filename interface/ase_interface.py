import numpy as np
from ase.io import read, write
from ase.visualize import view
from ase.build import add_adsorbate
from ase import atoms


# All angles seem to be at 60*
# Will take interface file as input: for simplicity might just use isolated Pd and O layers to start


# view(atoms)

# avector_length=(np.linalg.norm(cellvectors[0]))
# bvectorlength=(np.linalg.norm(cellvectors[1]))
# half_avector=(avector_length / 2)
# xy_centre=np.array([half_avector, half_avector])

def atom_set_finder(atoms, isymbol):
    cell_vectors = atoms.get_cell()
    print(cell_vectors)

    a_vector = cell_vectors[0]
    b_vector = cell_vectors[1]
    xy_plane_vector_1 = a_vector[:2]
    xy_plane_vector_2 = b_vector[:2]
    xy_plane_midpoint = 0.5 * (xy_plane_vector_1 + xy_plane_vector_2)

    #print("Midpoint of the XY plane is", xy_plane_midpoint)

    closest_distance = float(9999)
    closest_index = None

    for iatom in range(0, len(atoms)):
        #print(f"%d %s %10.6f %10.6f %10.6f" \
        #      % (iatom, atoms.symbols[iatom], atoms.positions[iatom][0], \
        #         atoms.positions[iatom][1], atoms.positions[iatom][2]))
        if isymbol in atoms.symbols[iatom]:
            distance_to_centre = np.linalg.norm(atoms.positions[iatom][:2] - xy_plane_midpoint)
            if distance_to_centre < closest_distance:
                closest_index = iatom
                closest_distance = distance_to_centre
    hip_hip_array=[[atoms.positions[closest_index][0], \
                    atoms.positions[closest_index][1], \
                    atoms.positions[closest_index][2]]]
    for secatom in range(0, len(atoms)):
        if isymbol in atoms.symbols[secatom]:
           dist=atoms.get_distance(closest_index, secatom, mic=True)
           if 2.5 < dist < 3.5:
             hip_hip_array.append([atoms.positions[secatom][0], \
                                   atoms.positions[secatom][1], \
                                   atoms.positions[secatom][2]])
                
    return hip_hip_array


atoms = read("Pd_O_isolated.cif")
Pd_index = atom_set_finder(atoms, "Pd")
print(Pd_index)
O_index = atom_set_finder(atoms, "O")
print(O_index)
print(len(O_index))








#


#
#
#
#
#
# view(pd_structure)
# transformation_matrix = np.dot(pd_structure.cell.T, np.linalg.inv(zno_structure.cell.T))
#
# pd_aligned_cell = pd_structure.cell @ transformation_matrix.T
# pd_aligned = pd_structure.copy()
# pd_aligned.set_cell(pd_aligned_cell, scale_atoms=True)
# view(pd_aligned)
##Isn't quite working: add in tolerance?
# tolerance = 0.001
# for pd_atom in pd_aligned:
#    for zno_atom in zno_structure:
#        displacement_vector = pd_atom.position - zno_atom.position
#        displacement_norm = np.linalg.norm(pd_aligned.get_distance(pd_atom.index, zno_atom.index, mic=True))
#        if displacement_norm < tolerance:
#            # Shift Pd atom to achieve the desired tolerance
#            fractional_coords = pd_aligned.get_scaled_positions()[pd_atom.index]
#            fractional_coords += zno_atom.position
#            pd_atom.position = pd_aligned.get_positions(fractional=True)[pd_atom.index]
#            break
#
#
# epitaxial_system = zno_structure.copy()
# for atom in pd_aligned:
#    epitaxial_system.append(atom)
#
# write('heterointerface.in', epitaxial_system, format='xyz')
# view(epitaxial_system)
#
##Not quite working as desired yet!
#
#
#
