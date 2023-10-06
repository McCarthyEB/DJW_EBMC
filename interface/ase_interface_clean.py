import numpy as np
from ase.io import read, write
from ase.visualize import view
from ase import Atom
from ase.build import make_supercell
from ase import Atoms


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
            # projection might be useful: corresponds to opposite side of right angled triangle formed by bond and direction vector
            projection = np.dot(vector, norm_direction_vector)
            if -perpendicular_range <= projection <= perpendicular_range:
                parallel_atoms.append(iatom)

    return parallel_atoms


def find_hexagonal_strip(atoms, seed_indices, seed_index, seed_symbol, perpendicular_range, tolerance):
    direction_vector = (atoms.positions[19] - atoms.positions[16])
    # Tidy this up later and generalize if it works!
    # Need to set it as the vector perpendicular to the desired strip for some reason
    parallel_atoms = find_parallel_atoms(atoms, seed_indices, direction_vector, seed_symbol, perpendicular_range,
                                         tolerance)

    return parallel_atoms


def calculate_relative_positions(Pd_atoms, O_atoms):
    relative_positions = []
    for pd_atom in Pd_atoms:
        for o_atom in O_atoms:
            relative_position = pd_atom.position - o_atom.position
            relative_positions.append(relative_position)
    return relative_positions


def cell_surface_finder(atoms):
    zmax = atoms.positions[0][2]
    zmin = zmax
    for iatom in range(0, len(atoms)):
        print(f"%d %s %10.6f %10.6f %10.6f charge: %10.2f" \
              % (iatom, atoms.symbols[iatom], atoms.positions[iatom][0], atoms.positions[iatom][1], atoms.positions[iatom][2],
                 atoms[iatom].charge))
        #
        # Identify top and bottom layers
        #
        if atoms.positions[iatom][2] > zmax:
            zmax = atoms.positions[iatom][2]
        #
        if atoms.positions[iatom][2] < zmin:
            zmin = atoms.positions[iatom][2]

    zdelta = 0.1
    top_layer = []
    bot_layer = []
    for iatom in range(0, len(atoms)):
        if ((atoms.positions[iatom][2] > (zmax - zdelta)) \
                & (atoms.positions[iatom][2] < (zmax + zdelta))):
            top_layer.append(iatom)

        if ((atoms.positions[iatom][2] < (zmin + zdelta)) \
                & (atoms.positions[iatom][2] > (zmin - zdelta))):
            bot_layer.append(iatom)
    #
    num_top_layer = len(top_layer)
    num_bot_layer = len(bot_layer)
    #
    print(f"Identified %d top layer atoms:" % len(top_layer))
    for iii in range(0, num_top_layer):
        iatom = top_layer[iii]
        print(f"%d %s %10.6f %10.6f %10.6f charge: %10.2f" \
              % (iatom, atoms.symbols[iatom], atoms.positions[iatom][0], atoms.positions[iatom][1],
                 atoms.positions[iatom][2],
                 atoms[iatom].charge))

    print(f"Identified %d bottom layer atoms:" % len(bot_layer))
    for iii in range(0, num_bot_layer):
        iatom = bot_layer[iii]
        print(f"%d %s %10.6f %10.6f %10.6f charge: %10.2f" \
              % (iatom, atoms.symbols[iatom], atoms.positions[iatom][0], atoms.positions[iatom][1], atoms.positions[iatom][2],
                 atoms[iatom].charge))


    return top_layer, bot_layer


def dist_list(atoms):
    atomdists = []
    atomvecs = []
    natoms = len(atoms)

    for iatom in range(0, natoms):
        # print(iatom)
        for sec_atom in range(iatom + 1, natoms):
            dist = atoms.get_distance(iatom, sec_atom, mic=True)
            vec = atoms.get_distance(iatom, sec_atom, mic=True, vector=True)
            print(dist)
            ilen = len(atomdists)
            isdifferent = True
            for idist in range(0, ilen):
                if abs(dist - atomdists[idist]) <= 0.1:
                    isdifferent = False
                    break
            if isdifferent:
                atomdists.append(dist)
                atomvecs.append(vec)

    # Adding a second loop heck to pick up the distance between atom repeat and its corresponding atom in the first cell
    for iatom in range(0, natoms):
        for sec_atom in range(natoms, 2 * natoms):
            if sec_atom % natoms == iatom and sec_atom < len(atoms):
                dist = atoms.get_distance(iatom, sec_atom,
                                          mic=False)  # mic=False to not apply minimum image convention!
                vec = atoms.get_distance(iatom, sec_atom, mic=False, vector=True)
                atomdists.append(dist)
                atomvecs.append(vec)

    #
    return atomdists, atomvecs


atoms = read("Pd_O_isolated.cif")
supercell = read("supercell_1_slabheight_14.cif")
supervac = read("slab_atoms_7_65839.cif")
# view(atoms)
view(supercell)

top_cell, bottom_cell = cell_surface_finder(supercell)
print("Top cell:", top_cell)
top_atoms = supercell[top_cell]
write('supercell_top.cif', top_atoms)
isolate_top = read('supercell_top.cif')

#
#view(isolate_top)
#

o_dist, o_vecs = dist_list(isolate_top)
# o_dist, o_vecs = dist_list(oxygen_strip_atoms)
print("Distances from central oxygen:")
for idist in range(0, len(o_dist)):
    print(o_dist[idist], "Å", "vec: ", o_vecs[idist])

Pd_hip_hip, Pd_indices = atom_set_finder(atoms, "Pd")
print("Pd_indices:", Pd_indices)
O_hip_hip, O_indices = atom_set_finder(atoms, "O")

seed_Pd_index = Pd_indices[0]  # This is the centre
print("seed_Pd_index", seed_Pd_index)

strip_atoms_indices = find_hexagonal_strip(atoms, Pd_indices, seed_Pd_index, "Pd", perpendicular_range=3.0,
                                           tolerance=0.1)
combined_indices = Pd_indices + strip_atoms_indices
pd_strip_atoms = atoms[combined_indices]
write("Pd_strip_atoms.cif", pd_strip_atoms)
view_pd_strip = read("Pd_strip_atoms.cif")
#view(view_pd_strip)

pd_dist, pd_vecs = dist_list(pd_strip_atoms)
print("Distances from central palladium:")
for dist in pd_dist:
    print(dist, "Å")

print("The strip of overlapping Pd hexagons contains %d atoms." % len(strip_atoms_indices))

seed_O_index = O_indices[0]
oxygen_strip_atoms_indices = find_hexagonal_strip(atoms, O_indices, seed_O_index, "O", perpendicular_range=4.0,
                                                  tolerance=1)
oxygen_combined_indices = O_indices + oxygen_strip_atoms_indices
oxygen_strip_atoms = atoms[oxygen_combined_indices]
write("oxygen_strip_atoms.cif", oxygen_strip_atoms)
view_oxygen_strip = read("oxygen_strip_atoms.cif")
#view(view_oxygen_strip)

print("The strip of overlapping O atoms contains %d atoms." % len(oxygen_strip_atoms))


######
def map_pd_to_o(pd_atoms, o_atoms):
    pd_to_o_mapping = {}  # Create a dictionary to store the mapping
    pd_to_o_distances = {}  # Create a dictionary to store the distances

    for pd_atom in pd_atoms:
        min_distance = None
        closest_oxygen = None

        for o_atom in o_atoms:
            distance = np.linalg.norm(pd_atom.position - o_atom.position)

            if min_distance is None or distance < min_distance:
                min_distance = distance
                closest_oxygen = o_atom

        pd_to_o_mapping[pd_atom.index] = closest_oxygen.index  # Store the mapping
        pd_to_o_distances[pd_atom.index] = min_distance  # Store the distance

    return pd_to_o_mapping, pd_to_o_distances


# Map Pd atoms in 'strip_atoms' to their closest O atoms in 'oxygen_strip_atoms' with distances
# Changing oxygen_strip_atoms to isolate_top
pd_to_o_mapping, pd_to_o_distances = map_pd_to_o(pd_strip_atoms, isolate_top)

# Print the mapping and distances
for pd_index, o_index in pd_to_o_mapping.items():
    distance = pd_to_o_distances[pd_index]
    print(f"Pd {pd_index} -> O {o_index} (Distance: {distance} Å)")

# Let a and b be the bond-containing vector that intersects the cell at points (a,0), (0,b)


#
# Look for best match for each Pd distance with the oxygen lattice
#
good_matches = []
for idist in range(0, len(pd_dist)):
    min = 9999.0
    for jdist in range(0, len(o_dist)):
        ratio = o_dist[jdist] / pd_dist[idist]
        #
        # Work out remainder
        #
        remain = ratio - float(int(ratio))
        print("Ratio for Pd dist %d with oxygen dust %d: %10.6f (rem: %10.6f)" % (idist, jdist, ratio, remain))
        #
        if remain < 0.10:
            print("Found match within 10%")
            good_matches.append([idist, jdist])
    print("")
#
print("Good matches: ", good_matches)
#
for imatch in range(0, len(good_matches)):
    #
    pd_index = good_matches[imatch][0]
    o_index = good_matches[imatch][1]
    pd_dist_rep = pd_dist[pd_index]
    print(pd_index, o_index)
    print("Looking at match %d, pd_dist: %10.6f o_dist: %10.6f" % (imatch, pd_dist[pd_index], o_dist[o_index]))
    print("Pd_vec: ", pd_vecs[pd_index])
    print("O_vec : ", o_vecs[o_index])
    #
    # Make a unit vector in the direction of interest
    #
    o_best_vec = o_vecs[o_index].copy()
    #
    uni_rep = o_best_vec / np.linalg.norm(o_best_vec)
    #
    # Find fractional co-ordinates
    latt = supervac.get_cell()
    recip_latt = supervac.cell.reciprocal()
    # recip_latt=oxygen_strip_atoms.get_reciprocal_cell()
    print("Lattice           :", latt)
    print("Reciprocal Lattice:", recip_latt)
    #
    frac = []
    for iii in range(0, 3):
        frac.append(np.dot(recip_latt[iii], o_best_vec))
    #
    print("Fractional version of inter-neighbour vector:")
    print(frac)
    #
    #
    # The a component is longer than the b so use that
    if (frac[0] > frac[1]):
        num_in_cell = int(1.0 / frac[0])
        print("Longest component is in a-vec direction")
    # The b component is longer than the a so use that
    else:
        num_in_cell = int(1.0 / frac[1])
        print("Longest component is in b-vec direction")
    #
    dist_in_cell = num_in_cell * o_dist[o_index]
    vec_in_cell = dist_in_cell * uni_rep
    #
    print("vec_in_cell: ", vec_in_cell)
    #
    # Work out supercell required to get vec in cell to repeat properly.
    #
    frac = []
    for iii in range(0, 3):
        frac.append(np.dot(recip_latt[iii], vec_in_cell))
    #
    print("frac version of vec_in_cell: ", frac)
    #

    num_for_super = [int(abs(1.0 / frac[0])), int(abs(1.0 / frac[1]))]
    #
    # make a supercell big enough for the repeat:
    #
    mat = [[float(num_for_super[1]), 0, 0], [0, float(num_for_super[0]), 0], [0, 0, 1]]
    print("Matrix for supercell:")
    print(mat)
    #No longer need this: determinant is 0!
    if np.linalg.det(mat) == 0:
        super = supervac
    else:
        super = make_supercell(supervac, mat)
    #
    new_system = super.copy()
    #
    # Find fractional co-ordinates in the supercell
    latt = new_system.get_cell()
    recip_latt = new_system.cell.reciprocal()
    # recip_latt=new_system.get_reciprocal_cell()
    print("Super Lattice           :", latt)
    print("Super Reciprocal Lattice:", recip_latt)
    #
    frac = []
    for iii in range(0, 3):
        frac.append(np.dot(recip_latt[iii], o_best_vec))
    #
    print("Fractional version of inter-neighbour vector in supercell:")
    print(frac)
    #
    #
    # The a component is longer than the b so use that
    if (frac[0] > frac[1]):
        num_in_cell = int(1.0 / frac[0])
        print("Longest component is in a-vec direction")
    # The b component is longer than the a so use that
    else:
        num_in_cell = int(1.0 / frac[1])
        print("Longest component is in b-vec direction")
    #
    dist_in_cell = num_in_cell * o_dist[o_index]
    # Issue with this line: it's taking dist_in_cell as 22.622 instead of 28.278
    # Taking it as the diagonal across the cell!
    #
    # Need the co-ords of one of the oxygens to set the origin, offset a little above the oxygens.
    #
    origin = isolate_top.positions[0].copy()
    origin[2] = origin[2] + 1.5
    #
    num_pd_in_cell = int(dist_in_cell / pd_dist[pd_index])
    # HERE! I want to use o_dist because pd_dist is the old cell!
    # It's taking the diagonal here in the original cell as the 'longest distance'.
    #
    pd_line_dist = num_pd_in_cell * pd_dist[pd_index]
    #
    stretch_factor = dist_in_cell / pd_line_dist
    #
    print("longest dist in cell = %10.6f can fit %d Pd covering dist: %10.6f stretch factor %10.6f" \
          % (dist_in_cell, num_pd_in_cell, pd_line_dist, stretch_factor))

    for iii in range(0, num_pd_in_cell + 2):  # This works. But why???
        new_coords = origin + iii * stretch_factor * pd_dist[pd_index] * uni_rep
        new_atom = Atom("Pd", new_coords)
        new_system.append(new_atom)
        new_system.wrap()

    write("new_system_%s.cif" % imatch, new_system)
    read_system = read("new_system_%s.cif" % imatch)
    view(read_system)

###
