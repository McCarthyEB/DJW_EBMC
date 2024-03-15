#
# development makefile
#
CC = gcc
#
#
# jens test comment
 BIN = ~/bin
#  BIN = .
#
#
#CFLAGS =
#CFLAGS =   -fullwarn  -O1 -mips4
#CFLAGS =   -sopt -fullwarn -g
#CFLAGS =   -O2 -fullwarn
#CFLAGS =   -O2 -mips4
#CFLAGS =   -O2 -Wall
#CFLAGS =  -O2 -g -fno-stack-protector
#CFLAGS =    -Wno-unused-result -g 
CFLAGS =   -O2 

#CFLAGS_5000 = -r5000  -sopt -fullwarn  -O2

#Flags known to work are:
#Irix/R4000/R8000 -O2 -xansi -fullwarn
#Irix6.3/R10000 (o2) its rubbish! Use O1 at the most (and compare to others) and mips4
#AIX/RS6000 -O3
#Linux -O2

#
COBJS = main.o energy.o potent.o choose_fragment.o connect_new_fragment.o find_hydrogens.o \
       standard_bond.o generate_neighbours.o compare_strings.o \
       join_vectors.o size_vector.o print_neighbours.o print_neighbours_heavies.o print_molecule.o \
       move_molecule.o vec_cross.o vec_dot.o rotate.o unit_vector.o \
       animate.o print_biosym_header.o output_routines.o join_atoms.o \
       centre_of_mass.o atomic_mass_list.o bump_check.o \
       read_line.o get_int.o next_space.o rock_molecule.o shake_molecule.o \
       relabel.o rot_section.o print_mdf.o \
       minimize_molecule.o minimize_inpore.o \
       disco_minimize_molecule.o c2disco_minimize_molecule.o c2disco_minimize_inpore.o \
       disco_minimize_inpore.o mopac_minimize_molecule.o mopac_minimize_inpore.o \
       get_c2disco_minimized.o print_strategy.o anal_disco.o anal_c2disco.o \
       anal_mopac.o regroup.o  locate_string.o cart_latt_vecs.o  \
       min_image.o fract_to_cart.o cart_to_fract.o int_to_char.o \
       timer.o real_random.o rotate_with_flags.o find_chunk.o \
       set_twist_axis_flags.o do_the_twist.o test_the_twist.o get_doub.o calculate_box.o print_box.o \
       pick_frm_wgt_list.o box_test.o cost_function.o do_symmetry.o \
       matrix_vec_mult.o banner.o bump_check_with_flags.o \
       fix_potentials.o fix_potentials_cff91.o calculate_charge.o \
       check_element.o print_dist_mat.o \
       update_partition_list.o join_frag_to_temp.o align_bonds.o set_dihedral.o \
       shuffle_neighbours.o shuffle_links.o shuffle_hyd_list.o \
       add_types.o copy_types.o  copy_strings.o test_concentrations.o \
       bonds_possible.o find_hyd_neighbours.o forbid_bond.o print_types.o \
       setup_defaults.o initialise_variables.o enforce_defaults.o \
       assemble_molmol_list.o \
       make_a_template.o initialise_animation.o open_file.o int_from_string.o \
       string_from_int.o next_none_space.o get_nonbond_info.o \
       get_nonbond_params.o read_potentials.o non_bond_r_eps_sixth.o \
       get_stretch_params.o assign_stretch.o calculate_intra_energy.o \
       get_equivalences.o neighbour_order.o ring_maker.o shuffle_atom_list.o \
       atom_has_hyds.o atom_separation_squared.o sort_out_molecule.o \
       get_seed_molecule.o get_pre_anneal_template.o print_seed_info.o pbc_interactions.o \
       print_statistics.o print_molecule_xmol.o analyse_output.o \
       strip_spaces.o index_sort.o delete_file.o print_peek_file.o \
       find_zebedde_file.o locate_string_in_string.o \
       ewald_sum.o ewald_real_part.o ewald_recip_part.o own_erf.o erf_mac.o \
       erf_asymptotic.o isol_coul.o internal_minimize_inpore.o print_product.o \
       matrix_x_ket.o ket_bra_to_matrix.o bra_x_ket.o bra_mat_ket.o bfgs_update.o \
       apply_delta_rigid.o rigid_body_gradients.o matrix_x_vector.o \
       non_bond_a_b_geometric.o non_bond_a_b_none.o gen_recip_klist.o \
       find_line_atoms.o non_bond_r_eps_arithmetic.o get_hbond_params.o \
       order_double_pot.o assign_hbond.o hbond_amber.o get_bond_incs.o \
       index_double_pot.o assign_bond_inc.o print_mdf_subset.o\
       neutralise_slab.o write_pdb.o write_xyz.o get_angle_params.o get_torsion_params.o \
       assemble_angle_list.o assemble_torsion_list.o look_up_equiv_pots.o \
       assemble_vdw_list.o angle_bend_energy.o torsion_energy.o assemble_intra_lists.o \
       deactivate_groups.o anal_gulp.o gulp_minimise_molecule.o gulp_minimize_inpore.o \
       dock_orientate.o gulp_inpore_justinput.o \
       gulp_inpore_annealed.o find_mol.o generate_input.o car_molecule_output.o condor_submit.o \
       os_query.o cluster.o disco_anneal_inpore.o shift_coords.o count_vdw_pbc.o \
       gather_molecule.o reader.o demarcate_molecule.o test_pore_symmetry.o \
       do_symmetry_all.o check_allowed_torsions.o copy_molecule_data.o

#
# Removed routines Dave Willock Nov 2005
#
#      license.o
#
#
LFLAGS = -lm
#


main.o             : main.c make_a_template.c  real_random.c  \
                     calculate_box.c setup_defaults.c enforce_defaults.c \
                     initialise_animation.c initialise_variables.c \
		     print_biosym_header.c \
                     generate_neighbours.c neutralise_slab.c \
                     cart_latt_vecs.c cost_function.c do_symmetry.c add_types.c \
                     find_zebedde_file.c order_double_pot.c \
                     maxima.h global_values.h constants.h structures.h header.h \
                     data.h ewald.h own_maths.h 
output_routines.o  : output_routines.c print_types.c structures.h  header.h 
reader.o           : reader.c locate_string.c structures.h  \
					 header.h reader.h data.h 
energy.o           : energy.c output_routines.c ewald_sum.c isol_coul.c min_image.c \
                     non_bond_r_eps_sixth.c non_bond_a_b_geometric.c non_bond_r_eps_arithmetic.c \
                     maxima.h ewald.h global_values.h constants.h \
                     structures.h  header.h data.h
connect_new_fragment.o : connect_new_fragment.c
choose_fragment.o  : choose_fragment.c add_types.c copy_types.c pick_frm_wgt_list.c test_concentrations.c \
                     maxima.h global_values.h structures.h header.h data.h
potent.o           :  potent.c structures.h  header.h data.h
find_hydrogens.o   : find_hydrogens.c structures.h
compare_strings.o  : compare_strings.c  data.h
generate_neighbours.o: generate_neighbours.c atom_separation_squared.c \
                       standard_bond.c header.h structures.h maxima.h
standard_bond.o    : standard_bond.c data.h
join_vectors.o     : join_vectors.c   write_pdb.o
size_vector.o      : size_vector.c
print_neighbours.o : print_neighbours.c structures.h
print_neighbours_heavies.o : print_neighbours_heavies.c structures.h
print_mdf.o 	   : print_neighbours.c structures.h
print_molecule.o   : print_molecule.c structures.h
animate.o          : animate.c structures.h  header.h
move_molecule.o    : move_molecule.c structures.h
vec_cross.o        : vec_cross.c
vec_dot.o          : vec_dot.c
rotate.o           : rotate.c structures.h
unit_vector.o      : unit_vector.c
print_biosym_header.o : print_biosym_header.c 
join_atoms.o       : join_atoms.c structures.h
centre_of_mass.o   : centre_of_mass.c atomic_mass_list.c structures.h
atomic_mass_list.o : atomic_mass_list.c data.h
bump_check.o       : bump_check.c structures.h 
read_line.o        : read_line.c
get_int.o          : get_int.c
next_space.o       : next_space.c
shake_molecule.o   : shake_molecule.c  real_random.c cost_function.c structures.h header.h
rock_molecule.o    : rock_molecule.c real_random.c cost_function.c  structures.h header.h
int_to_char.o      : int_to_char.c data.h
relabel.o          : relabel.c int_to_char.c data.h structures.h
rot_section.o      : rot_section.c real_random.c rotate.c cost_function.c unit_vector.c \
                     join_atoms.c structures.h header.h
fix_potentials.o   : fix_potentials_cff91.c fix_potentials.c \
                     structures.h header.h
fix_potentials_cff91.o   : fix_potentials_cff91.c structures.h header.h
calculate_charge.o : calculate_charge.c structures.h
delete_file.o      : delete_file.c
minimize_molecule.o: minimize_molecule.c disco_minimize_molecule.c c2disco_minimize_molecule.c \
                     mopac_minimize_molecule.c structures.h header.h data.h
minimize_inpore.o  : minimize_inpore.c structures.h header.h data.h
mopac_minimize_molecule.o: mopac_minimize_molecule.c calculate_charge.c\
						delete_file.c	structures.h  header.h data.h
c2disco_minimize_inpore.o: c2disco_minimize_inpore.c delete_file.c anal_c2disco.c \
                             get_c2disco_minimized.c structures.h header.h
c2disco_minimize_molecule.o: c2disco_minimize_molecule.c delete_file.c anal_c2disco.c \
                             get_c2disco_minimized.c structures.h header.h
get_c2disco_minimized.o: get_c2disco_minimized.c maxima.h header.h
disco_minimize_molecule.o: disco_minimize_molecule.c delete_file.c \
							structures.h header.h
disco_minimize_inpore.o  : disco_minimize_inpore.c delete_file.c \
							structures.h header.h
mopac_minimize_inpore.o  : mopac_minimize_inpore.c calculate_charge.c\
							delete_file.c structures.h header.h
print_strategy.o   : print_strategy.c 
anal_disco.o       : anal_disco.c structures.h maxima.h global_values.h structures.h header.h
anal_c2disco.o     : anal_c2disco.c structures.h maxima.h global_values.h structures.h header.h
anal_mopac.o       : anal_mopac.c structures.h header.h
regroup.o          : regroup.c structures.h
locate_string.o    : locate_string.c
cart_latt_vecs.o   : cart_latt_vecs.c constants.h ewald.h own_maths.h
min_image.o        : min_image.c fract_to_cart.c cart_to_fract.c  header.h
fract_to_cart.o    : fract_to_cart.c
cart_to_fract.o    : cart_to_fract.c
timer.o            : timer.c
real_random.o      : real_random.c
set_twist_axis_flags.o : set_twist_axis_flags.c real_random.c unit_vector.c join_atoms.c \
                         find_chunk.c gather_molecule.c maxima.h global_values.h structures.h \
                         constants.h data.h header.h
do_the_twist.o     : do_the_twist.c real_random.c rotate_with_flags.c do_symmetry.c \
                     maxima.h global_values.h structures.h \
                     constants.h data.h header.h
test_the_twist.o   : test_the_twist.c cost_function.c \
                     maxima.h global_values.h structures.h \
                     constants.h data.h header.h
find_chunk.o       : find_chunk.c structures.h
rotate_with_flags.o: rotate_with_flags.c move_molecule.c structures.h
get_doub.o         : get_doub.c get_int.c
calculate_box.o    : calculate_box.c structures.h header.h
print_box.o        : print_box.c print_biosym_header.c structures.h header.h
pick_frm_wgt_list.o: pick_frm_wgt_list.c real_random.c 
box_test.o         : box_test.c structures.h header.h
cost_function.o    : cost_function.c energy.c bump_check.c box_test.c \
                     structures.h header.h constants.h data.h
do_symmetry.o      : do_symmetry.c fract_to_cart.c cart_to_fract.c matrix_vec_mult.c \
                     min_image.c centre_of_mass.c move_molecule.c structures.h header.h
matrix_vec_mult.o  : matrix_vec_mult.c
banner.o           : banner.c
check_element.o    : check_element.c data.h
bump_check_with_flags.o: bump_check_with_flags.c  structures.h header.h
update_partition_list.o : update_partition_list.c structures.h
join_frag_to_temp.o : join_frag_to_temp.c pick_frm_wgt_list.c shuffle_neighbours.c \
                      forbid_bond.c align_bonds.c shuffle_atom_list.c shuffle_hyd_list.c \
                      maxima.h global_values.h structures.h header.h data.h 
align_bonds.o : align_bonds.c set_dihedral.c move_molecule.c join_atoms.c vec_cross.c \
                vec_dot.c size_vector.c standard_bond.c unit_vector.c \
                rotate.c structures.h header.h data.h
set_dihedral.o       : set_dihedral.c rotate.c vec_dot.c vec_cross.c unit_vector.c \
                       join_atoms.c maxima.h constants.h global_values.h structures.h data.h
shuffle_neighbours.o : shuffle_neighbours.c structures.h
shuffle_links.o      : shuffle_links.c structures.h
shuffle_hyd_list.o   : shuffle_hyd_list.c structures.h
add_types.o          : add_types.c compare_strings.c structures.h header.h
copy_strings.o       : copy_strings.c 
copy_types.o         : copy_types.c structures.h
test_concentrations.o : test_concentrations.c compare_strings.c structures.h header.h
bonds_possible.o     : bonds_possible.c compare_strings.c copy_strings.c structures.h header.h
find_hyd_neighbours.o : find_hyd_neighbours.c compare_strings.c copy_strings.c structures.h \
                        header.h
forbid_bond.o        : forbid_bond.c compare_strings.c structures.h header.h
print_types.o        : print_types.c structures.h
setup_defaults.o     : setup_defaults.c real_random.c structures.h header.h \
                       constants.h
initialise_variables.o : initialise_variables.c structures.h header.h \
                       constants.h
enforce_defaults.o   : enforce_defaults.c structures.h header.h \
                       constants.h
print_peek_file.o  :      print_peek_file.c structures.h print_molecule.c \
						  print_biosym_header.c 
make_a_template.o    : make_a_template.c join_vectors.c move_molecule.c \
                       ring_maker.c rot_section.c do_the_twist.c test_the_twist.c \
                       centre_of_mass.c main.c real_random.c set_twist_axis_flags.c \
                       pick_frm_wgt_list.c box_test.c calculate_box.c \
                       print_box.c cost_function.c do_symmetry.c add_types.c \
                       calculate_intra_energy.c write_pdb.c write_xyz.c print_dist_mat.c \
                       assemble_intra_lists.c structures.h header.h \
                       constants.h output_routines.c print_peek_file.c
open_file.o          : open_file.c 
write_pdb.o          : write_pdb.c
write_xyz.o          : write_xyz.c
initialise_animation.o : initialise_animation.c 
int_from_string.o    : int_from_string.c
string_from_int.o    : string_from_int.c
next_none_space.o    : next_none_space.c
get_nonbond_info.o   : get_nonbond_info.c locate_string.c next_none_space.c \
                       next_space.c read_line.c string_from_int.o global_values.h data.h
get_nonbond_params.o : get_nonbond_params.c locate_string.c string_from_int.c \
                       next_none_space.c next_space.c read_line.c \
                       global_values.h data.h
read_potentials.o    : read_potentials.c open_file.c read_line.c \
                       locate_string.c string_from_int.c get_nonbond_info.c \
                       get_nonbond_params.c get_stretch_params.c get_angle_params.c\
                       get_torsion_params.c global_values.h data.h
non_bond_r_eps_sixth.o : non_bond_r_eps_sixth.c min_image.c maxima.h global_values.h \
                         structures.h data.h constants.h
non_bond_a_b_geometric.o : non_bond_a_b_geometric.c min_image.c maxima.h global_values.h \
                           structures.h data.h constants.h
non_bond_a_b_none.o : non_bond_a_b_none.c min_image.c maxima.h global_values.h \
                           structures.h data.h constants.h
non_bond_r_eps_arithmetic.o : non_bond_r_eps_arithmetic.c min_image.c maxima.h global_values.h \
                              structures.h data.h constants.h
get_stretch_params.o : get_stretch_params.c locate_string.c \
                               string_from_int.c \
                               next_none_space.c next_space.c read_line.c \
                               global_values.h data.h
get_angle_params.o : get_angle_params.c locate_string.c \
                               string_from_int.c \
                               next_none_space.c next_space.c read_line.c \
                               global_values.h data.h
get_torsion_params.o : get_torsion_params.c locate_string.c \
                               string_from_int.c \
                               next_none_space.c next_space.c read_line.c \
                               global_values.h data.h
calculate_intra_energy.o : calculate_intra_energy.c min_image.c maxima.h global_values.h \
                           structures.h data.h constants.h header.h
assign_stretch.o         : assign_stretch.c maxima.h global_values.h structures.h \
                           data.h header.h
get_equivalences.o       : get_equivalences.c locate_string.c string_from_int.c \
                           next_none_space.c next_space.c read_line.c maxima.h \
                           structures.h global_values.h data.h header.h
neighbour_order.o        : neighbour_order.c maxima.h global_values.h structures.h \
                           header.h
shuffle_atom_list.o      : shuffle_atom_list.c maxima.h global_values.h structures.h header.h 
atom_has_hyds.o          : atom_has_hyds.c structures.h global_values.h 
atom_separation_squared.o : atom_separation_squared.c min_image.c structures.h header.h
ring_maker.o             : ring_maker.c shuffle_neighbours.c shuffle_links.c \
                           shuffle_hyd_list.c shuffle_atom_list.c atom_has_hyds.c \
                           neighbour_order.c forbid_bond.c atom_separation_squared.c \
                           maxima.h global_values.h structures.h data.h header.h
sort_out_molecule.o:       sort_out_molecule.c generate_neighbours.c centre_of_mass.c \
                           assign_stretch.c assign_bond_inc.c assign_hbond.c  \
                           potent.c find_hydrogens.c assign_hbond.c \
                           structures.h data.h header.h maxima.h
get_seed_molecule.o:      get_seed_molecule.c structures.h data.h header.h
get_pre_anneal_template.o: get_pre_anneal_template.c structures.h data.h header.h
Gprint_seed_info.o  :      print_seed_info.c  structures.h data.h header.h
pbc_interactions.o :      pbc_interactions.c cart_to_fract.c fract_to_cart.c min_image.c \
                          header.h maxima.h global_values.h structures.h  
print_statistics.o :      print_statistics.c structures.h data.h
print_molecule_xmol.o   : print_molecule_xmol.c maxima.h global_values.h structures.h
analyse_output.o    :      analyse_output.c generate_neighbours.c print_mdf.c \
                          print_neighbours.c print_molecule.c output_routines.c\
			  relabel.c print_biosym_header.c\
                          structures.h header.h
strip_spaces.o     :      strip_spaces.c
index_sort.o       :      index_sort.c 
find_zebedde_file.o:      find_zebedde_file.c data.h
#license.o          :      license.c
locate_string_in_string.o : locate_string_in_string.c global_values.h
assemble_intra_lists.o : assemble_intra_lists.c assemble_angle_list.c assemble_torsion_list.c\
                         assemble_vdw_list.c data.h global_values.h structures.h maxima.h header.h
assemble_angle_list.o : assemble_angle_list.c data.h global_values.h structures.h maxima.h header.h
assemble_torsion_list.o : assemble_torsion_list.c look_up_equiv_pots.c data.h\
                          global_values.h structures.h maxima.h header.h
assemble_vdw_list.o : assemble_vdw_list.c look_up_equiv_pots.c min_image.c data.h\
                          global_values.h structures.h maxima.h header.h
assemble_molmol_list.o : assemble_molmol_list.c  \
                          global_values.h structures.h maxima.h header.h
look_up_equiv_pots.o : look_up_equiv_pots.c data.h global_values.h structures.h maxima.h header.h
angle_bend_energy.o  : angle_bend_energy.c data.h global_values.h structures.h maxima.h header.h
torsion_energy.o  : torsion_energy.c data.h global_values.h structures.h maxima.h header.h
find_mol.o : find_mol.c
count_vdw_pbc.o : count_vdw_pbc.c 
demarcate_molecule.o : demarcate_molecule.c maxima.h global_values.h structures.h
test_pore_symmetry.o : test_pore_symmetry.c do_symmetry.c min_image.c maxima.h global_values.h structures.h
do_symmetry_all.o : do_symmetry_all.c do_symmetry.c maxima.h global_values.h structures.h
check_allowed_torsions.o : check_allowed_torsions.c find_chunk.c structures.h


#Ewald summation bits
ewald_sum.o        :      ewald_sum.c ewald_real_part.c ewald_recip_part.c \
                          structures.h ewald.h
ewald_real_part.o  :      ewald_real_part.c min_image.c own_erf.c structures.h \
			  ewald.h maxima.h \
                          own_maths.h header.h
ewald_recip_part.o :      ewald_recip_part.c min_image.c structures.h ewald.h \
                          header.h
own_erf.o          :      own_erf.c erf_mac.c erf_asymptotic.c own_maths.h
erf_mac.o          :      erf_mac.c own_maths.h
erf_asymptotic.o   :      erf_asymptotic.c own_maths.h
isol_coul.o             : isol_coul.c structures.h ewald.h maxima.h own_maths.h
internal_minimize_inpore.o : internal_minimize_inpore.c print_molecule.c output_routines.c \
                             centre_of_mass.c bfgs_update.c bra_x_ket.c\
                             maxima.h global_values.h structures.h data.h
gen_recip_klist.o : gen_recip_klist.c structures.h maxima.h own_maths.h ewald.h
print_product.o  : print_product.c
bfgs_update.o    : bfgs_update.c matrix_x_ket.c ket_bra_to_matrix.c bra_x_ket.c \
                   bra_mat_ket.c
matrix_x_ket.o   : matrix_x_ket.c
ket_bra_to_matrix.o : ket_bra_to_matrix.c
bra_x_ket.o      : bra_x_ket.c
bra_mat_ket.o    : bra_mat_ket.c
apply_delta_rigid.o :  apply_delta_rigid.c move_molecule.c rotate.c structures.h
rigid_body_gradients.o : rigid_body_gradients.c structures.h
matrix_x_vector.o      : matrix_x_vector.c
get_hbond_params.o     : get_hbond_params.c
order_double_pot.o     : order_double_pot.c structures.h maxima.h global_values.h \
                         structures.h data.h header.h
assign_hbond.o         : assign_hbond.c maxima.h global_values.h structures.h data.h \
                         header.h 
assign_bond_inc.o      : assign_bond_inc.c maxima.h global_values.h structures.h data.h \
                         header.h
hbond_amber.o          : hbond_amber.c min_image.c pbc_interactions.c hbond_amber.c \
                         maxima.h global_values.h structures.h data.h constants.h \
                         header.h
get_bond_incs.o        : get_bond_incs.c get_int.c get_doub.c read_line.c next_space.c \
                         next_none_space.c string_from_int.c locate_string_in_string.c \
                         data.h global_values.h structures.h maxima.h header.h
index_double_pot.o     : index_double_pot.c global_values.h
print_mdf_subset.o     : print_mdf_subset.c
print_dist_mat.o       : print_dist_mat.c min_image.c 
deactivate_groups.o    : deactivate_groups.c
anal_gulp.o            : anal_gulp.c
gulp_minImise_molecule.o : gulp_minimise_molecule.c
gulp_minimize_inpore.o : gulp_minimize_inpore.c
gulp_inpore_justinput.o : gulp_inpore_justinput.c
gulp_inpore_annealed.o :  gulp_inpore_annealed.c
dock_orientate.o      :  dock_orientate.c 
generate_input.o       :generate_input.c
car_molecule_output.o  :car_molecule_output.c
condor_submit.o        :condor_submit.c
os_query.o	       :os_query.c
cluster.o : 	cluster.c
disco_anneal_inpore.o: disco_anneal_inpore.c
shift_coords.o :shift_coords.c cart_to_fract.c fract_to_cart.c header.h
gather_molecule.o : gather_molecule.c
copy_molecule_data.o : copy_molecule_data.c structures.h

#Restraint bits

find_line_atoms.o : find_line_atoms.c structures.h header.h global_values.h maxima.h

#Slab bits
neutralise_slab.o : neutralise_slab.c structures.h header.h


zebedde_2017: $(COBJS) makefile
		$(CC) $(CFLAGS)  -o $(BIN)/zebedde_2017 $(COBJS) $(LFLAGS)
		echo "Made 2017 version"

zebedde_mc: $(COBJS) makefile
		$(CC) $(CFLAGS)  -o $(BIN)/zebedde_mc    $(COBJS) $(LFLAGS)
		echo "Made MC  version"

zebedde_trial: $(COBJS) makefile
		$(CC) $(CFLAGS)  -o $(BIN)/zebedde_trial $(COBJS) $(LFLAGS)
		echo "Made trial version"

zebedde_test: $(COBJS) makefile
		$(CC) $(CFLAGS)  -o $(BIN)/zebedde_test $(COBJS) $(LFLAGS)
		echo "Made test version"

zebedde_5000: $(COBJS) makefile
	       $(CC) $(CFLAGS_5000)  -o $(BIN)/zebedde_5000 $(COBJS) $(LFLAGS)
		echo "Made R5000 version"

zebedde_lancaster: $(COBJS) $(COBJS_LINUX) makefile
		$(CC) $(CFLAGS)  -o $(BIN)/zebedde_lancaster $(COBJS) $(LFLAGS)
		echo "Made lancaster builder version"

zebedde_poly: $(COBJS) $(COBJS_LINUX) makefile
		$(CC) $(CFLAGS)  -o $(BIN)/zebedde_poly $(COBJS) $(LFLAGS)
		echo "Made polymer builder version"

zebedde_kent: $(COBJS) $(COBJS_LINUX) makefile
		$(CC) $(CFLAGS)  -o $(BIN)/zebedde_kent $(COBJS) $(LFLAGS)
		echo "Made polymer builder version for kent"

zebedde: $(COBJS) $(COBJS_LINUX) makefile
		$(CC) $(CFLAGS)  -o $(BIN)/zebedde $(COBJS) $(LFLAGS)
		echo "Made Latest version"

clean:
	rm -f *.o core 

tidy:
	rm -f *rst  *prm *cor *log *~

clobber:
	rm -f $(COBJS) *.o core

test: test.o size_vector.o  matrix_vec_mult.o fract_to_cart.c cart_to_fract.c  makefile
	$(CC) $(CFLAGS) -o test test.o size_vector.o  matrix_vec_mult.o \
          fract_to_cart.o cart_to_fract.o $(LFLAGS)

test_random : test_random.o pick_frm_wgt_list.o real_random.o
	 $(CC) $(CFLAGS) -o test_random test_random.o pick_frm_wgt_list.o \
         real_random.o $(LFLAGS)

test_str : test_str.o compare_strings.o copy_strings.o
	$(CC) $(CFLAGS) -o test_str test_str.o compare_strings.o copy_strings.o $(LFLAGS)

TEST_POTS_C = test_pots.c open_file.c read_line.c locate_string.c int_from_string.c \
              get_int.c get_doub.c string_from_int.c next_space.c next_none_space.c \
get_nonbond_info.c get_nonbond_params.c

TEST_POTS = test_pots.o open_file.o read_line.o locate_string.o int_from_string.o \
            get_int.o get_doub.o string_from_int.o next_space.o next_none_space.o \
get_nonbond_info.o get_nonbond_params.o

test_pots : $(TEST_POTS)
	$(CC) $(CFLAGS) -o test_pots $(TEST_POTS) 

touch_pots :
	touch $(TEST_POTS_C)
