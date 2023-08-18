#
# Convert fractional to cartessian co-ordinates
#
def frac_to_cart(fvec, cell):
#
   vec=np.zeros(3)
   vec[0]=fvec[0]*cell[0][0]+fvec[1]*cell[1][0]+fvec[2]*cell[2][0]
   vec[1]=fvec[0]*cell[0][1]+fvec[1]*cell[1][1]+fvec[2]*cell[2][1]
   vec[2]=fvec[0]*cell[0][2]+fvec[1]*cell[1][2]+fvec[2]*cell[2][2]
#
   return vec
#
# Read the lattice vectors and atom coordinates from a cell file 
# and make an atoms object
#
def read_cell_djw(fp):
#
   lines = fp.readlines()
   nlines = len(lines)
   icell=0
   cell_line=0
   contents_line=0
   cell=np.zeros((3,3))
#
   have_unit_cell=False
   have_atoms=False
   have_cell_cons=False
   have_kpoint_grid=False
   have_fix_com=False
   first_atom=True
   natoms=0
   num_atoms=-1
   num_intermed_atoms=-1
   num_prod_atoms=-1
   num_cons=-1
#
   spec_mass_lab=[]
   spec_mass=[]
#
   spec_pot_lab=[]
   spec_pot=[]
#
   spec_lcao_lab=[]
   spec_lcao=[]
#
#
   read_lines=False
   read_cell=False
   read_atoms=False
   read_cell_cons=False
   read_ion_cons=False
   read_spec_mass=False
   read_spec_pot=False
   read_spec_lcao=False
#
   have_cell=False
   have_atoms=False
   have_cell_cons=False
   have_ion_cons=False
   have_masses=False
   have_pots=False
   have_lcao=False
   have_kpoint_grid=False
   have_fix_com=False
#
   for iline in range(nlines):
#
# Look for the stress lines
#
      line=lines[iline]
      linesplit=line.split()
      nwords=len(linesplit)
#
      if ( nwords >= 2 and  "endblock" in linesplit[0].lower ):
         print("Found ENDBLOCK: %s" % line)
         read_lines=False
         read_cell=False
         read_atoms=False
         read_cell_cons=False
         read_ion_cons=False
         read_spec_mass=False
         read_spec_pot=False
         read_spec_lcao=False
#
      if read_cell:
         print("Reading cell line: %s" % line)
         cell[icell][0] = float(linesplit[0])
         cell[icell][1] = float(linesplit[1])
         cell[icell][2] = float(linesplit[2])
         icell+=1
#
      elif read_atoms:
#
         label=linesplit[0]
#
         if is_frac:
            fracs=[float(linesplit[1]),float(linesplit[2]),float(linesplit[3])]
            coords=frac_to_cart(fracs, cell)
#
         else:
            coords=[float(linesplit[1]),float(linesplit[2]),float(linesplit[3])]
#
         new_atom=Atom(label,coords)
#
         if first_atom:
            print("Reading atoms")
            print("First atom line:")
            print(line)
            first_atom=False
#
            if is_intermed:
              atoms_intermed=Atoms([new_atom], cell=cell)
              num_intermed_atoms+=1
            elif is_prod:
              atoms_prod=Atoms([new_atom], cell=cell)
              num_prod_atoms+=1
            else:
              atoms=Atoms([new_atom], cell=cell)
              num_atoms+=1
         else:
            if is_intermed:
              atoms_intermed.append(new_atom)
              num_intermed_atoms+=1
            elif is_prod:
              atoms_prod.append(new_atom)
              num_prod_atoms+=1
            else:
              atoms.append(new_atom)
              num_atoms+=1
#
      elif read_cell_cons:
         have_cell_cons=True
#     
         cell_cons.append(int(linesplit[0]))
         cell_cons.append(int(linesplit[1]))
         cell_cons.append(int(linesplit[2]))
#
      elif read_spec_mass:    
        have_masses=True
        print("Reading masses.....")
        spec_mass_lab.append(linesplit[0])
        spec_mass.append(float(linesplit[1]))
#
      elif read_spec_pot:
        have_pots=True
        spec_pot_lab.append(linesplit[0])
        spec_pot.append(linesplit[1])
#
      elif read_spec_lcao:
        have_lcao=True
        spec_lcao_lab.append(linesplit[0])
        spec_lcao.append(int(linesplit[1]))
#
#
# Check for BLOCK flags                
#
      if ( nwords >= 2 and 'block' in linesplit[0].lower and 'end' not in linesplit[0].lower ):
         read_lines=True
         print("Looking at BLOCK : %s" % line)
#
# Read second word    
#
         if ( read_lines and "lattice_cart" in linesplit[1].lower ):
           read_cell=True
           icell=0
#
         elif ( read_lines and "positions" in linesplit[1].lower ):
           read_atoms=True
           is_frac=False
           is_intermed=False
           is_prod=False
           if "frac" in linesplit[1].lower:
              is_frac=True
           if "intermediate" in linesplit[1].lower:
              is_intermed=True
              first_atom=True
           elif "PRODUCT" in linesplit[1]:
              is_prod=True
              first_atom=True
#
         elif ( read_lines and "cell_constraints" in linesplit[1].lower ):
           read_cell_cons=True
           cell_cons=[]
#
         elif ( read_lines and "ionic_constraints" in linesplit[1].lower ):
           read_ion_cons=True
#
         elif ( read_lines and "species_mass" in linesplit[1].lower ):
           print("Will read masses...........")
           read_spec_mass=True
#
         elif ( read_lines and "species_pot" in linesplit[1].lower ):
           read_spec_pot=True
#
         elif ( read_lines and "species_lcao_states" in linesplit[1].lower ):
           read_spec_lcao=True
#
# Pick up one liners
#
#kpoint_mp_grid 1 1 1
#or kpoint_mp_grid : 1 1 1
#FIX_COM : false
      elif ( nwords >= 2 and 'kpoint_mp_grid' in linesplit[0].lower):
         have_kpoint_grid=True
         if ":" in linesplit[1]:
            kpoint_mp=[int(linesplit[2]),int(linesplit[3]),int(linesplit[4])]
         else:
            kpoint_mp=[int(linesplit[1]),int(linesplit[2]),int(linesplit[3])]
#
      elif ( nwords >= 2 and 'fix_com' in linesplit[0].lower):
         have_fix_com=True
         fix_com=linesplit[2]
#
   return atoms
#
