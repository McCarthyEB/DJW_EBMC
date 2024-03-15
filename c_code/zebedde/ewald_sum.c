/****************************************************************************************/
/***** Overseer code for Ewald Sum ******************************************************/
/***** Started 12 June 1996 Dave Willock ************************************************/
/****************************************************************************************/

#include <stdio.h>
#include "maxima.h"
#include "structures.h"
#include "ewald.h"
#include "own_maths.h"

double ewald_real_part(atom *p_mol_atom, atom *p_crystal_cell, int num_cell_atoms, int self);

double ewald_recip_part(atom *p_crystal_cell, int num_cell_atoms, double *p_kvecs, 
                        double *p_kvec2, double *p_gvec2, int num_kvecs, atom *p_mol_atom, 
                        double *p_cos_sum, double *p_sin_sum, 
                        int self, int using_host);

void ewald_sum(atom *p_molecule, int num_molecule_atoms, 
               atom *p_crystal_cell, int num_cell_atoms, 
               double *p_kvecs, double *p_kvec2,
               double *p_gvec2, int num_kvecs,
               double *p_cos_sum, double *p_sin_sum,
               int zero_first, int self, int using_host)
{
#include "header.h"
int mol_index;
atom *p_mol_atom;

for (mol_index=0; mol_index <= num_molecule_atoms; mol_index++)
  {
    p_mol_atom= p_molecule+mol_index;

    if (zero_first) p_mol_atom->electrostatic_pot = 0;

    p_mol_atom->electrostatic_pot += ewald_real_part(p_mol_atom, p_crystal_cell, num_cell_atoms, self); 

    p_mol_atom->electrostatic_pot += ewald_recip_part(p_crystal_cell, num_cell_atoms, p_kvecs, 
                                                      p_kvec2, p_gvec2, num_kvecs, p_mol_atom, 
                                                      p_cos_sum, p_sin_sum, self, using_host); 

  }

return;
}
