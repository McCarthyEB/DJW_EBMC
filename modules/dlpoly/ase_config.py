from ase import Atom
from ase import Atoms
from ase.io import read, write, Trajectory
#
#         0         3
#       12.4000000000        0.0000000000        0.0000000000
#       -0.0000000010       12.4000000000        0.0000000000
#       -0.0000000010       -0.0000000010       12.4000000000
#ow      
#    0.0000000000        0.0000000000        1.5545667730
#hw      
#    0.7569503270        0.0000000000        2.1404490490
#hw      
#
#
def read_config_file(config_file):
#
   fptr=open(config_file, 'rt')
   lines=fptr.readlines()
#
   is_lab=True
   cell=[]
   symbols=[]
   pos=[]
   skip=False
   for iline in range(0,len(lines)):
     words=lines[iline].split()
     nwords=len(words)
#
# After two lines get cell vectors
#
     if (iline > 1 and iline < 5):
        cell.append([float(words[0]),float(words[1]),float(words[2])])
#
#
     if ( iline > 4 ):
       if is_lab:
         if ("ow" in words[0]):
            symbols.append("O")
         elif ("hw" in words[0]):
            symbols.append("H")
         else:
            skip=True
#
         is_lab=False
       else:
         if not skip:
            pos.append([float(words[0]),float(words[1]),float(words[2])])
         skip=False
         is_lab=True
#
   return cell,symbols,pos

#
# Main code begins
#
config_file="CONFIG_last_500K"
cell,symbols,pos=read_config_file(config_file)
#
print(cell)
#
atoms=Atoms(symbols,pos)
atoms.set_cell(cell)
#
for iatom in range(0,len(atoms)):
    print(f"%d %s %10.6f %10.6f %10.6f"                              \
          % (iatom, atoms.symbols[iatom], atoms.positions[iatom][0], \
             atoms.positions[iatom][1], atoms.positions[iatom][2]))
#
write('check_config.cif', atoms ,format='cif')
#


     
       

