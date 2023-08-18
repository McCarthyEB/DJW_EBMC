import numpy as np
from ase import Atoms
from ase.io.trajectory import Trajectory
from ase.calculators.aims import Aims

a = 4.0 # approximate lattice constant
b = a / 2


calc = Aims(xc_pre=['pbe', '10'], #we do here 10 steps with PBE to stabilize the system
           override_warning_libxc="true", #MBEEF is not impplemented directly with FHI, for that we need to call libxc library
           xc='libxc MGGA_X_MBEEF+GGA_C_PBE_SOL', #MBEEF
           spin='none',
           k_grid=(20,20,20),   
           relativistic=('atomic_zora', 'scalar'),
           compute_forces="true",
           sc_accuracy_etot=1e-3,
           sc_accuracy_eev=1e-3,
           sc_accuracy_rho=1e-6,
           sc_accuracy_forces=1e-4,
           )

bulk = Atoms('Pd',
           cell=[(0, b, b), (b, 0, b), (b, b, 0)],
           pbc=1,
           calculator=calc)  
cell = bulk.get_cell()
traj = Trajectory('bulk.traj', 'w')
for x in np.linspace(0.90, 1.05, 5):
    bulk.set_cell(cell * x, scale_atoms=True)
    bulk.get_potential_energy()
    traj.write(bulk)

########## EOS #############
from ase.io import read
from ase.units import kJ
from ase.eos import EquationOfState
from ase.io import write

configs = read('bulk.traj@0:5')  # read 5 configurations

# Extract volumes and energies:

volumes = [bulk.get_volume() for bulk in configs]
energies = [bulk.get_potential_energy() for bulk in configs]
eos = EquationOfState(volumes, energies, eos='birchmurnaghan')
v0, e0, B = eos.fit()

print("Volume: ", v0)
print("energy of the bulk", e0)
print("FCC Lattice parameter a0: ", (4*v0)**(1/3), "Angstrom")
print("Bulk modulus: ", B / kJ * 1.0e24, 'GPa')
#eos.plot('bulk-eos-Pd.png') # no ase gui in haak