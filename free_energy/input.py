from ase.io import read
from ase.calculators.aims import Aims
#from ase.optimize import QuasiNewton
#from ase.vibrations import Vibrations
from ase.vibrations import Infrared
from ase.thermochemistry import IdealGasThermo

# Set environment variables needed by ASE - check basis sets
import os
os.environ['ASE_AIMS_COMMAND']="mpirun -np "+os.environ['SLURM_NTASKS']+" /home/scw1057/software/fhi-aims/bin/aims."+os.environ['VERSION']+".scalapack.mpi.x"
os.environ['AIMS_SPECIES_DIR']="/home/scw1057/software/fhi-aims/species_defaults/light"   # Light settings
#os.environ['AIMS_SPECIES_DIR']="/home/c.sacal6/software/fhi-aims-species-defaults/tight" # Tight settings

molecule = read('geometry.in')
#molecule.set_pbc(False)

# Set calculator for FHI-aims
calc = Aims(xc='pbe',
            vdw_correction_hirshfeld=True,
            spin='collinear',
            default_initial_moment=0.01,
            relativistic=('atomic_zora','scalar'),
            charge=0.,
            occupation_type=('gaussian','0.1'),
            n_max_pulay=20,
            charge_mix_param=0.05,
            sc_accuracy_etot=1e-6,
            sc_accuracy_eev=1e-6,
            sc_accuracy_rho=1e-6,
            sc_accuracy_forces=1e-6,
            output=['dipole'],
            restart_read_only='wavefunction',
            restart_write_only='wavefunction')

molecule.set_calculator(calc)

#Disabled as structure pre-optimised with internal FHI-aims method
#dynamics = QuasiNewton(molecule, trajectory='molecule_opt.traj')
#dynamics.run(fmax=0.01)

potentialenergy = molecule.get_potential_energy()

#vib = Vibrations(molecule)
vib = Infrared(molecule)
vib.run()
vib.summary()
vib_energies = vib.get_energies()

thermo = IdealGasThermo(vib_energies=vib_energies,
                        potentialenergy=potentialenergy,
                        atoms=molecule,
                        geometry='nonlinear',
                        spin=0,
                        symmetrynumber=1)

G = thermo.get_gibbs_energy(temperature=298.15, pressure=101325.)

