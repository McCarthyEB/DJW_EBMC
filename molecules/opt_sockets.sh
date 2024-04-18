#!/bin/bash --login
#SBATCH -o dihed.o.%J         # Job output file
#SBATCH -e dihed.e.%J         # Job error file
#SBATCH -J dihed              # Job name
#SBATCH -p dev                # Queue partition
#SBATCH --ntasks=40           # number of parallel processes (tasks)
#SBATCH --ntasks-per-node=40  # tasks to run per node
#SBATCH --time=00:50:00       # time limit
#SBATCH --exclusive           # exclusive node access
#SBATCH -A scw1161	      # Project code

## Usage
# sbatch script # submit job
# squeue        # job status
ulimit -s unlimited
#module purge
# load the VASP module, it loads required compiler, mpi and mkl libraries
module load mpi
# To use ase will need python
module load python/3.7.7-intel2020u1
#module load python
#module load ase/3.20.1 python/3.7.0 pymatgen/2018.11.30
module load ase/3.20.1
#
# Add carmm to python path to access socket versions
export PYTHONPATH=~/python/carmm:$PYTHONPATH
echo modules:
module list
#export VERSION=230612
export VERSION=231208
#
# Define directory job submitted from - change for your system
# Define MACHINE you are working on   - change for your system
# Define ASE_SCRIPT for the script that will control this job  
LAUNCH_DIR=$PWD
MACHINE='hawk'
JOBID=$SLURM_JOBID
ASE_SCRIPT='dihedral_scan_example.py'
#
# Define a list of stems for structures that need to be optimised. These should have
# corresponding "$stem".in files in $work_dir
#
stem_list="C5H4O2"
#
# Define the cores per task to use, this should be the closest integer to $SLURM_NTASKS/num_stems
#
CORES_PER_TASK=40
#
# run the simulation on /scratch
export work_dir=/scratch/"$USER"/"$JOBID"
mkdir -p $work_dir
#
#export ASE_AIMS_COMMAND="mpirun -np $SLURM_NTASKS /home/scw1057/software/fhi-aims/bin/aims.210513.scalapack.mpi.x"
export ASE_AIMS_COMMAND=mpirun -np $CORES_PER_TASK /home/scw1057/software/fhi-aims/bin/aims.$VERSION.scalapack.mpi.x
export AIMS_SPECIES_DIR="/home/scw1057/software/fhi-aims/species_defaults/defaults_2020/light"
#
export OMP_NUM_THREADS=1
export I_MPI_ADJUST_ALLTOALLV=2
export LOG_FILE="$LAUNCH_DIR"/"$ASE_SCRIPT"_"$JOBID".log  

NNODES=$SLURM_NNODES
NCPUS=$SLURM_NTASKS
PPN=$SLURM_NTASKS_PER_NODE

echo Running on machine: $MACHINE   >  $LOG_FILE 
echo Running on host `hostname`     >> $LOG_FILE 
echo Time is `date`                 >> $LOG_FILE  
echo Launch directory is $LAUNCH_DIR    >> $LOG_FILE
echo Working directory is $work_dir >> $LOG_FILE
echo Directory is `pwd`             >> $LOG_FILE  
echo Job ID is $JOBID               >> $LOG_FILE  
echo Using ase script $ASE_SCRIPT   >> $LOG_FILE
echo This jobs runs on the following machine: `echo $SLURM_JOB_NODELIST | uniq` >> $LOG_FILE 
echo Number of Processing Elements is $NCPUS >> $LOG_FILE 
echo Number of mpiprocs per node is $PPN >> $LOG_FILE 
echo >> $LOG_FILE 
echo FHI-aims species path set to $AIMS_SPECIES_DIR >> $LOG_FILE 
echo Start Time is `date` running NCPUs=$NCPUS PPN=$PPN >> $LOG_FILE 
echo Cores per task set to $CORES_PER_TASK >> $LOG_FILE 
start="$(date +%s)"
#
# ... change to the scratch working directory       
cd $work_dir
pwd >> $LOG_FILE
#
# Loop over the struct_list
#
for stem in $stem_list; do
#
     echo Processing "$stem"  >> $LOG_FILE
#
# Check if the sub_dirs already exists
#
     if  ! [[ -e "$stem" ]]; then
        mkdir "$stem"
     fi
#
     struct_name="$stem".in
     echo copying "$struct_name" to "$stem" >> $LOG_FILE
     cd "$stem"
     echo Current directory:  $PWD >> $LOG_FILE
#
     cp "$LAUNCH_DIR"/"$struct_name" .
     cp "$LAUNCH_DIR"/"$ASE_SCRIPT" .
     echo "Directory listing....."
     ls >> $LOG_FILE

     echo Running python script $ASE_SCRIPT >> $LOG_FILE
     echo  >> $LOG_FILE
     echo With command: >> $LOG_FILE
     echo python3 $ASE_SCRIPT  $struct_name $CORES_PER_TASK "> $stem"_ase.out >> $LOG_FILE
#
#  run the ase script
#
#     python3 $ASE_SCRIPT $struct_name > aims.out 2>&1  &
#
# Running sequentially to avoid socket clashes
#    python3 $ASE_SCRIPT $struct_name $CORES_PER_TASK > "$stem"_ase.out 2>&1 &
     python3 $ASE_SCRIPT $struct_name $CORES_PER_TASK > "$stem"_ase.out 2>&1 &

     echo fhi_aims set running using ase script $ASE_SCRIPT for structure $struct_name. >> $LOG_FILE

     cd $work_dir
done
#
# Wait until all background jobs complete
#
stdbuf -i0 -o0 -e0 jobs >> $LOG_FILE


wait
