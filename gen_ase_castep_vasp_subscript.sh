#!/bin/bash --login
#SBATCH -o Scaletest.o         # Job output file
#SBATCH -e Scaletest.e         # Job error file
#SBATCH -J Scaletest          # Job name
#SBATCH -p dev                # Queue partition
#SBATCH --ntasks=40           # number of parallel processes (tasks)
#SBATCH --ntasks-per-node=40  # tasks to run per node
#SBATCH --time=01:00:00       # time limit
#SBATCH --exclusive           # exclusive node access
#SBATCH -A scw1161            # Project code
# 
# Specify the code you are using: CASTEP or VASP
# Define MACHINE you are working on   - change for your system
# Define ASE_SCRIPT for the script that will control this job
#
# Define list of structure-directories that contain jobs to run and cores to be used for
# each. This can be used to run several jobs in parallel in a single job submission.
# For multiple runs with this script set CORES_PER_TASK for the number to assign to
# each instance of the code.
# If you need all the cores for one job just enter one item in the list
#
# Also define a sub_dir where this job has to be done. This will allow definitions of
# opt, bond_scan, dimer, freq etc for the series of jobs in our work flow.
#
# zfix gives the Z-co-ordinate below which atoms should be frozen for slabs with fixed lower sections.
# If this is not such a calculation set zfix="NoFix"
#
export CODE="VASP"
MACHINE='hawk'
ASE_SCRIPT='eos_Ag_vasp.py'
struct_list="Ag_bulk"
CORES_PER_TASK=80
#
struct_name="POSCAR_0"
#param_name="silica_standard_TS.param"
sub_dir="eos_PBE"
zfix="NoFix"
#
# CASTEP set up for each machine 
#
if [[ $CODE == "CASTEP" ]]; then
#
   echo "This is a CASTEP run"
#
   if [[ $MACHINE == "archer2" ]]; then
     echo "Importing modules and setting CASTEP_COMMAND for archer2"
# Move to directory that script was submitted from
     JOBID=$SLURM_JOBID
     export CASTEP_COMMAND='srun --distribution=block:block --hint=nomultithread castep.mpi'
#
     LAUNCH_DIR=$PWD
#
# Load the CASTEP module, avoid any unintentional OpenMP threading by
# setting OMP_NUM_THREADS, and launch the code.
# Setup the batch environment
     module load epcc-job-env
     module load castep/20.1.1-gcc10-mkl-cpe2103
     module load cray-python
# Make sure we pick up local installation of ase
     export PYTHONUSERBASE=/work/e05/e05/$USER/.local
     export PATH=$PYTHONUSERBASE/bin:$PATH
     export PYTHONPATH=$PYTHONUSERBASE/lib/python3.8/site-packages:$PYTHONPATH
#
# Location of pseudo potentials and standard param files for castep
#
     export OMP_NUM_THREADS=1
     export CASTEP_PP_PATH=/work/e05/e05/$USER/progs/castep/pseudo_pots
     export PARAM_FILES=/work/e05/e05/$USER/progs/castep/param_files
#
     CORES_PER_TASK="Not_used_on_archer2"

#
   elif [[ $MACHINE == "young" ]]; then
     echo "Importing modules and setting CASTEP_COMMAND for young"
#
     export CASTEP_COMMAND='mpirun castep.mpi '
     module unload default-modules/2018
     module unload compilers mpi
     module load compilers/intel/2019/update4
     module load mpi/intel/2019/update4/intel
     module load castep/19.1.1/intel-2019
# To use ase will need python
     module load python/3.9.6

     CORES_PER_TASK="Not_used_on_Young"
#
   elif [[  $MACHINE == "hawk" ]]; then
     echo "Importing modules and setting CASTEP_COMMAND for hawk"

     export CASTEP_PP_PATH=$HOME/progs/castep/pseudo_pots
     export PARAM_FILES=$HOME/castep/param_files
     export CASTEP_COMMAND='mpirun castep.mpi '

     export OMP_NUM_THREADS=1
     export I_MPI_ADJUST_ALLTOALLV=2

     JOBID=$SLURM_JOBID
     LAUNCH_DIR=$PWD

     NNODES=$SLURM_NNODES
     NCPUS=$SLURM_NTASKS
     PPN=$SLURM_NTASKS_PER_NODE
   fi
#
# VASP set up for each machine 
#
elif [[ $CODE == "VASP" ]]; then
#
# ASE builds the POTCAR file itself from the set available from this path:
#               The pseudopotentials are expected to be in:
#               LDA:  $VASP_PP_PATH/potpaw/
#               PBE:  $VASP_PP_PATH/potpaw_PBE/
#               PW91: $VASP_PP_PATH/potpaw_GGA/
#
   echo "This is a VASP run"
#
   if [[ $MACHINE == "archer2" ]]; then
# Move to directory that script was submitted from
      export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)
      cd $PBS_O_WORKDIR
      JOBID=$SLURM_JOBID
#
      module load vasp/5
      module load python-compute/3.6.0_gcc6.1.0

      CORES_PER_TASK="Not_used_on_archer2"
#
# On Archer2 PAW are availble at:
# $VASP_PSPOT_DIR
#     
   elif [[ $MACHINE == "young" ]]; then
#
      LAUNCH_DIR=$PWD
      MACHINE='young'
      JOBID=$JOB_ID
      export VASP_PP_PATH=$HOME/progs/vasp
#
      CORES_PER_TASK="Not_used_on_young"
#
   elif [[  $MACHINE == "hawk" ]]; then
      export OMP_NUM_THREADS=1
      ulimit -s unlimited
      JOBID=$SLURM_JOBID
      LAUNCH_DIR=$PWD
#
      export VASP_PP_PATH=$HOME/progs/vasp

      NNODES=$SLURM_NNODES
      NCPUS=$SLURM_NTASKS
      PPN=$SLURM_NTASKS_PER_NODE
#
# Do the module loads
#
      module purge
# load the VASP module, it loads required compiler, mpi and mkl libraries
      module load  vasp/5.4.4
# To use ase will need python
      module load python
#
   fi
fi
#
echo LAUNCH_DIR set to $LAUNCH_DIR
cd $LAUNCH_DIR
#
#export OMP_NUM_THREADS=1
#export I_MPI_ADJUST_ALLTOALLV=2
export LOG_FILE="$LAUNCH_DIR"/"$ASE_SCRIPT"_"$JOBID".log  

###################################
#           WORKDIR               #
# Young also needs a tmp dir      #
###################################
# run the simulation on /scratch
# Make a directory for this run, this will be the top directory used
# for the working area of the run.
# Make a directory for this run
if [[ $MACHINE == "hawk" ]]; then
    export work_dir=/scratch/$USER/$JOBID
#    rm -rf $work_dir
    mkdir -p $work_dir
#
elif [[ $MACHINE == "young" ]]; then
    export work_dir=$HOME/Scratch/workspace/$JOBID
    rm -rf $work_dir
    mkdir -p $work_dir
# On Young also need a tmp directory for this run
    export tmp_dir=$work_dir/tmp
    mkdir -p $tmp_dir
    export FORT_TMPDIR=$tmp_dir
#
elif [[ $MACHINE == "archer2" ]]; then
    export work_dir=$LAUNCH_DIR/workspace/$JOBID
    rm -rf $work_dir
    mkdir -p $work_dir
# On Young also need a tmp directory for this run
    export tmp_dir=$work_dir/tmp
    mkdir -p $tmp_dir
    export FORT_TMPDIR=$tmp_dir
else
    echo ERROR MACHINE $MACHINE not recognised when setting up workspace
    exit 0
fi
#
#
echo Running on machine: $MACHINE     >  $LOG_FILE 
echo Running on host `hostname`       >> $LOG_FILE 
echo Time is `date`                   >> $LOG_FILE  
echo Launch directory is $LAUNCH_DIR    >> $LOG_FILE
echo Working directory is $work_dir   >> $LOG_FILE
echo Directory is `pwd`               >> $LOG_FILE  
echo Job ID is $JOBID                 >> $LOG_FILE  
echo Using ase script $ASE_SCRIPT     >> $LOG_FILE
echo Structure file is $struct_name   >> $LOG_FILE
if [[ $CODE == "CASTEP" ]]; then
  echo CASTEP_PP_PATH set to $CASTEP_PP_PATH >> $LOG_FILE 
  echo Castep param file is $PARAM_FILES/$param_name >> $LOG_FILE
fi
if [[ $MACHINE == "hawk" ]]; then
  echo This jobs runs on the following machine: `echo $SLURM_JOB_NODELIST | uniq` >> $LOG_FILE 
  echo Number of Processing Elements is $NCPUS >> $LOG_FILE 
  echo Number of mpiprocs per node is $PPN   >> $LOG_FILE 
fi
echo >> $LOG_FILE 
echo Start Time is `date` running NCPUs=$NCPUS PPN=$(ppn) >> $LOG_FILE 
start="$(date +%s)"
echo                                   >> $LOG_FILE
#
# ... change to the scratch working directory       
cd $work_dir
#
# Loop over the struct_list
#
for struct_dir in $struct_list; do
     echo job = $struct_dir >> $LOG_FILE
     echo sub_dir = $sub_dir >> $LOG_FILE
     export jobwork_dir="$work_dir"/"$struct_dir"/"$sub_dir"
     export results_dir="$LAUNCH_DIR"/"$struct_dir"/"$sub_dir"
     echo Will be working in directory $jobwork_dir >> $LOG_FILE
     echo results to go to directory   $results_dir >> $LOG_FILE
#
# Check if the sub_dirs already exists
#
     if [ ! -d $jobwork_dir ]; then
       echo Creating working directory $jobwork_dir >> $LOG_FILE
       mkdir -p $jobwork_dir
     fi
#
     if [ ! -d $results_dir ]; then
       echo Creating results directory $results_dir >> $LOG_FILE
       mkdir -p $results_dir
     fi
# We are already in the work area, now move into the sub_dir
     cd $jobwork_dir
#
     if [[ $CODE == "CASTEP" ]]; then
         export full_param="$PARAM_FILES"/"$param_name"
#
       if [ "$full_param" ]; then
          if [ -e  "$full_param" ]; then
            echo Copying param file $full_param to  $jobwork_dir >> $LOG_FILE
            cp "$full_param" struct.param
          else
            echo param file "$full_param" file missing or mis-named >> $LOG_FILE
          fi
       fi
     fi
#
# Test what type of file we have defined and then copy to the appropriate standard
#
     export full_struct="$LAUNCH_DIR"/"$struct_dir"/"$struct_name"
#
     if [ -e  "$full_struct" ]; then
       echo Copying structure file $full_struct to  $jobwork_dir >> $LOG_FILE
#
       if [[ $struct_name == *"POSCAR"* ]]; then
          cp "$full_struct" POSCAR
       elif  [[ $struct_name == *".cell" ]]; then
          cp "$full_struct" struct.cell
       elif  [[ $struct_name == *".car" ]]; then
          cp "$full_struct" struct.car
       else
          echo ERROR "$struct_name" not recognised as one of the standards for VASP or CASTEP
          exit 0
       fi
#
     else
      echo struct file "$full_struct" file missing or mis-named >> $LOG_FILE
     fi
# copy the script to work area and to results_dir for reference later
#
     if [ -e  "$LAUNCH_DIR"/"$ASE_SCRIPT" ]; then
       cp "$LAUNCH_DIR"/"$ASE_SCRIPT"  .
       cp "$LAUNCH_DIR"/"$ASE_SCRIPT"  $results_dir
     else
       echo ase python script "$ASE_SCRIPT" file missing from "$LAUNCH_DIR" >> $LOG_FILE
     fi
#
     echo $PWD Running python script $ASE_SCRIPT >> $LOG_FILE
     echo  >> $LOG_FILE
     echo With command: >> $LOG_FILE
     echo python3 $ASE_SCRIPT $jobwork_dir $CORES_PER_TASK $JOBID $LAUNCH_DIR $struct_dir $sub_dir  \
                    $MACHINE $zfix > "$LAUNCH_DIR"/"$struct_dir"/"$sub_dir"/"$ASE_SCRIPT"_"$JOBID".out >> $LOG_FILE
#
#  run the ase script
#
     python3 $ASE_SCRIPT $jobwork_dir $CORES_PER_TASK $JOBID $LAUNCH_DIR $struct_dir $sub_dir $MACHINE $zfix \
                                    > "$LAUNCH_DIR"/"$struct_dir"/"$sub_dir"/"$ASE_SCRIPT"_"$JOBID".out &

     if [[ $CODE == "VASP" ]]; then
        echo VASP  run using ase script $ASE_SCRIPT for job $struct_dir running. >> $LOG_FILE
     elif [[ $CODE == "CASTEP" ]]; then
        echo CASTEP run using ase script $ASE_SCRIPT for job $struct_dir running. >> $LOG_FILE
     fi

     cd ..
   done
#
# Wait until all background jobs complete
#
wait
#
# Copy files back
#
   for struct_dir in $struct_list; do
#
     export jobwork_dir="$work_dir"/"$struct_dir"/"$sub_dir"
     echo "copying back for " $jobwork_dir >> $LOG_FILE
     cd $jobwork_dir
#
# Remove the "Processing" lines from the out file
#
     export home_sub_dir=${LAUNCH_DIR}/$struct_dir/"$sub_dir"
     awk '$1 !~ /Proc/ {print}' "$home_sub_dir"/"$ASE_SCRIPT"_"$JOBID".out > temp.out
     mv temp.out "$home_sub_dir"/"$ASE_SCRIPT"_"$JOBID".out

     cp * $home_sub_dir
#
# Only copy particular files from the CASTEP directory
#
     if [[ $CODE == "CASTEP" ]]; then
        mkdir $home_sub_dir/castep_"$JOBIB"
        cp CASTEP/castep.castep $home_sub_dir/castep_"$JOBID"
        cp CASTEP/castep.bands  $home_sub_dir/castep_"$JOBID"
        cp CASTEP/castep.cell   $home_sub_dir/castep_"$JOBID"
        cp CASTEP/castep.params $home_sub_dir/castep_"$JOBID"
     fi
#
     cd ..
#
done
#
# clean scratch
#

#rm -r $top_scratch 


# Record Total time for the job
stop="$(date +%s)"
finish=$(( $stop-$start ))
echo  $JOBID  Job-Time  $finish seconds
