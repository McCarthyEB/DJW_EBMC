#!/bin/bash

# Request 2 nodes (256 MPI tasks at 128 tasks per node) for 20 minutes.   

#SBATCH --job-name=VASP_test
#SBATCH --nodes=4
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00

# Replace [budget code] below with your project code (e.g. t01)
#SBATCH --account=e05-free
#SBATCH --partition=standard
#SBATCH --qos=standard

# Setup the job environment (this module needs to be loaded before any other modules)
module load epcc-job-env

# Load the VASP module, avoid any unintentional OpenMP threading by
# setting OMP_NUM_THREADS, and launch the code.
export OMP_NUM_THREADS=1
module load vasp/5
module load cray-python
#
# Ensure python knows about ASE
#
# python location to have ASE installation on work for the queue nodes to see
#
export PYTHONUSERBASE=/work/e05/e05/$USER/.local
export PATH=$PYTHONUSERBASE/bin:$PATH
#
# ASE builds the POTCAR file itself from the set available from this path:
#               The pseudopotentials are expected to be in:
#               LDA:  $VASP_PP_PATH/potpaw/
#               PBE:  $VASP_PP_PATH/potpaw_PBE/
#               PW91: $VASP_PP_PATH/potpaw_GGA/
#
export VASP_PP_PATH=/work/e05/e05/$USER/progs/vasp
#
# Define directory job submitted from - change for your system
# Define MACHINE you are working on   - change for your system
# Define ASE_SCRIPT for the script that will control this job  
LAUNCH_DIR=$PWD
MACHINE='archer2'
JOBID=$SLURM_JOBID
ASE_SCRIPT='ase_vasp_opt.py'
#
# Define list of structure-directories that contain jobs to run and cores to be used for
# each. This can be used to run several jobs in parallel in a single job submission.
# If you need all the cores for one job just enter one item in the list
#
# Also define a sub_dir where this job has to be done. This will allow definitions of 
# opt, bond_scan, dimer, freq etc for the series of jobs in our work flow.
#
struct_list="BLANK"
poscar_name="POSCAR_latest"
sub_dir="opt"

CORES_PER_TASK="Not_used_on_archer2"

export I_MPI_ADJUST_ALLTOALLV=2
export LOG_FILE="$LAUNCH_DIR"/"$ASE_SCRIPT"_"$JOBID".log  

NNODES=$SLURM_NNODES
NCPUS=$SLURM_NTASKS
PPN=$SLURM_NTASKS_PER_NODE

###################################
#           WORKDIR               #
###################################
# On ARCHER we are already in the work area
# when we launch but still define work_dir for
# consistency with other machines
# 
export work_dir=$LAUNCH_DIR

echo Running on machine: $MACHINE     >  $LOG_FILE 
echo Running on host `hostname`       >> $LOG_FILE 
echo Time is `date`                   >> $LOG_FILE  
echo Launch directory is $LAUNCH_DIR  >> $LOG_FILE
echo Working directory is $work_dir   >> $LOG_FILE
echo Current directory is `pwd`       >> $LOG_FILE  
echo Job ID is $JOBID                 >> $LOG_FILE  
echo Using ase script $ASE_SCRIPT     >> $LOG_FILE
echo This jobs runs on the following machine: `echo $SLURM_JOB_NODELIST | uniq` >> $LOG_FILE 
echo Number of Processing Elements is $NCPUS >> $LOG_FILE 
echo Number of mpiprocs per node is $PPN     >> $LOG_FILE 
echo >> $LOG_FILE 
echo VASP_PP_PATH set to $VASP_PP_PATH >> $LOG_FILE 
echo Script Start Time is `date` running NCPUs=$NCPUS PPN=$PPN >> $LOG_FILE 
start="$(date +%s)"
#
# ... change to the scratch working directory       
cd $work_dir
#
# Loop over the struct_list
#
for struct_dir in $struct_list; do
     echo job = $struct_dir >> $LOG_FILE
     echo sub_dir = $sub_dir >> $LOG_FILE
#
# Each structure will be run in ots own sub-directory 
# to hold the name of this jobwork_dir is the full path of the 
# location. 
# results_dir is the place that the results will be sent in the /home space
#
     export jobwork_dir="$work_dir"/"$struct_dir"/"$sub_dir"
     export results_dir="$LAUNCH_DIR"/"$struct_dir"/"$sub_dir"
     echo  >> $LOG_FILE
     echo Will be working in directory $jobwork_dir >> $LOG_FILE
     echo results to go to directory   $results_dir >> $LOG_FILE
#
# Check if the sub_dirs already exists
#
     echo  >> $LOG_FILE
     if [ ! -d "$jobwork_dir" ]; then
       echo Creating working directory $jobwork_dir >> $LOG_FILE
       mkdir -p $jobwork_dir
     fi
#
     echo  >> $LOG_FILE
     if [ ! -d "$results_dir" ]; then
       echo Creating results directory $results_dir >> $LOG_FILE
       mkdir -p $results_dir
     fi
#
     echo  >> $LOG_FILE
     export full_poscar="$LAUNCH_DIR"/"$struct_dir"/"$poscar_name"
     if [ -e  "$full_poscar" ]; then
       echo Copying POSCAR file $full_poscar to  $jobwork_dir >> $LOG_FILE
       cp "$full_poscar" "$jobwork_dir"/POSCAR
     else
      echo poscar "$full_poscar" file missing or mis-named >> $LOG_FILE
     fi
#
# copy the script to work area and to results_dir for reference later
#
     if [ -e  "$LAUNCH_DIR"/"$ASE_SCRIPT" ]; then
       cp "$LAUNCH_DIR"/"$ASE_SCRIPT"  $jobwork_dir
       cp "$LAUNCH_DIR"/"$ASE_SCRIPT"  $results_dir
     else
       echo ase python script "$ASE_SCRIPT" file missing from "$LAUNCH_DIR" >> $LOG_FILE
     fi
#
# Now move into the sub_dir and run the task in the background
     cd $jobwork_dir
#
     echo Moved to directory : $PWD 
     echo To run python script $ASE_SCRIPT >> $LOG_FILE
     echo  >> $LOG_FILE
     echo With script arguements : >> $LOG_FILE
     echo Arg. 1, jobwork_dir    : $jobwork_dir    >> $LOG_FILE
     echo Arg. 2, CORES_PER_TASK : $CORES_PER_TASK >> $LOG_FILE
     echo Arg. 3, JOBID          : $JOBID          >> $LOG_FILE
     echo Arg. 4, LAUNCH_DIR     : $LAUNCH_DIR     >> $LOG_FILE
     echo Arg. 5, struct_dir     : $struct_dir     >> $LOG_FILE
     echo Arg. 6, sub_dir        : $sub_dir        >> $LOG_FILE
     echo Arg. 7, MACHINE        : $MACHINE        >> $LOG_FILE
     echo  >> $LOG_FILE
     echo using command: >> $LOG_FILE
     echo python3 $ASE_SCRIPT $jobwork_dir $CORES_PER_TASK $JOBID $LAUNCH_DIR $struct_dir $sub_dir  \
                    $MACHINE > "$results_dir"/"$ASE_SCRIPT"_"$JOBID".out >> $LOG_FILE
#
#  run the ase script
#
     python3 $ASE_SCRIPT $jobwork_dir $CORES_PER_TASK $JOBID $LAUNCH_DIR $struct_dir $sub_dir $MACHINE \
                                    > "$results_dir"/"$ASE_SCRIPT"_"$JOBID".out &

     echo  >> $LOG_FILE
     echo vasp run using ase script $ASE_SCRIPT for job $struct_dir set running. >> $LOG_FILE
     echo  >> $LOG_FILE

   done

#
# Wait until all background jobs complete

#
wait
#
# Copy files back
#
   for struct_dir in $struct_list; do

     
     export jobwork_dir="$work_dir"/"$struct_dir"/"$sub_dir"
     export results_dir="$LAUNCH_DIR"/"$struct_dir"/"$sub_dir"
     echo "copying back from " $jobwork_dir >> $LOG_FILE
     echo "copying back to   " $results_dir >> $LOG_FILE
     echo  >> $LOG_FILE
#
# Move into the work area
#
     cd $jobwork_dir
     
     cp *.py   "$results_dir"

     ls >> "$LOG_FILE"

     cp OUTCAR  "$results_dir"/OUTCAR_ase_"$JOBID"_"$sub_dir"
     cp CONTCAR "$results_dir"/CONTCAR_ase_"$JOBID"_"$sub_dir"
     cp INCAR  "$results_dir"/INCAR_ase_"$JOBID"_"$sub_dir"
     cp DOSCAR  "$results_dir"/DOSCAR_ase_"$JOBID"_"$sub_dir"
     cp POS*   "$results_dir"
#
done
#
# clean scratch
#

#rm -r $top_scratch 


# Record Total time for the job
stop="$(date +%s)"
finish=$(( $stop-$start ))
echo VASP $JOBID  Job-Time  $finish seconds
