#!/bin/bash
#PBS -q standard
#PBS -N Au2O_3          
#PBS -o Au2O_3-out.txt
#PBS -e Au2O_3-err.txt
#PBS -l select=1 
#PBS -l walltime=24:00:00

#PBS -A e05-react-wil

export NPROC=24 

export OMP_NUM_THREADS=1
ulimit -s unlimited

#module purge
module add vasp5
# To use ase will need python
module load python-compute/3.6.0_gcc6.1.0

# Move to directory that script was submitted from
export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)
cd $PBS_O_WORKDIR
#
# ASE builds the POTCAR file itself from the set available from this path:
#               The pseudopotentials are expected to be in:
#               LDA:  $VASP_PP_PATH/potpaw/
#               PBE:  $VASP_PP_PATH/potpaw_PBE/
#               PW91: $VASP_PP_PATH/potpaw_GGA/
#
export VASP_PP_PATH=$HOME/progs/vasp

# Define directory job submitted from - change for your system
# Define MACHINE you are working on   - change for your system
# Define ASE_SCRIPT for the script that will control this job  
LAUNCH_DIR=$PBS_O_WORKDIR
MACHINE='archer'
JOBID=$PBS_JOBID
ASE_SCRIPT='ase_vasp_opt.py'
#
# Define list of structure-directories that contain jobs to run and cores to be used for
# each. This can be used to run several jobs in parallel in a single job submission.
# If you need all the cores for one job just enter one item in the list
#
# Also define a sub_dir where this job has to be done. This will allow definitions of 
# opt, bond_scan, dimer, freq etc for the series of jobs in our work flow.
#
struct_list="Au38_2O_111_111_CH3_3"
poscar_name="POSCAR_latest"
sub_dir="opt"

CORES_PER_TASK="Not_used_on_archer"

export OMP_NUM_THREADS=1
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

echo Running on machine: $machine >  $LOG_FILE 
echo Running on host `hostname`   >> $LOG_FILE 
echo Time is `date`               >> $LOG_FILE  
echo Directory is `pwd`           >> $LOG_FILE  
echo Job ID is $JOBID             >> $LOG_FILE  
echo This jobs runs on the following machine: `echo $SLURM_JOB_NODELIST | uniq` >> $LOG_FILE 
echo Number of Processing Elements is $NCPUS >> $LOG_FILE 
echo Number of mpiprocs per node is $PPN >> $LOG_FILE 
echo >> $LOG_FILE 
echo VASP_PP_PATH set to $VASP_PP_PATH >> $LOG_FILE 
echo VASP Start Time is `date` running NCPUs=$NCPUS PPN=$PPN >> $LOG_FILE 
start="$(date +%s)"
#
# ... change to the scratch working directory       
cd $work_dir
#
# Loop over the job_list
#
for struct_dir in $struct_list; do
     echo job = $struct_dir >> $LOG_FILE 
     echo sub_dir = $sub_dir >> $LOG_FILE 
     export jobwork_dir="$work_dir"/"$struct_dir"/"$sub_dir"
     echo working in directory $jobwork_dir >> $LOG_FILE 
#
# Check if the sub_dir already exists
#
     if [ ! -d $jobwork_dir ]; then
       echo Creating working directory $jobwork_dir >> $LOG_FILE
       mkdir -p $jobwork_dir
     fi
# We are already in the work area, now move into the sub_dir   
     cd $jobwork_dir
#
     export full_poscar="$LAUNCH_DIR"/"$struct_dir"/"$poscar_name"
     if [ -e  "$full_poscar" ]; then
       echo Copying POSCAR file $full_poscar to  $jobwork_dir >> $LOG_FILE
       cp "$full_poscar" POSCAR
     else
      echo poscar "$full_poscar" file missing or mis-named >> $LOG_FILE
     fi 
#
# copy the script to work area
#
     if [ -e  "$LAUNCH_DIR"/"$ASE_SCRIPT" ]; then
       cp "$LAUNCH_DIR"/"$ASE_SCRIPT"  .
     else
       echo ase python script "$ASE_SCRIPT" file missing from "$LAUNCH_DIR" >> $LOG_FILE
     fi 
     

     echo $PWD Running python script $ASE_SCRIPT >> $LOG_FILE 
     echo  >> $LOG_FILE 
     echo With command: >> $LOG_FILE 
     echo python3 $ASE_SCRIPT $jobwork_dir $CORES_PER_TASK $JOBID $LAUNCH_DIR $struct_dir $sub_dir  \
                    $MACHINE > "$LAUNCH_DIR"/"$struct_dir"/"$sub_dir"/"$ASE_SCRIPT"_"$JOBID".out >> $LOG_FILE
#
#  run the ase script                                                             
#
     python3 $ASE_SCRIPT $jobwork_dir $CORES_PER_TASK $JOBID $LAUNCH_DIR $struct_dir $sub_dir $MACHINE \
                                    > "$LAUNCH_DIR"/"$struct_dir"/"$sub_dir"/"$ASE_SCRIPT"_"$JOBID".out &

     echo vasp run using ase script $ASE_SCRIPT for job $struct_dir running. >> $LOG_FILE 

     cd ..
   done

#
# Wait until all background jobs complete
#
wait

   for struct_dir in $job_list; do
  
     export jobwork_dir="$work_dir"/"$struct_dir"/"$sub_dir"
#     export jobwork_dir="$struct_dir"/work_"$JOBID"
     echo "copying back for " $jobwork_dir >> $LOG_FILE 
     cd $jobwork_dir

     cp OUTCAR  ${LAUNCH_DIR}/$struct_dir/OUTCAR_ase_"$JOBID"_$i"sub_dir"
     cp CONTCAR  ${LAUNCH_DIR}/$struct_dir/CONTCAR_ase_"$JOBID"_$i"sub_dir"
     cp INCAR  ${LAUNCH_DIR}/$struct_dir/INCAR_ase_"$JOBID"_$i"sub_dir"

     cd ..
#
done
#
# clean scratch
#

#rm -r $top_scratch 


