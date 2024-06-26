#!/bin/bash --login
#SBATCH -o CsPbI3.o           # Job output file
#SBATCH -e CsPbI3.e           # Job error file
#SBATCH -J CsPbI3_100         # Job name
#SBATCH -p dev                         # Queue partition
#SBATCH --ntasks=40                    # number of parallel processes (tasks)
#SBATCH --ntasks-per-node=40           # tasks to run per node
#SBATCH --time=00:50:00                # time limit
#SBATCH --exclusive                    # exclusive node access
#SBATCH -A scw1161                     # Project code

## Usage
# sbatch script # submit job
# squeue        # job status
# 
echo "Starting script...."
#
export CODE="VASP"
MACHINE='hawk'
ASE_SCRIPT='ase_vasp_slab_opt.py'
struct_list="CsPbI3"    
CORES_PER_TASK=40
#
struct_name="POSCAR_CsPbI3"       
sub_dir="opt"
zfix=40
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
# To use ase will need python
     module load ase/3.20.1 python/3.7.0 pymatgen/2018.11.30

     JOBID=$SLURM_JOBID
     LAUNCH_DIR=$PWD
# Record machine set up info
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
     LAUNCH_DIR=$PWD
     JOBID=$SLURM_JOBID
#
     module load vasp/5
     module load cray-python
# Make sure we pick up local installation of ase
     export PYTHONUSERBASE=/work/e05/e05/$USER/.local
     export PATH=$PYTHONUSERBASE/bin:$PATH
     export PYTHONPATH=$PYTHONUSERBASE/lib/python3.8/site-packages:$PYTHONPATH

     CORES_PER_TASK="Not_used_on_archer2"
#
# On Archer2 PAW are availble at:
# $VASP_PSPOT_DIR make a link so that ase can find these
#
     export VASP_PP_PATH=$LAUNCH_DIR/pot_link
#
     if [ -e $VASP_PP_PATH ];then
        echo link for VASP_PP_PATH directory already in place
     else
        mkdir -p $VASP_PP_PATH
#
        export ARCHER2_PP_PBE="$VASP_PSPOT_DIR"/PBE
        export    LINK_PP_PBE="$VASP_PP_PATH"/potpaw_PBE
        ln -s $ARCHER2_PP_PBE $LINK_PP_PBE
     fi
#     
   elif [[ $MACHINE == "young" ]]; then
#     
      module unload compilers mpi
      module load compilers/intel/2017/update1
      module load mpi/intel/2017/update1/intel
      module load vasp/5.4.4-18apr2017/intel-2017-update1
# To use ase will need python
      module load python3/recommended
#
# Locate PAWs
#
      export VASP_PP_PATH=$HOME/progs/vasp
#
      LAUNCH_DIR=$PWD
      JOBID=$JOB_ID
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
#
# Do the module loads
#
      module purge
# load the VASP module, it loads required compiler, mpi and mkl libraries
      module load  vasp/5.4.4
# To use ase will need python
      module load ase/3.20.1 python/3.7.0 pymatgen/2018.11.30
#
# Record machine set up info
     NNODES=$SLURM_NNODES
     NCPUS=$SLURM_NTASKS
     PPN=$SLURM_NTASKS_PER_NODE
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
#    rm -rf $work_dir
    mkdir -p $work_dir
# On Young CASTEP also needs a tmp directory for this run
    if [[ $CODE == "CASTEP" ]]; then
       export tmp_dir=$work_dir/tmp
       mkdir -p $tmp_dir
       export FORT_TMPDIR=$tmp_dir
    fi 
#
elif [[ $MACHINE == "archer2" ]]; then
    export work_dir=$LAUNCH_DIR/workspace/$JOBID
#    rm -rf $work_dir
    mkdir -p $work_dir
# On ARCHER2 CASTEP also needs a tmp directory for this run
    if [[ $CODE == "CASTEP" ]]; then
       export tmp_dir=$work_dir/tmp
       mkdir -p $tmp_dir
       export FORT_TMPDIR=$tmp_dir
    fi 
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
echo Start Time is `date`  >> $LOG_FILE 
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
     export home_sub_dir=${LAUNCH_DIR}/$struct_dir/"$sub_dir"
     echo "copying back for : " $jobwork_dir >> $LOG_FILE
     echo "to home sub dir  : " $home_sub_dir >> $LOG_FILE
     cd $jobwork_dir

     cp *.py $home_sub_dir
#
# Only copy particular files from the CASTEP directory
#
     if [[ $CODE == "CASTEP" ]]; then
#
# Remove the "Processing" lines from the out file
#
        awk '$1 !~ /Proc/ {print}' "$home_sub_dir"/"$ASE_SCRIPT"_"$JOBID".out > temp.out
        mv temp.out "$home_sub_dir"/"$ASE_SCRIPT"_"$JOBID".out
#
        mkdir $home_sub_dir/castep_"$JOBID"
#
        cp CASTEP/castep.castep $home_sub_dir/castep_"$JOBID"
        cp CASTEP/castep.bands  $home_sub_dir/castep_"$JOBID"
        cp CASTEP/castep.cell   $home_sub_dir/castep_"$JOBID"
        cp CASTEP/castep.params $home_sub_dir/castep_"$JOBID"
#
     elif [[ $CODE == "VASP" ]]; then
        cp OUTCAR  "$home_sub_dir"/OUTCAR_"$JOBID"
        cp CONTCAR "$home_sub_dir"/CONTCAR_"$JOBID"
        cp INCAR   "$home_sub_dir"/INCAR_"$JOBID"
        cp POS*    "$home_sub_dir"
# Files for occasional jobs
        if [ -e LOCPOT ]; then
           cp LOCPOT  "$results_dir"/LOCPOT_"$JOBID"
        fi
        if [ -e DOSCAR ]; then
           cp DOSCAR  "$results_dir"/DOSCAR_"$JOBID"
        fi
        if  [ -e *".csv" ]; then
          cp *".csv" "$home_sub_dir"
        fi
#
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
