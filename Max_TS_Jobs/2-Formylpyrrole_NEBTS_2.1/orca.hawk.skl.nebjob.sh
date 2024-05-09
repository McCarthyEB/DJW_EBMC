#!/bin/bash
#SBATCH -e %x.e.%j
#SBATCH -o %x.o.%j
#SBATCH -N 2
#SBATCH --tasks-per-node=40
#SBATCH --time=72:00:00
#SBATCH -p compute
#SBATCH --account=scw1039
#SBATCH --exclusive

set -eu
# Jobname below is set automatically when submitting like this: sbatch orca.hawk.skl.h2o.sh
#Can alternatively be set manually below. job variable should be the name of the inputfile without extension (.inp)
job=${SLURM_JOB_NAME}
job=$(echo ${job%%.*})

module purge
#Set OPENMPI paths here:
module load compiler/gnu/9.2.0
module load orca/5.0.0
module list

# Creating scratch folder for the user
WDPATH=/scratch/$USER/orca.$SLURM_JOB_ID
rm -rf $WDPATH
mkdir -p $WDPATH
 
NNODES=$SLURM_NNODES
NCPUS=$SLURM_NTASKS
PPN=$SLURM_NTASKS_PER_NODE

# Copy only the necessary stuff in submit directory to scratch directory. Add more here if needed.
CASE=nebjob
LOG=ORCA.SKL.$job.out.$NCPUS.$SLURM_JOB_ID

cp  $SLURM_SUBMIT_DIR/$CASE.inp $WDPATH/$job.inp
cp  $SLURM_SUBMIT_DIR/prod.xyz $WDPATH/prod.xyz
cp  $SLURM_SUBMIT_DIR/react.xyz $WDPATH/react.xyz

# Creating nodefile in scratch
scontrol show hostnames "$SLURM_JOB_NODELIST" | sort | uniq > $WDPATH/$job.nodes

# cd to scratch
cd $WDPATH

# Copy job and node info to beginning of outputfile
echo "Job execution start: $(date)" >>  $SLURM_SUBMIT_DIR/$LOG
echo "Shared library path: $LD_LIBRARY_PATH" >>  $SLURM_SUBMIT_DIR/$LOG
echo "Slurm Job ID is: ${SLURM_JOB_ID}" >>  $SLURM_SUBMIT_DIR/$LOG
echo "Slurm Job name is: ${SLURM_JOB_NAME}" >>  $SLURM_SUBMIT_DIR/$LOG
echo $SLURM_NODELIST >> $SLURM_SUBMIT_DIR/$LOG

#Start ORCA job. ORCA is started using full pathname (necessary for parallel execution).
#Output file is written directly to submit directory on frontnode.
start="$(date +%s)"
/apps/chemistry/orca/5.0.0/el7/orca $WDPATH/$job.inp >>  $SLURM_SUBMIT_DIR/$LOG
stop="$(date +%s)"
finish=$(( $stop-$start ))
echo ORCA case = $CASE $SLURM_JOBID  nproc = $NCPUS ppn = $PPN Job-Time $finish seconds
 
# ORCA has finished here. Now copy important stuff back (xyz files, GBW files etc.).
# Add or remove items as needed.
cp $WDPATH/* $SLURM_SUBMIT_DIR
