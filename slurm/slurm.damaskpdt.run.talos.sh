#!/bin/bash -l
#SBATCH -o ./DAMASKPDT.STDOUT.%j
#SBATCH -e ./DAMASKPDT.STDERR.%j
#SBATCH -D ./
#SBATCH -J damaskpdt
#SBATCH --partition=p.talos
#SBATCH --nodes=28
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mail-user=m.kuehbach@mpie.de
#SBATCH --time=24:00:00
#SBATCH --exclusive

export MYDAMASKFN=DAMASKPDT.Paper13.Bambach2019.DISTEX.xml
echo $MYDAMASKFN

##name scheme details simulation set, i) mpi processes, ii) threads, iii) trial id
export FLOWCU=1
export DISTEX=3
export DISLOU=4
export SDFTEX=5
export SDFLOU=6
export VORTEX=7
export VORLOU=8

export MP=28
export OM=40
export TRIAL=1

export SIMID=$DISTEX$MP$OM$TRIAL

echo $SIMID
mkdir -p $SIMID


module load cmake
module load intel
module load impi
##unset I_MPI_HYDRA_BOOTSTRAP
##unset I_MPI_PMI_LIBRARY
module load mkl
module load boost
module list

if [ ! -z $SLURM_CPUS_PER_TASK ] ; then
	export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
else
	export OMP_NUM_THREADS=1
fi
export MKL_NUM_THREADS=1
export OMP_PLACES=cores

echo $OMP_NUM_THREADS
echo $MKL_NUM_THREADS
echo $OMP_PLACES


export LIBRARY_PATH=/mpcdf/soft/SLE_15/packages/x86_64/intel_parallel_studio/2018.4/mkl/lib/intel64_lin/
echo $LIBRARY_PATH

ulimit -s unlimited
ulimit


##list allocated TALOS resources nodes
echo $SLURM_JOB_NODELIST
##echo $SLURM_NODELIST
echo $SLURM_JOB_NUM_NODES
##echo $SLURM_NNODES

echo "job $SLURM_JOB_NAME with job id $SLURM_JOB_ID is running on $SLURM_JOB_NUM_NODES node(s): $SLURM_JOB_NODELIST"

srun damaskpdt_vorocomposer_final $SIMID $MYDAMASKFN 500g256_256_256_compressionZ

mv *.$SIMID.* $SIMID/
