#!/bin/bash -l
#SBATCH -o ./DAMASKPDT.Make.STDOUT.%j
#SBATCH -e ./DAMASKPDT.Make.STDERR.%j
#SBATCH -D ./
#SBATCH -J damaskpdt
#SBATCH --partition=p.talos
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=m.kuehbach@mpie.de
#SBATCH --time=10:00:00

module load cmake
module load intel
module load impi
unset I_MPI_HYDRA_BOOTSTRAP
unset I_MPI_PMI_LIBRARY
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


#-DCMAKE_C_COMPILER=icc ..

srun cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=icpc .. 
srun make

