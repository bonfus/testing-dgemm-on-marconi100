#!/bin/bash
#SBATCH --nodes=2               # number of nodes
#SBATCH --ntasks-per-node=4     # number of tasks per node
#SBATCH --ntasks-per-socket=2   # number of tasks per socket
#SBATCH --cpus-per-task=32      # number of HW threads per task (equal to OMP_NUM_THREADS*4)
#SBATCH --mem=230000MB
#SBATCH --time 01:00:00         # format: HH:MM:SS
#SBATCH -A cin_preM100 
#SBATCH -p system

cat run_me.sh

module load profile/global pgi/19.10--binary gnu/8.4.0 fftw/3.3.8--gnu--8.4.0 cuda/10.1


# PGI + Internal Blas
mpif90 -mp test.f90 -lblas -o pgi_blas.x

for i in {1..16} {18..32..2} {36..128..4}
do
    export OMP_NUM_THREADS=$i
    mpirun --bind-to none  ./pgi_blas.x 2000 2000  | tr '\n' '\t'
    echo
done 

echo '=== OpenBLAS ==='
mpif90 -mp test.f90 -o pgi_openblas.x -L./openblas/lib -lopenblas
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`/openblas/lib

for i in {1..16} {18..32..2} {36..128..4}
do
    export OMP_NUM_THREADS=$i
    mpirun --bind-to none ./pgi_openblas.x 2000 2000  | tr '\n' '\t'
    echo
done 

