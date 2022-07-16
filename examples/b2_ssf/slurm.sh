#!/bin/bash
#SBATCH --export=ALL
#SBATCH --job-name=TEST
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -t 24:00:00
#SBATCH --requeue

echo "===================================================="
echo " submitted to local node: " $HOSTNAME
echo "===================================================="


module purge
module load openmpi/intel-13.0/1.6.3/64
module load fftw
module load intel/13.0/64/13.0.1.117
module load intel-mkl/11.0/1/64

VASP=/home/EAC/VASP/vasp.5.3.5
PROFESS=/home/EAC/PROFESS/pPROFESS_AMD_2722
ABINIT=/home/EAC/abinit/abinit-7.0.5/bins/abinit
ABINIT_BLPS=/home/EAC/BLPS_gen/InvertKS_abinit/abinit
PROFESS_TEST=/home/mohanc/7_Important/PROFESS/trunk/OFDFT/pOFDFT
CORRELATION=./curr_corr_functs.exe
SPLIT=./splitXCAR.exe
D310=/home/mohanc/3_Mohan_Codes/D310/src/D310.exe

cd $SLURM_SUBMIT_DIR

$D310 > OUTPUT

