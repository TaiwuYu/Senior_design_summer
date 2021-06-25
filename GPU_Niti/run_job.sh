#PBS -A PAS0971 
#PBS -N NiTi_rand_515
#PBS -j oe
#PBS -l nodes=1:ppn=3:gpus=1
#PBS -l walltime=2:00:00
#PBS -S /bin/bash
#PBS -m abe

cd $PBS_O_WORKDIR
#export OMP_NUM_THREADS=28
time ./HEA_model_v9_update.exe > output_rand_515
