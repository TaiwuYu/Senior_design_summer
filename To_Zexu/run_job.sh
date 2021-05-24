#PBS -A PAS0971 
#PBS -N void_test_T240
#PBS -j oe
#PBS -l nodes=1:ppn=1:gpus=1
#PBS -l walltime=2:00:00
#PBS -S /bin/bash
#PBS -m abe

cd $PBS_O_WORKDIR
#export OMP_NUM_THREADS=28
time ./Void_model.exe > output_test_T240
