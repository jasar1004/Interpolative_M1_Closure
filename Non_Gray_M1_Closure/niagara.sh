#!/bin/bash                                                                                
## I want a simulation of 2 nodes                                                          
#SBATCH -N 2                                                                               
## Each node has 80 cpus. In total, I'll have 80x2=160 cpus                                
#SBATCH --ntasks-per-node=80                                                               
## Maximum allowed: 24 hours                                                               
#SBATCH -t 02:00:00                                                                        
## Job's name                                                                              
#SBATCH -J M1NONGRAYLCHI2
cd $SLURM_SUBMIT_DIR/

##Load here the modules you used to compile the executable. Here it was with ICC and MPICH2
module load NiaEnv/2019b intel/2019u4 intelmpi/2019u4

mpirun -np 100 ./Optimal_Mobius_Transformation > Output.log
