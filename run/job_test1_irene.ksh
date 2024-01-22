#!/bin/bash

#MSUB -r  qgw                 # Request name
#MSUB -n  2                   # Number of tasks to run
#MSUB -N  1                   # Number of nodes to use
#MSUB -T 3600                 # Elapsed time limit in seconds
#MSUB -o qgw.o%I              # Standard output. %I is the job id
#MSUB -e qgw.e%I              # Error output. %I is the job id
#MSUB -q rome                 # partition  name
#MSUB -A gen12020             # Project ID for accounting
#MSUB -x                      # exclusve compute node
#MSUB -m store,work,scratch   # list the file system used by the job

ulimit -s unlimited

NPROC=2

mpirun -n $NPROC qg.e
