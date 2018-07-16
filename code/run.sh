#!/bin/bash

export PETSC_DIR=/PATH/TO/PERSCDIR #i.e: export PETSC_DIR=/home/user/petsc-3.8.4
export PETSC_ARCH=arch-used #i.e: #export PETSC_ARCH=arch-linux2-c-debug
export SLEPC_DIR=/PATH/TO/PERSCDIR #i.e: export SLEPC_DIR=/home/user/slepc-3.8.3 


NUMPROCS=$1
#hostfile should be a file with the list of nodes in the cluster, i.e:
#node01
#node02
#...
#iface is the network interface that will be used to communicate with other nodes
time mpirun -n $NUMPROCS -iface eth0 -f ~/hostfile python3 shell_matrix_example.py -eps_nev 5 -log_view -eps_view  -eps_max_it 4000 -eps_tol 1e-14 -eps_ncv 11 -eps_monitor_all


