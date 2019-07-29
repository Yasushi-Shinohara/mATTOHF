#!/bin/sh
ulimit -s unlimited
export OMP_STACKSIZE=256G
export OMP_NUM_THREADS=1
#mpirun -np 1 ../src/mATTOHF <<<'"./input.dat"' >log.log
#mpirun -np 16 ../src/mATTOHF <<<'"./input.dat"' >log.log
mpirun -np 2 ../src/mATTOHF <<< '"./input.dat"' > log.log
