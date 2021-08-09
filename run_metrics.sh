#!/bin/sh

make clean
make teragen terametrics
mpirun -np 1 ./teragen -c 250000
mpirun -np 1 ./terametrics -c 20