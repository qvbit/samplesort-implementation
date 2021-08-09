#!/bin/sh

make clean
make teragen teravalidate
make terasort
mpirun -np 8 ./teragen -c 250000
mpirun -np 8 ./terasort
mpirun -np 8 ./teravalidate