#!/bin/bash

make;
mpirun ./simpar-mpi 1 3 10 1;
mpirun ./simpar-mpi 1 3 1000000 20;
mpirun ./simpar-mpi 1 10 2000000 10;
mpirun ./simpar-mpi 1 30 20000000 10;