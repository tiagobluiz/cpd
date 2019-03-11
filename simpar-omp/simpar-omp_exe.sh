#!/bin/bash

make;
./simpar-omp 1 3 10 1;
./simpar-omp 1 3 1000000 20;
./simpar-omp 1 10 2000000 10;
./simpar-omp 1 30 20000000 10;
