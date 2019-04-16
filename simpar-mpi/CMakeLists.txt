cmake_minimum_required(VERSION 3.13)
project(simpar_mpi C)

set(CMAKE_C_STANDARD 11)

add_executable(simpar_mpi
        simpar-mpi.c)