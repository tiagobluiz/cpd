cmake_minimum_required(VERSION 3.3)
project(simpar-mpi)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpicxx)

set(SOURCE_FILES simpar-mpi.c testing-mpi.c)
add_executable(simpar-mpi ${SOURCE_FILES})