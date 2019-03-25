cmake_minimum_required(VERSION 3.13)
project(simpar_omp_1 C)

set(CMAKE_C_STANDARD 11)

add_executable(simpar_omp_1
        simpar-omp.1.c)

# added -fopenmp
set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -D_DEBUG")

set(SOURCE_FILES simpar-omp.1.c)
add_executable(openmp_test_clion ${SOURCE_FILES})
