cmake_minimum_required(VERSION 3.13)
project(Particles_Simulation C)

set(CMAKE_C_STANDARD 11)

include_directories(.)

add_executable(Particles_Simulation
        simpar.c)

target_link_libraries(Particles_Simulation m)
