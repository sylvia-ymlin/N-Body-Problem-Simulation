cmake_minimum_required(VERSION 3.10)
project(NBodySimulation C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(OpenMP)
if(OpenMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_C_FLAGS}")
endif()

set(SOURCES
    io.c
    morton.c
    kmeans.c
    naive.c
    barnes_hut.c
)

add_library(core_lib ${SOURCES})
if(OpenMP_FOUND)
    target_link_libraries(core_lib OpenMP::OpenMP_C)
endif()

add_executable(nbody_simulate main.c)
target_link_libraries(nbody_simulate PRIVATE core_lib m)

if(OpenMP_FOUND)
    target_link_libraries(nbody_simulate PRIVATE OpenMP::OpenMP_C)
endif()
