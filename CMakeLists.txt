cmake_minimum_required(VERSION 3.10)
project(NBodySimulation C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(OpenMP)
if(OpenMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_C_FLAGS}")
endif()

# Core Sources
set(CORE_SOURCES
    src/core/io.c
)

# Utility Sources
set(UTILS_SOURCES
    src/utils/morton.c
    src/utils/kmeans.c
)

# Kernel Sources
set(KERNEL_SOURCES
    src/kernels/v1_naive.c
    src/kernels/v2_barnes_hut.c
    src/kernels/v3_arena.c
    src/kernels/v4_morton.c
    src/kernels/v5_parallel.c
)

add_library(core_lib ${CORE_SOURCES} ${UTILS_SOURCES} ${KERNEL_SOURCES})
if(OpenMP_FOUND)
    target_link_libraries(core_lib OpenMP::OpenMP_C)
endif()

add_executable(nbody_simulate src/core/main.c)
target_link_libraries(nbody_simulate PRIVATE core_lib m)

if(OpenMP_FOUND)
    target_link_libraries(nbody_simulate PRIVATE OpenMP::OpenMP_C)
endif()
