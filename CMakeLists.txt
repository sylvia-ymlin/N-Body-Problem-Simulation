cmake_minimum_required(VERSION 3.10)
project(NBodySimulation C)

# Set C standard
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)

# Find OpenMP (Specific Logic for macOS/Homebrew)
if(APPLE)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Xpreprocessor -fopenmp")
    set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp")
    set(OpenMP_C_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY "/opt/homebrew/opt/libomp/lib/libomp.dylib")
    
    # Include directories for Homebrew libomp
    include_directories("/opt/homebrew/opt/libomp/include")
    link_directories("/opt/homebrew/opt/libomp/lib")
endif()

find_package(OpenMP REQUIRED)

# Add optimization flags
if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_C_COMPILER_ID MATCHES "Clang")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wall -Wextra")
endif()

# Define the executable
add_executable(simulate
    src/galsim.c
    src/barnes_hut.c
    src/kmeans.c
    src/morton.c
)

# Link OpenMP and Math library
target_link_libraries(simulate PUBLIC OpenMP::OpenMP_C m)

# Include directories
target_include_directories(simulate PUBLIC src)
