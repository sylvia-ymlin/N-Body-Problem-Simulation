cmake_minimum_required(VERSION 3.10)
project(NBodySimulation C)

# Set C standard
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)

# Generate compile_commands.json for IDE IntelliSense
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Find OpenMP (Modern approach for macOS/Homebrew and others)
if(APPLE)
    # Check common Homebrew paths for libomp
    list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew/opt/libomp" "/usr/local/opt/libomp")
endif()

find_package(OpenMP REQUIRED)

# Add optimization flags
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    add_compile_options(-O3 -Wall -Wextra -march=native -ffast-math)
endif()

# ---------------------------------------------------------
# Define Versioned Executables for Benchmarking
# ---------------------------------------------------------

# v1: Naive O(N^2)
add_executable(v1_naive v1_naive/galsim.c)
target_link_libraries(v1_naive PRIVATE OpenMP::OpenMP_C m)

# v2: Barnes-Hut with malloc
add_executable(v2_barnes_hut 
    v2_barnes_hut/galsim.c 
    v2_barnes_hut/barnes_hut.c
)
target_include_directories(v2_barnes_hut PRIVATE v2_barnes_hut)
target_link_libraries(v2_barnes_hut m)

# v3: Barnes-Hut with Arena Allocator
add_executable(v3_arena 
    v3_arena/galsim.c 
    v3_arena/barnes_hut.c
)
target_include_directories(v3_arena PRIVATE v3_arena)
target_link_libraries(v3_arena m)

# v4: Barnes-Hut with Morton Code
add_executable(v4_morton 
    v4_morton/galsim.c 
    v4_morton/barnes_hut.c
    v4_morton/morton.c
)
target_include_directories(v4_morton PRIVATE v4_morton)
target_link_libraries(v4_morton m)

# v5: Final Parallel Version
add_executable(v5_parallel 
    v5_parallel/galsim.c 
    v5_parallel/barnes_hut.c
    v5_parallel/morton.c
    v5_parallel/kmeans.c
)
target_include_directories(v5_parallel PRIVATE v5_parallel)
target_link_libraries(v5_parallel PUBLIC OpenMP::OpenMP_C m)

# ---------------------------------------------------------
# End of Configuration
# ---------------------------------------------------------
