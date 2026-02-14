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
    add_compile_options(-O3 -Wall -Wextra)
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

# v2x: Barnes-Hut with Arena Allocator (Experimental)
add_executable(v2x_arena 
    v2x_arena_experimental/galsim.c 
    v2x_arena_experimental/barnes_hut.c
)
target_include_directories(v2x_arena PRIVATE v2x_arena_experimental)
target_link_libraries(v2x_arena m)

# v3: Barnes-Hut with Arena + Morton Code
add_executable(v3_morton 
    v3_morton/galsim.c 
    v3_morton/barnes_hut.c
    v3_morton/morton.c
)
target_include_directories(v3_morton PRIVATE v3_morton)
target_link_libraries(v3_morton m)

# v4: Final Parallel Version
add_executable(v4_parallel 
    v4_parallel/galsim.c 
    v4_parallel/barnes_hut.c
    v4_parallel/morton.c
    v4_parallel/kmeans.c
)
target_include_directories(v4_parallel PRIVATE v4_parallel)
target_link_libraries(v4_parallel PUBLIC OpenMP::OpenMP_C m)

# ---------------------------------------------------------
# End of Configuration
# ---------------------------------------------------------
