#!/bin/bash
PROJECT_ROOT=$(pwd)
BUILD_DIR=$PROJECT_ROOT/build
DATA_DIR=$PROJECT_ROOT/data

echo "=== Remote Server Benchmark (1-Step Scaling) ==="
cd $BUILD_DIR
for t in 1 2 4 8 16; do
    echo "----------------------------------------"
    echo "Running with $t threads..."
    # Using 1 step for the Million Particle benchmark to get fast results
    ./nbody_simulate 5 1000000 ../data/input_1M.gal 1 0.001 $t 0.5 $t 1
done
echo "----------------------------------------"
echo "Benchmark Complete."
