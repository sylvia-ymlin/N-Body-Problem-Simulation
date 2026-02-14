#!/bin/bash
set -e

# Setup
PROJECT_ROOT=$(pwd)
DATA_DIR=$PROJECT_ROOT/data
BUILD_DIR=$PROJECT_ROOT/build
LOG_FILE=$PROJECT_ROOT/full_benchmark_mac.log

# Ensure data dir exists
mkdir -p $DATA_DIR
cd $PROJECT_ROOT

# Helper function for timing (Mac compatible)
run_with_time() {
    local start_time=$(python3 -c 'import time; print(time.time())')
    "$@"
    local end_time=$(python3 -c 'import time; print(time.time())')
    local duration=$(echo "$end_time - $start_time" | bc)
    echo "Total time: ${duration}s"
}

# Clear log
echo "Starting Benchmark on $(hostname) at $(date)" > $LOG_FILE
echo "Hardware Info:" >> $LOG_FILE
sysctl -n machdep.cpu.brand_string >> $LOG_FILE
sysctl -n hw.memsize >> $LOG_FILE

# Generate input data
echo "Generating input data..." | tee -a $LOG_FILE
python3 scripts/generate_data.py 1000 $DATA_DIR/input_1k.gal
python3 scripts/generate_data.py 2000 $DATA_DIR/input_2k.gal
python3 scripts/generate_data.py 5000 $DATA_DIR/input_5k.gal
python3 scripts/generate_data.py 10000 $DATA_DIR/input_10k.gal
python3 scripts/generate_data.py 20000 $DATA_DIR/input_20k.gal
python3 scripts/generate_data.py 50000 $DATA_DIR/input_50k.gal

# Build
echo "Building project..." | tee -a $LOG_FILE
mkdir -p $BUILD_DIR
cd $BUILD_DIR
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(sysctl -n hw.logicalcpu)
cd $PROJECT_ROOT

# --- Phase I: v1 vs v2 (N=50k) ---
echo "--- Phase I: v1 vs v2 (N=50k) ---" | tee -a $LOG_FILE
echo "Running v1_naive (10 steps) -> Extrapolating to 100 steps..." | tee -a $LOG_FILE
# Run 10 steps, multiply by 10 for estimation
START=$(python3 -c 'import time; print(time.time())')
./build/v1_naive 50000 $DATA_DIR/input_50k.gal 10 0.001 1 0.5 1 2>&1 | tee -a $LOG_FILE
END=$(python3 -c 'import time; print(time.time())')
DURATION=$(echo "$END - $START" | bc)
ESTIMATED=$(echo "$DURATION * 10" | bc)
echo "v1_naive (10 steps) took ${DURATION}s. Estimated 100 steps: ${ESTIMATED}s" | tee -a $LOG_FILE

echo "Running v2_barnes_hut (100 steps)..." | tee -a $LOG_FILE
run_with_time ./build/v2_barnes_hut 50000 $DATA_DIR/input_50k.gal 100 0.001 1 0.5 1 2>&1 | tee -a $LOG_FILE

# --- Phase II: Arena Allocator (N=1k-50k) ---
echo "--- Phase II: Arena Allocator Impact (v2 vs v3_arena) ---" | tee -a $LOG_FILE
for N in 1000 2000 5000 10000 20000 50000; do
    INPUT_FILE=$DATA_DIR/input_${N/000/k}.gal
    if [ "$N" -eq "1000" ]; then INPUT_FILE="$DATA_DIR/input_1k.gal"; fi
    
    echo "N=$N: v2_barnes_hut (malloc)..." | tee -a $LOG_FILE
    ./build/v2_barnes_hut $N $INPUT_FILE 100 0.001 1 0.5 1 | grep "took" | tee -a $LOG_FILE
    
    echo "N=$N: v3_arena (arena)..." | tee -a $LOG_FILE
    ./build/v3_arena $N $INPUT_FILE 100 0.001 1 0.5 1 | grep "took" | tee -a $LOG_FILE
done

# --- Phase III: Morton Code & Optional Arena (v4) ---
echo "--- Phase III: Morton Code & Optional Arena (N=50k) ---" | tee -a $LOG_FILE
INPUT_FILE=$DATA_DIR/input_50k.gal

echo "v4_morton (use_arena=0)..." | tee -a $LOG_FILE
run_with_time ./build/v4_morton 50000 $INPUT_FILE 100 0.001 1 0.5 1 0 2>&1 | tee -a $LOG_FILE

echo "v4_morton (use_arena=1)..." | tee -a $LOG_FILE
run_with_time ./build/v4_morton 50000 $INPUT_FILE 100 0.001 1 0.5 1 1 2>&1 | tee -a $LOG_FILE

# --- Phase IV: Parallel Scaling & Optional Arena (v5) ---
echo "--- Phase IV: Parallel Scaling & Optional Arena (N=50k) ---" | tee -a $LOG_FILE
# Strong Scaling (Threads=1,2,4,8,12,16)
for T in 1 2 4 8 16; do
    echo "v5_parallel Threads=$T (use_arena=0)..." | tee -a $LOG_FILE
    ./build/v5_parallel 50000 $INPUT_FILE 100 0.001 $T 0.5 $T 0 | grep "Simulation took" | tee -a $LOG_FILE
    
    echo "v5_parallel Threads=$T (use_arena=1)..." | tee -a $LOG_FILE
    ./build/v5_parallel 50000 $INPUT_FILE 100 0.001 $T 0.5 $T 1 | grep "Simulation took" | tee -a $LOG_FILE
done

echo "Benchmark Complete." | tee -a $LOG_FILE
