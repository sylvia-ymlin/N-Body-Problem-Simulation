#!/bin/bash

# Unified Benchmark Suite
# Usage: ./scripts/unified_benchmark.sh <N> <steps>

N=${1:-5000}
STEPS=${2:-10}
DT=0.001
THREADS=8
THETA=0.5
K=8
INPUT="data/inputs/input_${N}.bin"

if [ ! -f "$INPUT" ]; then
    echo "Generating input data for N=$N..."
    mkdir -p data/inputs
    python3 scripts/generate_data.py $N $INPUT disk
fi

echo "--- N-Body Performance Journey: v1 to v4 ---"
echo "Configuration: N=$N, Steps=$STEPS, Threads=$THREADS"
echo "--------------------------------------------"

VERSIONS=("v1_naive" "v2_barnes_hut" "v2x_arena" "v3_morton" "v4_parallel")

for VER in "${VERSIONS[@]}"; do
    EXE="./build/$VER"
    if [ -f "$EXE" ]; then
        echo "Testing $VER..."
        $EXE $N $INPUT $STEPS $DT $THREADS $THETA $K
        echo ""
    else
        echo "Skip $VER (Binary not found: $EXE)"
    fi
done

echo "Benchmark Complete."
