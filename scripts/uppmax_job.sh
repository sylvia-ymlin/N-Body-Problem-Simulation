#!/bin/bash
set -euo pipefail

PROJECT_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BUILD_DIR="${BUILD_DIR:-$PROJECT_ROOT/build-uppmax}"
BENCHMARK_MODE="${BENCHMARK_MODE:-full}"
RESULT_FILE="${RESULT_FILE:-$PROJECT_ROOT/data/metrics/uppmax_${BENCHMARK_MODE}_results.json}"
THREADS="${THREADS:-1,2,4,8,16}"
N_BASE="${N_BASE:-100000}"
NSTEPS="${NSTEPS:-1}"
THETA="${THETA:-0.5}"

if [ -n "${SLURM_CPUS_PER_TASK:-}" ] && [ -z "${OMP_NUM_THREADS:-}" ]; then
  export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
fi

load_first_module() {
  local name
  for name in "$@"; do
    if module load "$name" >/dev/null 2>&1; then
      echo "Loaded module: $name"
      return 0
    fi
  done
  return 1
}

if command -v module >/dev/null 2>&1; then
  module purge >/dev/null 2>&1 || true
  load_first_module GCC/13.3.0 GCC/13.2.0 GCC/12.3.0 GCC/11.3.0 gcc >/dev/null 2>&1 || true
  load_first_module CMake/3.29.3 CMake/3.27.7 CMake/3.26.3 CMake/3.25.2 CMake >/dev/null 2>&1 || true
fi

CC_BIN="${CC:-$(command -v gcc || command -v cc || true)}"
CMAKE_BIN="${CMAKE:-$(command -v cmake || true)}"
PYTHON_BIN="${PYTHON:-$(command -v python3 || command -v python || true)}"

if [ -z "$CC_BIN" ] || [ -z "$CMAKE_BIN" ] || [ -z "$PYTHON_BIN" ]; then
  echo "Missing one of gcc/cc, cmake, or python3 on PATH."
  exit 1
fi

echo "Compiler: $CC_BIN"
echo "CMake:    $CMAKE_BIN"
echo "Python:   $PYTHON_BIN"
echo "Build dir: $BUILD_DIR"

"$CMAKE_BIN" -S "$PROJECT_ROOT" -B "$BUILD_DIR" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER="$CC_BIN"
"$CMAKE_BIN" --build "$BUILD_DIR" -j "${SLURM_CPUS_PER_TASK:-4}"

case "$BENCHMARK_MODE" in
  scaling)
    "$PYTHON_BIN" "$PROJECT_ROOT/scripts/benchmark_suite.py" \
      --binary "$BUILD_DIR/nbody_simulate" \
      --threads "$THREADS" \
      --n-base "$N_BASE" \
      --steps "$NSTEPS" \
      --theta "$THETA" \
      --output "$RESULT_FILE"
    ;;
  full)
    "$PYTHON_BIN" "$PROJECT_ROOT/scripts/benchmark_full_suite.py" \
      --binary "$BUILD_DIR/nbody_simulate" \
      --threads "$THREADS" \
      --n-base "$N_BASE" \
      --theta "$THETA" \
      --output "$RESULT_FILE"
    ;;
  ablation)
    "$PYTHON_BIN" "$PROJECT_ROOT/scripts/benchmark_ablation.py" \
      --binary "$BUILD_DIR/nbody_simulate" \
      --output "$RESULT_FILE"
    ;;
  all)
    "$PYTHON_BIN" "$PROJECT_ROOT/scripts/benchmark_full_suite.py" \
      --binary "$BUILD_DIR/nbody_simulate" \
      --threads "$THREADS" \
      --n-base "$N_BASE" \
      --theta "$THETA" \
      --output "$PROJECT_ROOT/data/metrics/uppmax_full_results.json"
    "$PYTHON_BIN" "$PROJECT_ROOT/scripts/benchmark_ablation.py" \
      --binary "$BUILD_DIR/nbody_simulate" \
      --output "$PROJECT_ROOT/data/metrics/uppmax_ablation_results.json"
    RESULT_FILE="$PROJECT_ROOT/data/metrics/uppmax_full_results.json and uppmax_ablation_results.json"
    ;;
  *)
    echo "Unknown BENCHMARK_MODE: $BENCHMARK_MODE"
    exit 1
    ;;
esac

echo "UPPMAX benchmark results written to $RESULT_FILE"
