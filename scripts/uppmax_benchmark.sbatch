#!/bin/bash
#SBATCH -A uppmax2025-2-316
#SBATCH -p pelle
#SBATCH -J nbody_pelle_1node
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH -t 00:30:00
#SBATCH -o nbody_pelle_1node_%j.out
#SBATCH -e nbody_pelle_1node_%j.err

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-$PWD}"

export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

bash "$PROJECT_ROOT/scripts/uppmax_job.sh"
