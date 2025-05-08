#!/bin/bash
#SBATCH --job-name=team-10-2025-huge
#SBATCH --nodes=12
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=huge-job_%j.out
#SBATCH --error=huge-job_%j.err

ITERATIONS=50
SAVE_FREQ=10
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

mkdir -p runs-output

echo "Job ID:           $SLURM_JOB_ID"
echo "Job name:         $SLURM_JOB_NAME"
echo "Running on nodes: $SLURM_JOB_NUM_NODES"
echo "Tasks per node:   $SLURM_NTASKS_PER_NODE"
echo "CPUs per task:    $SLURM_CPUS_PER_TASK"
echo "Total tasks:      $SLURM_NTASKS"
echo "------------------------------"
date
echo "Starting simulation..."

srun -N $SLURM_JOB_NUM_NODES build/team_10_2025   \
    --config     runs/res10000-ppp100.json      \
    --iterations $ITERATIONS                  \
    --save-freq  $SAVE_FREQ                   \
    --out-prefix runs-output/res10000-ppp100

echo "Simulation finished."
date
