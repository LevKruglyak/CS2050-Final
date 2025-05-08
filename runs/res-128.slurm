#!/bin/bash
#SBATCH --job-name=team-10-res-1024
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=01:00:00
#SBATCH --output=res-1024_%j.out
#SBATCH --error=res-1024_%j.err

ITERATIONS=100
SAVE_FREQ=10

echo "Job ID:           $SLURM_JOB_ID"
echo "Job name:         $SLURM_JOB_NAME"
echo "Running on nodes: $SLURM_JOB_NUM_NODES"
echo "Tasks per node:   $SLURM_NTASKS_PER_NODE"
echo "CPUs per task:    $SLURM_CPUS_PER_TASK"
echo "Total tasks:      $SLURM_NTASKS"
echo "------------------------------"
date
echo "Starting simulation..."

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun -N $SLURM_JOB_NUM_NODES -n $SLURM_NTASKS \
    build/team_10_2025                        \
    --config     runs/res-128.json            \
    --iterations $ITERATIONS                  \
    --save-freq  $SAVE_FREQ                   \
    --out-prefix runs-output/res-128

echo "Simulation finished."
date
