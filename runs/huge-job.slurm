#!/bin/bash
#SBATCH --job-name=team-10-2025-huge
#SBATCH --nodes=12
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --output=huge-job.out
#SBATCH --error=huge-job.err

ITERATIONS=500
SAVE_FREQ=50
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
    --config     runs/res8192-ppp50.json      \
    --iterations $ITERATIONS                  \
    --save-freq  $SAVE_FREQ                   \
    --out-prefix runs-output/res8192-ppp50

echo "Simulation finished."
date
