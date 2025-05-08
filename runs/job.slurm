#!/bin/bash
#SBATCH --job-name=team-10-2025
#SBATCH --nodes=13
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --output=job.out
#SBATCH --error=job.err

ITERATIONS=25
SAVE_FREQ=0 # Don't output anything
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
echo "Building..."
cmake --build build
echo "Starting simulations..."

for i in 4 6 8 10 12; do
echo "RUN res1024-$i"
srun -N $((i+1)) build/lkxpm --config runs/res1024-ppp50.json --iterations $ITERATIONS --save-freq $SAVE_FREQ 
echo "RUN res2048-$i"
srun -N $((i+1)) build/lkxpm --config runs/res2048-ppp50.json --iterations $ITERATIONS --save-freq $SAVE_FREQ 
echo "RUN res4096-$i"
srun -N $((i+1)) build/lkxpm --config runs/res4096-ppp50.json --iterations $ITERATIONS --save-freq $SAVE_FREQ 
done

echo "Simulations finished."
date
