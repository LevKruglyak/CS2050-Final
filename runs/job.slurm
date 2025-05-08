#!/bin/bash
#SBATCH --job-name=team-10-2025
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err

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
echo "Starting simulations..."

for i in 1 4 8; do
echo "RUN res1024-ppp50-$i"
srun -n $((i+1)) build/team_10_2025   \
    --config     runs/res1024-ppp50.json      \
    --iterations $ITERATIONS                  \
    --save-freq  $SAVE_FREQ                   \
    --out-prefix runs-output/res1024-ppp50-$i

echo "RUN res2048-ppp50-$i"
srun -n $((i+1)) build/team_10_2025   \
    --config     runs/res2048-ppp50.json      \
    --iterations $ITERATIONS                  \
    --save-freq  $SAVE_FREQ                   \
    --out-prefix runs-output/res2048-ppp50-$i
done

# Too big for ram
for i in 4 8; do
echo "RUN res4096-ppp50-$i"
srun -n $((i+1)) build/team_10_2025   \
    --config     runs/res4096-ppp50.json      \
    --iterations $ITERATIONS                  \
    --save-freq  $SAVE_FREQ                   \
    --out-prefix runs-output/4096-ppp50-$i
done

echo "Simulations finished."
date
