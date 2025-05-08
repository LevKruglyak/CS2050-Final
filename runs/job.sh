#!/bin/bash

ITERATIONS=50
SAVE_FREQ=10

mkdir -p runs-output

for i in 1 4 8; do
echo "RUN res1024-ppp50-$i"
mpirun -n $((i+1)) build/team_10_2025           \
    --config     runs/res1024-ppp50.json      \
    --iterations $ITERATIONS                  \
    --save-freq  $SAVE_FREQ                   \
    --out-prefix runs-output/res1024-ppp50-$i

echo "RUN res2048-ppp50-$i"
mpirun -n $((i+1)) build/team_10_2025           \
    --config     runs/res2048-ppp50.json      \
    --iterations $ITERATIONS                  \
    --save-freq  $SAVE_FREQ                   \
    --out-prefix runs-output/res2048-ppp50-$i

echo "RUN res4096-ppp50-$i"
mpirun -n $((i+1)) build/team_10_2025           \
    --config     runs/res4096-ppp50.json      \
    --iterations $ITERATIONS                  \
    --save-freq  $SAVE_FREQ                   \
    --out-prefix runs-output/4096-ppp50-$i
done
