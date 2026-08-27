#!/bin/bash

for seed in {1..24}
do
    echo "sbatch run_SPARK_sim.sh $seed"
    sbatch run_SPARK_sim.sh $seed
done