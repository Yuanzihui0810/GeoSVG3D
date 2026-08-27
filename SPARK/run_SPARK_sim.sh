#!/bin/bash
#SBATCH -J SPARK_sim
#SBATCH -p amd_512
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 10
#SBATCH --mem=50G
#SBATCH -o SPARK_sim_%j.out
#SBATCH -e SPARK_sim_%j.err

source activate zhyuanR44

arg1=$1
Rscript --vanilla /public3/home/scg5453/zhyuan/GeoSVG_SPARK_output/sim2/SPARK_sim.R $arg1
