#!/bin/bash
#SBATCH -J SPARK_real
#SBATCH -p amd_512
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 10
#SBATCH --mem=20G
#SBATCH -o SPARK_real_%j.out
#SBATCH -e SPARK_real_%j.err

source activate zhyuanR44

Rscript --vanilla /public3/home/scg5453/zhyuan/GeoSVG_SPARK_output/SPARK_real_server.R