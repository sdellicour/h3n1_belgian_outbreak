#!/bin/bash
#SBATCH --job-name=Analysis_2_j_1000
#SBATCH --time=235:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem-per-cpu=10240
#SBATCH --qos=gpu_prio

module load beagle-lib/3.0.2-fosscuda-2018b

cd
cd H3N1
java -jar beast_1105_290920.jar -beagle_gpu -beagle_double -beagle_order 1 -overwrite Analysis_2_j_1000.xml
