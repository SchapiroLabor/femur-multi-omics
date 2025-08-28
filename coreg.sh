#!/bin/bash
#SBATCH --job-name=macsima
#SBATCH --partition=cpu-single
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=300gb
#SBATCH --array=1

module load system/singularity/3.11.3
input_dir=/gpfs/bwfor/work/ws/hd_hi296-multiplexed_analysis/Soehnlein/coreg/beta_test/data
output_dir=/gpfs/bwfor/work/ws/hd_hi296-multiplexed_analysis/Soehnlein/coreg/beta_test/data/output

tools_container=/gpfs/bwfor/home/hd/hd_hd/hd_hi296/containers/femur-multi-omics_beta_v1.4.sif
singularity_file=/gpfs/bwfor/home/hd/hd_hd/hd_hi296/singularity.config

#PATHS OF SCRIPTS INSIDE CONTAINER

coreg_script=/tools/register.py

singularity exec --bind $input_dir:/mnt,$output_dir:/media  --no-home $tools_container python $coreg_script -o /media -inmi mnt/"mics.tif" -idxmi 0 1 -mppmi 0.170 -inma /mnt/"maldi_stack.tif" -idxma 0  -mppma 5