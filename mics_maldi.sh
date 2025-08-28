#!/bin/bash
#SBATCH --job-name=macsima
#SBATCH --partition=cpu-single
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=300gb
#SBATCH --array=1
module load system/singularity/3.11.3
#INPUTS
input_dir=/gpfs/bwfor/work/ws/hd_hi296-multiplexed_analysis/Soehnlein/coreg/beta_test/bspline_test/step_04

output_dir=/gpfs/bwfor/work/ws/hd_hi296-multiplexed_analysis/Soehnlein/coreg/beta_test/output/registration
dir1=/gpfs/bwfor/work/ws/hd_hi296-multiplexed_analysis/Soehnlein/coreg/beta_test/affine_transformation_results/pyramid
dir2=/gpfs/bwfor/work/ws/hd_hi296-multiplexed_analysis/Soehnlein/coreg/beta_test/data

mics_img="C-002_S-000_S_PE_R-01_W-A01_ROI-02_A-S100A9_C-Poly.tif"
maldi_img="maldi_neutros.tif"
#sample_array=/gpfs/bwfor/home/hd/hd_hd/hd_hi296/tsvs/mics_maldi.tsv
#sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $sample_array)

tools_container=/gpfs/bwfor/home/hd/hd_hd/hd_hi296/containers/femur-multi-omics_beta_v1.0.sif

singularity_file=/gpfs/bwfor/home/hd/hd_hd/hd_hi296/singularity.config


#END INPUTS

#PATHS OF SCRIPTS INSIDE CONTAINER

coreg_mics_maldi=/tools/scripts/coreg/mics_maldi/mics_maldi_reg.py
pyramidal_tiff=/tools/scripts/img_wr/pyramidal_tiff/pyramidal_tiff.py
append_pyramid=/tools/scripts/img_wr/append_pyramids/append_pyramids.py

#singularity exec --bind $dir1:/mnt,$dir2:/media --no-home $tools_container python $append_pyramid -i1 /mnt/$mics_img -i2 /mnt/$maldi_img -o /media  -n "mics_maldi_stack" -pl 8

singularity exec --bind $input_dir:/mnt,$output_dir:/media  --no-home $tools_container python $pyramidal_tiff -i /mnt -o /mnt/"pyramid" -n "pyramid_step_04" -res -rdim 38000 61314 -pl 8

#singularity exec --bind $input_dir:/mnt,$output_dir:/media --no-home $tools_container python $append_pyramid -i1 /mnt/$mics_img -i2 /mnt/$maldi_img -o /media  -n "mics_maldi_stack" -pl 8




