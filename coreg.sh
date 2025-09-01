#Absolute paths to your local directories
input_dir=/home/local_data_disk/data_folder
output_dir=/home/local_data_disk/output_folder
tools_container=/home/local_data_disk/containers/femur-multi-omics_beta_v1.5.sif

#Paths of scripts inside the container
coreg_script=/tools/register.py

#Execution of register.py script

singularity exec --bind $input_dir:/mnt,$output_dir:/media  --no-home $tools_container python $coreg_script -o /media -inmi mnt/"mics.tif" -idxmi 0 4 -mppmi 0.170 -markmi /mnt/"markers.csv" -inma /mnt/"maldi_stack.tif" -idxma 10  -mppma 5 -markma /mnt/"maldi_markers.csv"