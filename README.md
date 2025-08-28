# Python module to register image stacks of fluorescence microscopy and maldi images from the same sample

The main script is found in ./scripts/register.py.  This has a CLI documented below.  The workflow of the registration requires the selection of one MALDI channel and two channels of the fluorescence image,i.e. a "dense" marker, usually DAPI, and a "fiducial" marker,i.e. a more disperse and particle-like marker.  The MALDI channel should co-localize as much as possible with the fiducial marker in the fluorescence stack.
# Download container

- Singularity
```
singularity pull docker://ghcr.io/schapirolabor/femur-multi-omics:v1.5.0
```
## CLI
The main script is to be found in ./scripts/register.py
You can visualize in your terminal the CLI documentation using the following command
```bash
python register.py -h
```


### Required arguments
| Argument|Long name|Type|Description|Default value|
|---------|---------|----|-----------|-------------|
| -o | --output_dir | pathlib.Path | output directory where the registered maldi & mics image will be saved | NA |
| -inmi | --input_mics | pathlib.Path | Path to the macsima image stack (.tif) | NA |
|-mppmi|--pixel_size_microns_mics | float | microns per pixel (mpp) of the macsima image|NA|
|-idxmi|--reference_indices_mics | int  | Two integers separated by a space. The integers are the 0-based indices of the channels to be used in the mics stack for an initial coarse registration (usually DAPI) and a fine one respectively.  For the fine registration the following recommendations are given:(1) a marker that is sparse and (2) whose expression is highly correlated with the reference metabolite from the maldi image (see -idxma argument).|NA|
|-markmi|--marker_names_mics | pathlib.Path |.csv file with a column named "marker_name" containing the name of the channels of the mics stack (-inmi). The occurrence of the names in this column should match the order of occurrence in the mics stack. | NA |
|-inma|--input_maldi | pathlib.Path |Path to the maldi image stack (.tif)| NA |
|-mppma|--pixel_size_microns_maldi | float |microns per pixel (mpp) of the maldi image |NA|

|-qc|boolean flag | --qc_metrics | measure features of contrast, intensity and sharpness of each tile in the cycle and appends them to a table |FALSE|
|-wt|boolean flag | --write_table | writes a table with the acquisition parameters, metadata and,if enabled, qc metrics of each tile. Table will be saved in --output/cycle_info   |FALSE|

