# Python scripts to register multi-channel fluorescence microscopy and maldi images

The main script is found in ./scripts/register.py.  This has a CLI documented below.  The workflow of the registration requires the selection of one MALDI channel and two channels of the fluorescence image,i.e. a "dense" marker, usually DAPI, and a "fiducial" marker,i.e. a more disperse and particle-like marker.  The MALDI channel should co-localize as much as possible with the fiducial marker in the fluorescence stack.
The output of this tool is an ome.tif file with the registered MACsima and MALDI images.  The MALDI images in this file are upscaled to the same pixel size of the MACsima data.

# Usage
## via python
- (1) Download the [./scripts](https://github.com/SchapiroLabor/femur-multi-omics/tree/main/scripts) folder
- (2) Download the [environment.yml](https://github.com/SchapiroLabor/femur-multi-omics/blob/main/environment.yml) file
- (3) Create conda environment:
```
conda env create -f *local_directory*/environment.yml -n *custom_environment_name*
```
- (4) Activate environment:
```
conda activate *custom_environment_name*
```
- (5) Check the CLI of the register.py script by executing:
```
python *local_directory*/register.py --help
```
or check the CLI documentation below.

- (6) Execute register.py with your correspondent entries for the CLI, for example:
```
python *local_directory*/register.py -o *output_directory*  -inmi *macsima_stack.tif* \
-mppmi *0.170* -idxmi *0* *4* -markmi *markers.csv* -inma *maldi_stack.tif* \
-mppma *5* -idxma *10* -markma *maldi_markers.csv* 
```
In the example above we suppose the dense marker in the MACSima images is in the channel with index 0 and the fiducial marker is in channel 4.  The latter marker is associated to channel 10 of the maldi stack.

## via container

**Singularity**
- (1) Download container:

```
singularity pull docker://ghcr.io/schapirolabor/femur-multi-omics:v1.5.0
```
This will download a container named femur-multi-omics:v1.5.0.sif

- (2) Execute the register.py script inside the /tools folder in the .sif container.  Be aware of two things in the example below,(a) the variables mnt, media and tools are names of folders inside the container and need no editing in the example below. (b) We suppose that the macsima and maldi stacks and associated .csv information is in a local directory *input_directory*, while the output of the script will be written in a different local directory *output_directory*.

Example:

```
singularity exec --bind *input_directory*:/mnt,*output_directory*:/media  --no-home femur-multi-omics:v1.5.0.sif python /tools/register.py -o /media  -inmi mnt/*macsima_stack.tif* \
-mppmi *0.170* -idxmi *0* *4* -markmi mnt/*markers.csv* -inma mnt/*maldi_stack.tif* \
-mppma *5* -idxma *10* -markma mnt/*maldi_markers.csv* 
```

See also an example of the implementation above inside a bash script in [coreg.sh](https://github.com/SchapiroLabor/femur-multi-omics/blob/main/coreg.sh).

## CLI
The main script is ./scripts/register.py
You can visualize in your terminal the CLI documentation using the following command
```
python register.py --help
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
|-idxma|--reference_indices_maldi | int |0-based index of the reference channel in the maldi stack |NA|
|-markma|marker_names_maldi | pathlib.Path | .txt file with a column named "m/z" containing the mass-to-charge ratio of the maldi stack (-inma).The occurrence of the ratio in this column should match the order of occurrence in the maldi stack. |NA|

