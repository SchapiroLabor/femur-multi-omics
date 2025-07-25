#public libraries
from pathlib import Path
from skimage.filters import threshold_otsu
from skimage import transform
import tifffile as tifff
import numpy as np
from skimage.transform import resize
#local libraries
import CLI
import register_funcs


def is_pyramid(img_path):
    """
    Checks if image has pyramidal levels

    """
    
    with tifff.TiffFile(img_path) as tif:
        pyramid=len(tif.series[0].levels) > 1
    return pyramid



def extract_and_resize(tiff_stack_path,ch_index,mpp_src=None,mpp_target=None):
    """
    This function extracts a channel from a pyramidal tifff_stack and resizes the image
    as to match the target microns_per_pixel

    Args:
        list_of_dicts (list): list of dictionaries with common keys
    Returns:
        merged_dict (dict): dictionary with the values stored in lists
    """
    if is_pyramid(tiff_stack_path):
        with tifff.TiffFile(tiff_stack_path) as tif:
            if len(tif.pages)>1:
                coll_xy_dim=[lev.shape[1:3] for lev in tif.series[0].levels]
            else:
                coll_xy_dim=[lev.shape for lev in tif.series[0].levels]
    

        resize_factor=mpp_src/mpp_target
        target_dim=np.rint( resize_factor*np.array(coll_xy_dim[0]) )
        target_dim=[int(element) for element in target_dim]
        nearest_lvl_index=np.argmin( np.abs( [target_dim[0]-element[0] for element in coll_xy_dim] ) )

        img_aux=tifff.imread(tiff_stack_path,series=0,key=ch_index,level=nearest_lvl_index)
        img_type=img_aux.dtype.name
        #output_img=resize(img_aux.asarray(),output_shape=target_dim,order=1,preserve_range=True).astype(img_type)
        output_img=resize(img_aux,output_shape=target_dim,order=1,preserve_range=True).astype(img_type)

    return output_img





def main():

    args = CLI.get_args()
    mpp_mics=args.pixel_size_microns_mics
    mpp_maldi=args.pixel_size_microns_maldi
    chan_mics=args.reference_indices_mics
    chan_maldi=args.reference_indices_maldi

    mics_resized=extract_and_resize(args.input_mics,
                               chan_mics,
                               mpp_mics,
                               mpp_maldi
                               )
    #tifff.imwrite(Path(r"C:\Users\VictorP\repos\femur-multi-omics\output") / "resized_mics.tif",mics_resized,photometric="minisblack")
 
    maldi_ref_img=tifff.imread(args.input_maldi,key=chan_maldi)
    maldi_registered,trfm_log=register_funcs.register_arrays(mics_resized,maldi_ref_img,mpp_maldi)
    maldi_stack=tifff.imread(args.input_maldi)
    register_funcs.register_from_params(maldi_stack,trfm_log,mpp_maldi,args.output_dir)
 


if __name__ == '__main__':
    main()









