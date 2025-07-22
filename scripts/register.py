#public libraries
from pathlib import Path
from skimage import transform
from skimage.transform import pyramid_gaussian
import tifffile as tifff
import numpy as np
from skimage.transform import resize
from skimage.util import img_as_float32
from skimage.exposure import rescale_intensity
import itk
import os
#local libraries
import CLI


#HELPERS
def is_pyramid(img_path):
    """
    Checks if image has pyramidal levels

    """
    
    with tifff.TiffFile(img_path) as tif:
        pyramid=len(tif.series[0].levels) > 1
    return pyramid

def get_moment_features(img_arr):
    features={}
    moments = itk.ImageMomentsCalculator.New(itk.GetImageFromArray(img_arr))
    moments.Compute()
    principal_axis=itk.array_from_matrix(moments.GetPrincipalAxes())
    minor_axis=principal_axis[0]
    major_axis=principal_axis[1]

    x_min=minor_axis[0]
    y_min=minor_axis[1]
    x_max=major_axis[0]
    y_max=major_axis[1]
    

    features["centroid"]=tuple(moments.GetCenterOfGravity())
    features["major_ax"]=major_axis
    features["minor_ax"]=minor_axis
    features["major_ax_or"]=np.arctan(x_max/y_max)#[radians]relative to y+ axis clockwise (/) negative, counterclowise (\) positive
    features["minor_ax_or"]=np.arctan(x_min/y_min)
    return features



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
                #collection of xy_dimensions of each level in tiff_stack
                coll_xy_dim=[lev.shape[1:3] for lev in tif.series[0].levels]
            else:
                coll_xy_dim=[lev.shape for lev in tif.series[0].levels]

        resize_factor=mpp_src/mpp_target
        target_dim=np.rint( resize_factor*np.array(coll_xy_dim[0]) )
        target_dim=[int(element) for element in target_dim]
        nearest_lvl_index=np.argmin( np.abs( [target_dim[0]-element[0] for element in coll_xy_dim] ) )

        img_aux=tifff.imread(tiff_stack_path,series=0,key=ch_index,level=nearest_lvl_index)
        img_type=img_aux.dtype.name
        output_img=resize(img_aux,output_shape=target_dim,order=1,preserve_range=True).astype(img_type)

    return output_img

def register_references(fixed_1,fixed_2,moving,mpp,outdir,ch_index):
    transform_log={}
    workdir=Path( os.path.dirname(__file__) )
    trf_params_dir= workdir / "maps" /"transforms"
    reg_params_dir= workdir / "maps" /"registrations"
    qc_dir=outdir / "qc"
    qc_dir_imgs=qc_dir/"imgs"
    (qc_dir_imgs).mkdir(exist_ok=True,parents=True)
    tifff.imwrite(qc_dir_imgs/f"mics_01_src-ch_{ch_index[0]:03d}.tif",fixed_1,photometric="minisblack")
    tifff.imwrite(qc_dir_imgs/f"mics_02_src-ch_{ch_index[1]:03d}.tif",fixed_2,photometric="minisblack")
    tifff.imwrite(qc_dir_imgs/f"maldi_01_src-ch{ch_index[2]:03d}.tif",moving,photometric="minisblack")
    #Step(0) Initial transform: align principal axis
    mofix_1=get_moment_features(fixed_1)
    momov=get_moment_features(moving)
    d_theta=np.rad2deg(mofix_1["major_ax_or"])-np.rad2deg(momov["major_ax_or"])
    d_theta=np.rint(d_theta).astype("int")

    #Write this initial angle in the transforms folder
    transform_log["00_rotation"]=d_theta

    init_transform_dir=trf_params_dir / "00_rotation"
    init_transform_dir.mkdir(exist_ok=True,parents=True)

    #with open(init_transform_dir /'TransformParameters.0.txt', 'w', encoding='utf-8') as f:
    #   f.write(f"angle ({d_theta}) [deg]")

    moving_rot=transform.rotate(img_as_float32(moving),angle=d_theta,resize=True,preserve_range=True,order=0)
    #tifff.imwrite(init_transform_dir/ "result.tif",moving_rot,photometric="minisblack")
    #moving_rot=rescale_intensity(moving_rot,out_range=(np.min(moving),np.max(moving))).astype(moving.dtype.name)
    

    fixed_1itk= itk.GetImageFromArray(img_as_float32(fixed_1))
    fixed_2itk= itk.GetImageFromArray(img_as_float32(fixed_2))
    moving_itk = itk.GetImageFromArray(img_as_float32(moving_rot))
    fixed_1itk.SetSpacing([mpp,mpp])
    fixed_2itk.SetSpacing([mpp,mpp])
    moving_itk.SetSpacing([mpp,mpp])
    #Step(1&2) rigid and affine transform
    
    
    reg_params_files= ["01_rigid.txt","02_affine.txt"]

    reg_params_files=[reg_params_dir/file for file in reg_params_files]
    trf_params_outdir=[trf_params_dir/file.stem for file in reg_params_files]
    for directory in trf_params_outdir:
        directory.mkdir(exist_ok=True,parents=True)

    fixed_itk=[fixed_1itk,fixed_2itk]
    #allocate variable for transformed moving img
    
    for fix,reg_file,trf_outdir in zip(fixed_itk,reg_params_files,trf_params_outdir):
        params_obj=itk.ParameterObject.New()
        params_obj.AddParameterFile(str(reg_file))
        
        moving_itk,result_trf_params = itk.elastix_registration_method(
        fix, 
        moving_itk,
        parameter_object=params_obj,
        output_directory=str(trf_outdir),
        log_file_name="log.txt",
        log_to_console=False
        )
        #save transforms for later implementation on the entire maldi stack
        transform_log[reg_file.stem]=result_trf_params

    #transform_log={00_rotation,01_rigid,02_affine}

    return transform_log
#TODO: re-write this whole function as an python generator
def apply_transform(img,transform_log,mpp_mics,mpp_maldi,upscaled_shape=None):
    # Save original intensity values
    minimum=np.min(img)
    maximum=np.max(img)
    ref_dtype=img.dtype.name
    # Apply initial rotation
    img_upt=transform.rotate(img_as_float32(img),
                            angle=transform_log["00_rotation"],#given in degrees
                            resize=True,
                            preserve_range=True,
                            order=0)
    # Apply rigid transform
    img_upt = itk.transformix_filter(img_upt, transform_log["01_rigid"])
    # Upscale in preparation for affine transform
    # TODO: use upscaled_shape as boolean to decide if upscaled is applied or not
    img_upscaled=transform.resize(img_upt,output_shape=upscaled_shape,order=0,preserve_range=True)
    # Adapt spacing,origin and size of img_upscaled to the source dimensions/space
    src_mpp=np.array([mpp_maldi,mpp_maldi])
    src_origin=np.array([0,0])
    target_mpp=np.array([mpp_mics,mpp_mics])

    target_origin=(0.5*target_mpp)+src_origin-(0.5*src_mpp)
    itk_upscaled=itk.GetImageFromArray(img_upscaled)
    itk_upscaled.SetSpacing(target_mpp.tolist())
    itk_upscaled.SetOrigin(target_origin.tolist())
    affine_params_file=Path().cwd() / "maps" /"transforms"/"02_affine"/"TransformParameters.0.txt"
    #affine_params_object=itk.ParameterObject.New()
    #affine_params_object.ReadParameterFile(str(affine_params_file))
    affine_params_object=transform_log["02_affine"]
    affine_params_object.SetParameter(0,"Size", [str(upscaled_shape[1]),str(upscaled_shape[0])])
    affine_params_object.SetParameter(0,"Spacing", [str(value) for value in target_mpp])
    affine_params_object.SetParameter(0,"Origin", [str(value) for value in target_origin])

    result = itk.transformix_filter(itk_upscaled, affine_params_object)
    result =rescale_intensity(itk.GetArrayFromImage(result),out_range=(minimum,maximum)).astype(ref_dtype)

    return result

def create_pyramid(mics_path,
                    maldi_path,
                    outdir,
                    file_name,
                    transform_log,
                    mpp_mics,
                    mpp_maldi,
                    upscaled_shape):

    outdir.mkdir(parents=True, exist_ok=True)
    out_file_path= outdir / f'{file_name}.tif'
    with tifff.TiffFile(mics_path) as tif:
        sublayers=len(tif.series[0].levels)-1
    sublayers_idx=range(1,sublayers+1)
    

    with tifff.TiffWriter(out_file_path, ome=False, bigtiff=True) as tif:
        #write first the original resolution image,i.e. first layer
        for n,img in enumerate([mics_path,maldi_path]):
            with tifff.TiffFile(img) as stack:
                no_of_channels=len(stack.pages)
                no_of_levels=len(stack.series[0].levels)

            for ch in range(0,no_of_channels):
            #Read first layer and write it 
                if n==0:
                    first_layer=tifff.imread(img,series=0,key=ch,level=0)
                elif n==1:
                    first_layer=apply_transform(
                        tifff.imread(img,series=0,key=ch,level=0),
                        transform_log,
                        mpp_mics,
                        mpp_maldi,
                        upscaled_shape
                        )
                ref_dtype=first_layer.dtype.name
           
                tif.write(
                first_layer,
                description="",
                subifds=sublayers,
                metadata=False,  # do not write tifffile metadata
                tile=(256, 256),
                photometric='minisblack'
                    )
                #Create pyramid if needed
                if no_of_levels>1:
                    pass
                else:
                    pyramid=pyramid_gaussian( first_layer, max_layer=sublayers, preserve_range=True,order=0,sigma=0)
                    next(pyramid)

                # Write sub-layers in file
                for index in sublayers_idx:
                    if n==0:
                        sublayer_img=tifff.imread(img,series=0,key=ch,level=index)
                    elif n==1:
                        sublayer_img=next(pyramid).astype(ref_dtype)
                    tif.write(
                    sublayer_img,
                    subfiletype=1,
                    metadata=False,
                    tile=(256, 256),
                    photometric='minisblack',
                    compression="lzw"#lzw works better when saving channel-by-channel and jpeg 2000 when saving the whole stack at once
                        )
                
    #tifff.tiffcomment(out_file_path, ome_xml)




def main():

    args = CLI.get_args()
    mics_stack_path=args.input_mics
    maldi_stack_path=args.input_maldi
    output_dir=args.output_dir
    mpp_mics=args.pixel_size_microns_mics
    mpp_maldi=args.pixel_size_microns_maldi
    chan_mics=args.reference_indices_mics
    chan_maldi=args.reference_indices_maldi

    #Extract MICS (fixed) and MALDI(moving) channels to be used for registration
    mics=[]
    maldi=tifff.imread(maldi_stack_path,series=0,key=chan_maldi)
    with tifff.TiffFile(mics_stack_path) as tif:
        mics_dimensions=tif.series[0].pages[0].shape
    
    for index in chan_mics:
        mics_resized=extract_and_resize(mics_stack_path,
                                    index,
                                    mpp_mics,
                                    mpp_maldi
                                    )
        mics.append(mics_resized)
    #Extract transforms by using only the reference channels in MICS and MALDI
    transf_log=register_references(mics[0],
                        mics[1],
                        maldi,
                        mpp_maldi,
                        output_dir,
                        chan_mics+[chan_maldi]
                        )
    
    create_pyramid(mics_stack_path,
                    maldi_stack_path,
                    output_dir,
                    "mics_maldi",
                    transf_log,
                    mpp_mics,
                    mpp_maldi,
                    mics_dimensions)
    
    """
    with tifff.TiffFile(maldi_stack_path) as tif:
        for n,page in enumerate(tif.pages):
            result=apply_transform(page.asarray(),transform_log,mpp_mics,mpp_maldi,upscaled_shape=mics_dimensions)
            tifff.imwrite(output_dir /f"{n}.tif",result,photometric="minisblack")
    """

if __name__ == '__main__':
    main()









