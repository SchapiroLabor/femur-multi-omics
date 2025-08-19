#public libraries
from pathlib import Path
from skimage import transform
import itertools
import pandas as pd
from skimage.transform import pyramid_gaussian
import tifffile as tifff
import numpy as np
from skimage.transform import resize
from skimage.util import img_as_float32
from skimage.exposure import rescale_intensity
import itk
import os
from types import GeneratorType
#local libraries
import CLI


#HELPERS
def is_pyramid(img_path):
    """
    Checks if image has pyramidal levels

    """
    
    with tifff.TiffFile(img_path) as tif:
        levels=len(tif.series[0].levels)
        pyramid=levels > 1
    return pyramid,levels


def extract_img_props(img_path,pixel_size):

    #Checks if image has pyramidal levels
    with tifff.TiffFile(img_path) as tif:
        pyr_levels=len(tif.series[0].levels)
        is_pyramid=pyr_levels > 1
        data_type=tif.series[0].dtype.name
        height,width=tif.series[0].shape[-2::]
    
    img_props={"pixel_size":pixel_size,
               "pixel_size_unit":"µm",
               "data_type":data_type,
               "pyramid":is_pyramid,
               "levels":pyr_levels,
               "size_x":width,
               "size_y":height ,
               "bits": "16"
               }
    
    return img_props

def extract_maldi_names(maldi_info):
    
    with open(maldi_info, 'r') as file:
        content=file.readlines()
        total_lines=len(content)
        for n,line in enumerate(content):
            if "m/z" in line:
                init_line=n+1
                break
    print("init_line:",init_line)
    maldi_ch_names=[ element.split(";")[0] for element in content[init_line:total_lines] ]

    return maldi_ch_names
        



def get_moment_features(img_arr,mpp):
    features={}
    img=itk.GetImageFromArray(img_arr)
    img.SetSpacing([mpp,mpp])
    width_units=mpp*np.array(itk.size(img))[0]
    height_units=mpp*np.array(itk.size(img))[1]
    width_pix=np.array(itk.size(img))[0]
    height_pix=np.array(itk.size(img))[1]
    moments = itk.ImageMomentsCalculator.New(img)
    moments.Compute()
    principal_axis=itk.array_from_matrix(moments.GetPrincipalAxes())
    minor_axis=principal_axis[0]
    major_axis=principal_axis[1]

    x_min=minor_axis[0]
    y_min=minor_axis[1]
    x_max=major_axis[0]
    y_max=major_axis[1]

    
    centroid=tuple(moments.GetCenterOfGravity())
    #centroid x,y, with origin on top left corner (y inverted)
    features["centroid_XinvY"]=(centroid[0],centroid[1])
    features["centroid_pixels"]=(centroid[0]//mpp,centroid[1]//mpp)
    #centroid x,y with origin on bottom left corner (standard axis)
    features["centroid_XY"]=(centroid[0],height_units-centroid[1])
    features["major_ax"]=major_axis
    features["minor_ax"]=minor_axis
    features["major_ax_or"]=np.arctan(x_max/y_max)#[radians]relative to y+ axis clockwise (/) negative, counterclowise (\) positive
    features["minor_ax_or"]=np.arctan(x_min/y_min)
    features["size"]=[width_pix,height_pix]
    return features

def get_crop_params(width,height,Transform):

    corners = np.array([[0, 0], 
                    [0, height- 1], 
                    [width - 1, height - 1], 
                    [width - 1, 0]
                    ]
                )
    
    corners_ = np.rint(Transform.inverse(corners))
    min_x = corners_[:, 0].min()
    min_y = corners_[:, 1].min()
    max_x = corners_[:, 0].max()
    max_y = corners_[:, 1].max()
    out_rows = max_y - min_y + 1
    out_cols = max_x - min_x + 1
    output_shape = np.rint([out_cols,out_rows]).astype("int")
    translation_offset=[int(min_x), int(min_y)]

    return output_shape,translation_offset 



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
        output_img=resize(img_aux,output_shape=target_dim,order=0,preserve_range=True).astype(img_type)

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

def apply_init_transform(fix,moving,mpp,template_file):
    mofix=get_moment_features(fix,mpp)
    momov=get_moment_features(moving,mpp)

    theta=mofix["major_ax_or"]-momov["major_ax_or"]                                 
    xc,inv_yc=momov["centroid_XinvY"]#centroid with inverted-y
    xc_pixels,inv_yc_pixels=momov["centroid_pixels"]
    #Calculate dimensions of output image after initial transformation
    width,height=momov["size"]
    tr1=transform.EuclideanTransform(rotation=theta,translation=(xc_pixels,inv_yc_pixels))
    tr2=transform.EuclideanTransform(translation=(-xc_pixels,-inv_yc_pixels))
    Transform=tr2+tr1 #T is the matrix producing the rotation around centroid
    outsize,offset=get_crop_params(width,height,Transform)
    x_offset= mpp*offset[0]
    y_offset= -mpp*offset[1]
    #Use Transformix parameter map to apply the initial transform
    transform_map=itk.ParameterObject.New()
    transform_map.ReadParameterFile(str(template_file))

    transform_map.SetParameter(0,"Transform",["EulerTransform"])
    transform_map.SetParameter(0,"NumberOfParameters",["3"])
    #CenterofRotation arguments are x,y in units with y value in inverted y-axis
    transform_map.SetParameter(0,"CenterOfRotationPoint",[ str(xc) , str(inv_yc) ])
    transform_map.SetParameter(0,"TransformParameters", [ str(theta),str( y_offset), str(x_offset)])
    transform_map.SetParameter(0,"Size", [str(val) for val in outsize ])
    transform_map.SetParameter(0,"Spacing",[str(mpp),str(mpp)])
    transform_map.SetParameter(0,"Origin",["0","0"])

    mov_itk=itk.GetImageFromArray( img_as_float32(moving) )
    mov_itk.SetSpacing([mpp,mpp])
    transformed_image=itk.transformix_filter(mov_itk,transform_map,log_to_console=False)


    return transformed_image,transform_map


def register_references2(fixed_1,fixed_2,moving,mpp,qc_dir):
    #Define variables
    transform_scheme=["00_init","01_rigid","02_affine"]
    global_trf_map=itk.ParameterObject.New()
    workdir=Path( os.path.dirname(__file__) )
    transform_template=workdir /"maps"/"transforms" /"template"/"TransformParameters.0.txt"
    registration_maps=sorted( (workdir / "maps" / "registrations").glob("*.txt") )
    out_transform_dir=qc_dir / "transforms"

    qc_out=[]

    for element in transform_scheme:
        aux=out_transform_dir / element
        aux.mkdir(exist_ok=True,parents=True)
        qc_out.append(aux)
        
    fix_itk=[]
    for img in [fixed_1,fixed_2]:
        aux=itk.GetImageFromArray(img_as_float32(img))
        aux.SetSpacing([mpp,mpp])
        fix_itk.append(aux)
    #Initial transformation to align the principal axis
    #mov_updated below is already an itk image with the mpp spacing
    mov_updated,transform_map=apply_init_transform(fixed_1,moving,mpp,transform_template)

    transform_map.WriteParameterFile(transform_map.GetParameterMap(0),
                                     str(qc_out[0] / "TransformParameters.0.txt")
                                     )
    
    itk.imwrite(mov_updated,qc_out[0] / "init.tif")
    
    global_trf_map.AddParameterMap(transform_map.GetParameterMap(0))

    for Fix,Reg,Out in zip(fix_itk,registration_maps,qc_out[1::]):
        reg_map=itk.ParameterObject.New()
        reg_map.AddParameterFile(str(Reg))

        mov_updated,result_trf_params = itk.elastix_registration_method(
            Fix, 
            mov_updated,
            parameter_object=reg_map,
            output_directory=str(Out),
            log_file_name="log.txt",
            log_to_console=False
        )

        global_trf_map.AddParameterMap(result_trf_params.GetParameterMap(0))

    return global_trf_map

def apply_transform2(img_path,transform_map,mpp_mics,mpp_maldi,upscaled_size=None):
    #index of last transformation map

    no_of_transforms=transform_map.GetNumberOfParameterMaps()

    if upscaled_size==None:
        pass
    else:

        last_index=no_of_transforms-1
        src_mpp=np.array([mpp_maldi,mpp_maldi])
        src_origin=np.array([0,0])
        target_mpp=np.array([mpp_mics,mpp_mics])
        target_origin=(0.5*target_mpp)+src_origin-(0.5*src_mpp)
        transform_map.SetParameter(last_index,"Size", [str(upscaled_size[1]),str(upscaled_size[0])])
        transform_map.SetParameter(last_index,"Spacing", [str(value) for value in target_mpp])
        transform_map.SetParameter(last_index,"Origin", [str(value) for value in target_origin])
    
    itk_imgs=[]

    with tifff.TiffFile(img_path) as tif:
        for channel in tif.pages:
            aux=itk.GetImageFromArray( channel.asarray() )
            aux.SetSpacing([mpp_maldi,mpp_maldi])
            itk_imgs.append(aux)

    for ch in itk_imgs:
        result=ch
        for n in range(no_of_transforms):
            single_transform=itk.ParameterObject.New()
            single_transform.AddParameterMap(transform_map.GetParameterMap(n))
            result=itk.transformix_filter(result,single_transform,log_to_console=False)
        yield itk.GetArrayFromImage( result )

    

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

def create_pyramid2(img_instances,
                    levels,
                    outdir,
                    file_name,
                    img_data_type
                    ):

    outdir.mkdir(parents=True, exist_ok=True)
    out_file_path= outdir / f'{file_name}.tif'

    #
    types=[]
    for element in img_instances:

        if isinstance(element, GeneratorType):

            types.append("generator")

        elif isinstance(element, Path):

            types.append("path")

        else:

            types.append("other")
    
    #
    pyramid_levels=[]
    for path,element in zip(img_instances,types):

        if element=="path":

            pyramid_levels.append( is_pyramid(path)[1] )

        else:

            pyramid_levels.append(1)

    #
    unfolded_instances=[]  
    for INST,TYPE,LEVL in zip(img_instances,types,pyramid_levels):
        
        if TYPE=="path":
            deficit=levels-LEVL
            aux_path=INST
            print(deficit)
            with tifff.TiffFile(aux_path) as tif:
                no_channels=len(tif.pages)

            for ch_idx in range(no_channels):
                if deficit==0:
                    unfolded_instances.append( (tifff.imread(aux_path,series=0,key=ch_idx,level=L) for L in range( abs(LEVL) ) ) )

                elif deficit<0:
                    unfolded_instances.append( (tifff.imread(aux_path,series=0,key=ch_idx,level=L) for L in range(levels) ) )

                elif deficit>0:
                    aux_1=(tifff.imread(aux_path,series=0,key=ch_idx,level=L) for L in range(LEVL) )
                    aux_2=pyramid_gaussian( tifff.imread(aux_path,series=0,key=ch_idx,level=LEVL-1), max_layer=deficit, preserve_range=True,order=1,sigma=1)
                    next(aux_2)
                    unfolded_instances.append(itertools.chain(aux_1,aux_2))

        if TYPE=="generator":
            for channel in INST:
                unfolded_instances.append( pyramid_gaussian( channel  , max_layer=levels-1, preserve_range=True,order=1,sigma=1) )

    
    sublayers=levels-1
    with tifff.TiffWriter(out_file_path, ome=False, bigtiff=True) as tif:
        #write first the original resolution image,i.e. first layer
        for img_generator in unfolded_instances:
            for layer,img_layer in enumerate(img_generator):
                if layer==0:
                    tif.write(
                        img_layer.astype(img_data_type),
                        description="",
                        subifds=sublayers,
                        metadata=False,  # do not write tifffile metadata
                        tile=(256, 256),
                        photometric='minisblack'
                        )
                elif layer>0:
                    tif.write(
                        img_layer.astype(img_data_type),
                        subfiletype=1,
                        metadata=False,
                        tile=(256, 256),
                        photometric='minisblack',
                        compression="lzw"#lzw works better when saving channel-by-channel and jpeg 2000 when saving the whole stack at once
                        )
                
    #tifff.tiffcomment(out_file_path, ome_xml)


def main():

    #Collect arguments
    args = CLI.get_args()
    mics_stack_path=args.input_mics
    maldi_stack_path=args.input_maldi
    mpp_mics=args.pixel_size_microns_mics
    mpp_maldi=args.pixel_size_microns_maldi
    chan_mics=args.reference_indices_mics
    chan_maldi=args.reference_indices_maldi
    file_info_mics=args.marker_names_mics
    file_info_maldi=args.marker_names_maldi
    output_dir=args.output_dir
    #Initialize variables & functions
    qc_dir= output_dir / "qc_coreg"

    qc_dir.mkdir(parents=True,exist_ok=True)

    dst_info=extract_img_props(mics_stack_path,mpp_mics)

    mics_dimensions=(dst_info["size_y"],dst_info["size_x"])

    mics_names=pd.read_csv(file_info_mics)["marker_name"].tolist()

    maldi_names=extract_maldi_names(file_info_maldi)

    channel_names=mics_names + maldi_names
    #Extract MICS (fixed image(s)) and MALDI(moving image) channels to be used for registration
    maldi=tifff.imread(maldi_stack_path,series=0,key=chan_maldi)
    mics=[]
        #Resize of MICS image to match pixel size of MALDI (downscale expected)
    for index in chan_mics:
        mics_resized=extract_and_resize(mics_stack_path,
                                    index,
                                    mpp_mics,
                                    mpp_maldi
                                    )
        print(mics_resized.shape)
        mics.append(mics_resized)
    #Extract transforms (transformation map) by using only the reference channels in MICS and MALDI
    transformations_map=register_references2(
                        mics[0],
                        mics[1],
                        maldi,
                        mpp_maldi,
                        qc_dir
                        )
    
    generator=apply_transform2(maldi_stack_path,transformations_map,mpp_mics,mpp_maldi)#,mics_dimensions)
    out_test=Path(r"C:\Users\VictorP\repos\femur_data\output")
    for n,img in enumerate(generator):
        #itk.imwrite(img,out_test / f"maldi_{n:03d}.tif")
        tifff.imwrite(out_test / f"maldi_{n:03d}.tif",img,photometric="minisblack")

    """

    create_pyramid2([mics_stack_path,generator],
                    9,
                    out_test,
                    "integrated",
                    "uint16"
                    )

    """

    """
    create_pyramid(mics_stack_path,
                    maldi_stack_path,
                    output_dir,
                    "mics_maldi",
                    transf_log,
                    mpp_mics,
                    mpp_maldi,
                    mics_dimensions)
    
    
    with tifff.TiffFile(maldi_stack_path) as tif:
        for n,page in enumerate(tif.pages):
            result=apply_transform(page.asarray(),transform_log,mpp_mics,mpp_maldi,upscaled_shape=mics_dimensions)
            tifff.imwrite(output_dir /f"{n}.tif",result,photometric="minisblack")

    #Write metadata in OME format into the pyramidal file
    channel_names=markers_updated["marker_name"].tolist()
    ome_xml=ome_writer.create_ome(channel_names,src_props,version)
    tifff.tiffcomment(pyramid_abs_path, ome_xml.encode("utf-8"))

    #Write updated markers.csv
    markers_updated = markers_updated.drop(columns=['keep','ind','processed','factor','bg_idx'])
    markers_updated .to_csv(args.markerout, index=False)
    """

if __name__ == '__main__':
    main()









