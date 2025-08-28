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
import itk
import os
from types import GeneratorType
import tracemalloc
import time
#local libraries
import CLI
import ome_writer


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

def update_markers_info(mics_info,maldi_info):
    df_mics=pd.read_csv(mics_info)
    df_maldi=pd.read_csv(maldi_info,comment="#",sep=";")
    maldi_ch_names=[str(value) for value in df_maldi["m/z"].tolist() ]
    
    column_names=df_mics.columns.tolist()
    aux= pd.DataFrame(index=range(len(maldi_ch_names)),columns=column_names)
    aux["marker_name"]=maldi_ch_names
    aux["channel_number"]=[1+df_mics["channel_number"].max()]*len(maldi_ch_names)
    aux["cycle_number"]=[1+df_mics["cycle_number"].max()]*len(maldi_ch_names)
    aux.fillna("",inplace=True)
    markers_updated=pd.concat([df_mics,aux],ignore_index=True)

    """
    with open(maldi_info, 'r') as file:
        content=file.readlines()
        total_lines=len(content)
        for n,line in enumerate(content):
            if "m/z" in line:
                init_line=n+1
                break
    print("init_line:",init_line)
    maldi_ch_names=[ element.split(";")[0] for element in content[init_line:total_lines] ]
    """
    return markers_updated
        


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
    if is_pyramid(tiff_stack_path)[0]:
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


def register_references(fixed_1,fixed_2,moving,mpp,qc_dir):
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

def apply_transform(img_path,transform_map,mpp_mics,mpp_maldi,upscaled_size=None):
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

def extract_levels_from_tiff(path,ch,levs):
    with tifff.TiffFile(path) as tif:
        for l in range(levs):
            yield tif.series[0].levels[l].pages[ch].asarray()

def create_pyramid(img_instances,
                    levels,
                    outdir,
                    file_name,
                    img_data_type
                    ):

    outdir.mkdir(parents=True, exist_ok=True)
    out_file_path= outdir / f'{file_name}.ome.tif'

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
                if (deficit==0 or deficit<0):

                    unfolded_instances.append(extract_levels_from_tiff(aux_path,ch_idx,levels))

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
                        photometric='minisblack',
                        compression="lzw"
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
                
    return out_file_path


def main(version):

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

    mics_props=extract_img_props(mics_stack_path,mpp_mics)

    mics_dimensions=(mics_props["size_y"],mics_props["size_x"])

    #mics_names=pd.read_csv(file_info_mics)["marker_name"].tolist()

    #maldi_names=extract_maldi_names(file_info_maldi)
    markers_updated=update_markers_info(file_info_mics,file_info_maldi)
    channel_names=markers_updated["marker_name"].tolist()

    #Extract MICS (fixed image(s)) and MALDI(moving image) channels to be used for registration
    maldi=tifff.imread(maldi_stack_path,series=0,key=chan_maldi)
    mics=[]

    #Resize of MICS image to match pixel size of MALDI (downscale expected)
    for index in chan_mics:
        mics_resized=extract_and_resize(mics_stack_path,index,mpp_mics,mpp_maldi)
        mics.append(mics_resized)

    #Extract transforms (transformation map) by using only the reference channels in MICS and MALDI
    transformations_map=register_references(
                        mics[0],
                        mics[1],
                        maldi,
                        mpp_maldi,
                        qc_dir
                        )
    #Apply transforms to maldi image and upscale to the dimensions of the mics data set
    maldi_generator=apply_transform(
                                    maldi_stack_path,
                                    transformations_map,
                                    mpp_mics,
                                    mpp_maldi,
                                    mics_dimensions
                                    )

    #Create pyramidal image of the MICS and MALDI channels
    out_img_path=create_pyramid(
                    [mics_stack_path,maldi_generator],
                    mics_props["levels"],
                    output_dir,
                    "integrated",
                    mics_props["data_type"]
                    )

    #Write metadata in OME format into the pyramidal file
    ome_xml=ome_writer.create_ome(channel_names,mics_props,version)
    tifff.tiffcomment(out_img_path, ome_xml.encode("utf-8"))

    #Write updated markers file containg mics and maldi channels
    markers_updated.to_csv(output_dir/"markers_bs_updated.csv",index=False)


if __name__ == '__main__':
    _version = 'v1.5.0'

    tracemalloc.start()
    st = time.time()

    main(_version)

    print("Memory peak:",((10**(-9))*tracemalloc.get_traced_memory()[1],"GB"))
    rt = time.time() - st
    tracemalloc.stop()
    print(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")









