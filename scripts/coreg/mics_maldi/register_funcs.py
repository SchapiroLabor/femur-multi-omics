import itk
import numpy as np
from skimage import transform
from skimage.util import img_as_float
from skimage.exposure import rescale_intensity
from skimage.util import img_as_float
from tifffile import imwrite

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

def reg_params():
    parameter_object_rigid = itk.ParameterObject.New()
    rigid_parameter_map = parameter_object_rigid.GetDefaultParameterMap('rigid')
    rigid_parameter_map["AutomaticTransformInitialization "]=["true"]
    rigid_parameter_map['ResultImageFormat']= ["tif"]
    rigid_parameter_map['AutomaticTransformInitializationMethod'] = ['CenterOfGravity']
    rigid_parameter_map['NumberOfSpatialSamples']= ["2000"]
    rigid_parameter_map['MaximumNumberOfIterations']=["2000"]
    rigid_parameter_map["Interpolator"]=["NearestNeighborInterpolator"]
    rigid_parameter_map["FinalBSplineInterpolationOrder"]=["0"]
    rigid_parameter_map["ResultImagePixelType"]=["unsigned short"]
    rigid_parameter_map["FixedInternalImagePixelType"]=["unsigned short"]
    rigid_parameter_map["MovingInternalImagePixelType"]=["unsigned short"]
    rigid_parameter_map["FixedInternalImagePixelType"]= ["float"]
    rigid_parameter_map["MovingInternalImagePixelType"]=["float"]
    parameter_object_rigid.AddParameterMap(rigid_parameter_map)

    return parameter_object_rigid

def register_arrays(fix,mov,mpp):
    transform_log={}

    mofix=get_moment_features(fix)
    momov=get_moment_features(mov)
    d_theta=np.rad2deg(mofix["major_ax_or"])-np.rad2deg(momov["major_ax_or"])
    d_theta=np.rint(d_theta).astype("int")

    transform_log["rotation"]=d_theta

    mov_rot=transform.rotate(img_as_float(mov),angle=d_theta,resize=True,preserve_range=True,order=0)
    mov_rot=rescale_intensity(mov_rot,out_range=(np.min(mov),np.max(mov))).astype(mov.dtype.name)

    fix_itk = itk.GetImageFromArray(fix)
    mov_itk = itk.GetImageFromArray(mov_rot)
    fix_itk.SetSpacing([mpp,mpp])
    mov_itk.SetSpacing([mpp,mpp])

    registered_mov_img, result_transform_parameters = itk.elastix_registration_method(
        fix_itk, 
        mov_itk,
        parameter_object=reg_params(),
        log_to_console=False
        )
    
    transform_log["rigid"]=result_transform_parameters


    return registered_mov_img,transform_log

def register_from_params(img_stack,transf_params,mpp,out_dir):
    no_of_channels=img_stack.shape[0]
    for ch in range(0,no_of_channels):
        channel=img_stack[ch,:,:]
        rot=transform.rotate(img_as_float(channel),
                             angle=transf_params["rotation"],#given in degrees
                             resize=True,
                             preserve_range=True,
                             order=0)
        
        rot=rescale_intensity(rot,out_range=(np.min(channel),np.max(channel))).astype(channel.dtype.name)
        rot_itk=itk.GetImageFromArray(rot)
        rot_itk.SetSpacing([mpp,mpp])
        registered=itk.transformix_filter(rot_itk,transf_params["rigid"])

        """
        registered_arr=itk.GetArrayFromImage(registered)
        if ch==0:
            out_stack=np.zeros((no_of_channels,)+registered_arr.shape,dtype=registered_arr.dtype)
        else:
            pass
        out_stack[ch,:,:]=registered_arr
        """
        itk.imwrite( registered,out_dir/ "ch-{x}_maldi_reg.tif".format(x=f"{ch:03d}"))
    #imwrite(out_dir/"maldi_stack_reg.tif",out_stack,photometric="minisblack")
    


