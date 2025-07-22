import itk
import numpy as np
from skimage import transform
from skimage.util import img_as_float
from skimage.exposure import rescale_intensity
from skimage.util import img_as_float
from tifffile import imwrite
import pathlib








def register_from_params(img_stack,transf_params,mpp,out_dir):
    no_of_channels=img_stack.shape[0]
    out_dir.mkdir(parents=True, exist_ok=True)
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
    


