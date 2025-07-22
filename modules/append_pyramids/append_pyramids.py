import tifffile as tifff
from skimage.transform import pyramid_gaussian
from pathlib import Path
import numpy as np
import CLI



"""
def create_pyramid(img,sub_layers):
    if isinstance(workdir,Path):
        img=tifff.imread(img)
    elif isinstance(img,np.ndarray):
        pass
    pyramid=pyramid_gaussian( img, max_layer=sub_layers, preserve_range=True,order=1,sigma=1)
    return pyramid
"""
def append_pyramids(imgs,sublayers,outdir,file_name):

    outdir.mkdir(parents=True, exist_ok=True)
    out_file_path= outdir / f'{file_name}.tif'
    sublayers_idx=range(1,sublayers+1)

    with tifff.TiffWriter(out_file_path, ome=False, bigtiff=True) as tif:
        #write first the original resolution image,i.e. first layer
        for img in imgs:
            with tifff.TiffFile(img) as pyramid:
                no_of_channels=len(pyramid.pages)

            for ch in range(0,no_of_channels):
            #Read first layer and write it 
                first_layer=tifff.imread(img,series=0,key=ch,level=0)
                ref_dtype=first_layer.dtype.name
           
                tif.write(
                first_layer,
                description="",
                subifds=sublayers,
                metadata=False,  # do not write tifffile metadata
                tile=(256, 256),
                photometric='minisblack'
                    )

                # Write sub-layers in file
                for index in sublayers_idx:
                    sublayer_img=tifff.imread(img,series=0,key=ch,level=index)
                    
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
    imgs=[args.input_1,args.input_2]
    append_pyramids(imgs,
                    args.pyramid_levels, 
                    args.output,
                    args.stack_name
                    )
    
if __name__ == '__main__':
    main()



