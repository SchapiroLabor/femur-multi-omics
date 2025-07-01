import tifffile as tifff
from skimage.transform import pyramid_gaussian
from pathlib import Path
import numpy as np
import ome_writer
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
def write_pyramidal_tiff(imgs,sublayers,outdir,file_name,ome=False,ome_xml="",mpp=1):
    ref_img=tifff.imread(imgs[0])
    width=ref_img.shape[1]
    height=ref_img.shape[0]
    if ome:
        if not ome_xml:
            #create minimal ome-metadata
            marker_filter_map=[ (file.stem.split("_")[2],file.stem.split("_")[3]) for file in imgs]
            info={"size_x":width,
                  "size_y":height,
                  "pixel_size":float(mpp),
                  }
                  
            ome_xml_str=ome_writer.create_ome(marker_filter_map,info)
        else:
            ome_xml_str=ome_xml
    else:
        ome_xml_str=ome_xml


    file_format=np.where([ome],"ome.tiff","tiff")[0]
    outdir.mkdir(parents=True, exist_ok=True)
    out_file_path= outdir / f'{file_name}.{file_format}'


    with tifff.TiffWriter(out_file_path, ome=False, bigtiff=True) as tif:
        #write first the original resolution image,i.e. first layer
        for img in imgs:
            #Read first layer and write it 
            first_layer=tifff.imread(img)
            ref_dtype=first_layer.dtype.name

            tif.write(
            first_layer,
            description="",
            subifds=sublayers,
            metadata=False,  # do not write tifffile metadata
            tile=(256, 256),
            photometric='minisblack'
                    )
            # Create pyramid
            pyramid=pyramid_gaussian( first_layer, max_layer=sublayers, preserve_range=True,order=1,sigma=1)
            # Skip first/original layer of the pyramid
            next(pyramid)
            # Write sub-layers in file
            for sub_level in pyramid:
                tif.write(
                sub_level.astype(ref_dtype),
                subfiletype=1,
                metadata=False,
                tile=(256, 256),
                photometric='minisblack',
                compression="lzw"#lzw works better when saving channel-by-channel and jpeg 2000 when saving the whole stack at once
                        )
                
    tifff.tiffcomment(out_file_path, ome_xml)
                





def main():
    args = CLI.get_args()
    
    #write_pyramidal_tiff(imgs,sublayers,outdir,file_name,ome=False,ome_xml="")
    write_pyramidal_tiff(
                        list( sorted( (args.input).glob("*.tif") ) ),
                        args.pyramid_levels, 
                        args.output,
                        args.stack_name,
                        args.ome_format,
                        args.ome_xml,
                        args.microns_per_pixel
                        )
    



if __name__ == '__main__':
    main()



