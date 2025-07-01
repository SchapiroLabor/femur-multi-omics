import argparse
import pathlib

def get_args():
    parser=argparse.ArgumentParser()
    #Required arguments
    
    parser.add_argument('-i',
                    '--input',
                    required=True,
                    type=pathlib.Path,
                    help='Path to folder with images'
                    )
    
    parser.add_argument('-o',
                    '--output',
                    required=True,
                    type=pathlib.Path,
                    help='Path where the stack will be saved. If directory does not exist it will be created.'
                    )
    
    parser.add_argument('-n',
                    '--stack_name',
                    required=True,
                    help='string specifying the name of the pyramidal image file'
                    )
    #Optional arguments
    parser.add_argument('-pl',
                    '--pyramid_levels',
                    type=int,
                    default=4,
                    help='integer specifying number of sub-levels in the pyramid'
                    )
    
    parser.add_argument('-mpp',
                    '--microns_per_pixel',
                    type=float,
                    default=None,
                    help='integer specifying number of sub-levels in the pyramid'
                    )
    
    parser.add_argument('-ome',
                    '--ome_format',
                    action='store_true',
                    help='flag to store the file as ome.tif'
                    )
    
    parser.add_argument('-omx',
                    '--ome_xml',
                    type=str,
                    default="",
                    help='custom ome-xml string containing the metadata. If nothing is passed a minimal metadata will be created'
                    )
    
    parser.add_argument('-res',
                    '--resize_input',
                    action='store_true',
                    help='flag to resize the input image'
                    )
    
    parser.add_argument('-rdim',
                    '--resize_dimensions',
                    action='extend',
                    nargs=2,
                    type=int,
                    help='two-entry tuple with yx dimensions, i.e. size of height & width '
                    )
    
    args=parser.parse_args()
    return args