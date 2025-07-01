import argparse
import pathlib

#---CLI-BLOCK---#
def get_args():
    """
    This function parses the command line arguments and returns them as a namespace object.

    returns: namespace object with the arguments.
    """
    parser=argparse.ArgumentParser()

    #Mandatory arguments

    parser.add_argument('-o',
                        '--output_dir',
                        required=True,
                        type=pathlib.Path,
                        help='output directory where the registered maldi images will be saved'
    )

    parser.add_argument('-inmi',
                    '--input_mics',
                    required=True,
                    type=pathlib.Path,
                    help='Path to the macsima image stack (.tif)'
                    )
    
    parser.add_argument('-mppmi',
                    '--pixel_size_microns_mics',
                    required=True,
                    type=float,
                    help='microns per pixel (mpp) of the macsima image'
                    )
    
    parser.add_argument('-idxmi',
                    '--reference_indices_mics',
                    required=True,
                    action="append",
                    type=int,
                    help="""Two integers separated by a space. 
                    The integers are the 0-based indices of the channels to 
                    be used in the stack for an initial coarse registration and a 
                    fine one respectively.
                    """
                    )


    parser.add_argument('-inma',
                    '--input_maldi',
                    required=True,
                    type=pathlib.Path,
                    help='Path to the maldi image stack (.tif)'
                    )
    
    
    parser.add_argument('-mppma',
                    '--pixel_size_microns_maldi',
                    required=True,
                    type=float,
                    help='microns per pixel (mpp) of the maldi image'
                    )
    
    parser.add_argument('-idxma',
                    '--reference_indices_maldi',
                    required=True,
                    action="append",
                    type=int,
                    help="""
                    Two integers separated by a space. 
                    The integers are the 0-based indices of the channels 
                    to be used in the stack for a coarse registration and a fine one respectively.
                    The channels selected here will be registered against the ones selected in the -idxmi
                    argument.
                    """
                    )
    

    
    args=parser.parse_args()

    return args
#---END_CLI-BLOCK---#