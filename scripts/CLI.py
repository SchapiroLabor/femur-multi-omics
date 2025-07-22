import argparse
import pathlib

#---CLI-BLOCK---#
#TODO: due to the large number of arguments I should try to parse them from a json file(jsonargparse)?
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
                        help='output directory where the registered maldi&mics image will be saved'
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
                    nargs='+',
                    type=int,
                    help="""Two integers separated by a space. 
                    The integers are the 0-based indices of the channels to 
                    be used in the mics stack for an initial coarse registration (usually DAPI) and a 
                    fine one respectively.  For the fine registration the following recommendations are given:
                    1) a marker that is sparse and 2) whose expression is highly correlated with the reference 
                    metabolite from the maldi image (see -idxma argument).
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
                    type=int,
                    help="0-based index of the reference channel in the maldi stack "
                    )
    

    
    args=parser.parse_args()

    return args
#---END_CLI-BLOCK---#