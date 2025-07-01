import argparse
import pathlib

def get_args():
    parser=argparse.ArgumentParser()
    #Required arguments
    
    parser.add_argument('-i1',
                    '--input_1',
                    required=True,
                    type=pathlib.Path,
                    help='Path to pyramidal tiff 1'
                    )
    
    parser.add_argument('-i2',
                    '--input_2',
                    required=True,
                    type=pathlib.Path,
                    help='Path to pyramidal tiff 1'
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

    

    args=parser.parse_args()
    return args