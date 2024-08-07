#Functions for parsing command line arguments for ome ilastik prep
import argparse
from . import __version__ # This still has to be adjusted with __init__.py function


def ParseInputDataExtract():
   """Function for parsing command line arguments for input to single-cell
   data extraction"""

#if __name__ == '__main__':
   parser = argparse.ArgumentParser()
   parser.add_argument('--masks',nargs='+', required=True)
   parser.add_argument('--image', required=True)
   parser.add_argument('--channel_names', required=True)
   parser.add_argument('--output', required=True)
   parser.add_argument(
      '--mask_props', nargs = "+",
      help="""
         Space separated list of additional metrics to be calculated for every mask.
         This is for metrics that depend only on the cell mask. If the metric depends
         on signal intensity, use --intensity-props instead.
         See list at https://scikit-image.org/docs/dev/api/skimage.measure.html#regionprops
      """
   )
   parser.add_argument(
      '--intensity_props', nargs = "+",
      help="""
         Space separated list of additional metrics to be calculated for every marker separately.
         By default only mean intensity is calculated.
         If the metric doesn't depend on signal intensity, use --mask-props instead.
         See list at https://scikit-image.org/docs/dev/api/skimage.measure.html#regionprops

         Additionally available is: gini_index, intensity_median, intensity_sum, intensity_std,
            contrast, dissimilarity, homogeneity, energy, correlation, ASM.
         Further information on the parameters:
         gini_index:
               which calculates a single number between 0 and 1, representing how unequal the signal is distributed in each region.
               Will be calculated for every marker
               See https://en.wikipedia.org/wiki/Gini_coefficient
         contrast, dissimilarity, homogeneity, energy, correlation, ASM:
               glcm/graycoprops features from scikit-image features. The default distance is set to 1 pixel and the default angle is set to 0 rad.
               Both parameters can be controlled via CLI inputs. However, both parameters are limited to 1 input each.
               Will be calculated for every marker
               See https://scikit-image.org/docs/stable/api/skimage.feature.html  
      """
   )
   parser.add_argument(
        '--glcm_angle', type=float, default=0,
        help="Angle in radians for GLCM calculation. Default is 0 radians. Currently limited to 1 angle"
    )
   parser.add_argument(
        '--glcm_distance', type=int, default=1,
        help="Distance in pixels for GLCM calculation. Default is 1 pixel."
    )
   #parser.add_argument('--suffix')
   parser.add_argument('--version', action='version', version=f'mcquant {__version__}')
   args = parser.parse_args()
   #Create a dictionary object to pass to the next function
   dict = {'masks': args.masks, 'image': args.image,\
    'channel_names': args.channel_names,'output':args.output,
    'intensity_props': set(args.intensity_props if args.intensity_props is not None else []).union(["intensity_mean"]),
    'mask_props': args.mask_props,
    'glcm_angle': args.glcm_angle,
    'glcm_distance': args.glcm_distance,
    
   }
   #Print the dictionary object
   print(dict)
   #Return the dictionary
   return dict
