#Functions for parsing command line arguments for ome ilastik prep
import argparse


def ParseInputDataExtract():
   """Function for parsing command line arguments for input to single-cell
   data extraction"""

#if __name__ == '__main__':
   parser = argparse.ArgumentParser()
   parser.add_argument('--masks')
   parser.add_argument('--z_stacks',nargs='*')
   parser.add_argument('--channel_names')
   parser.add_argument('--output',nargs='*')
   args = parser.parse_args()
   #Create a dictionary object to pass to the next function
   dict = {'masks': args.masks, 'z_stacks': args.z_stacks,\
    'channel_names': args.channel_names,'output':args.output}
   #Print the dictionary object
   print(dict)
   #Return the dictionary
   return dict
