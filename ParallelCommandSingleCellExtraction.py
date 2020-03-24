#Functions for parallel processing single-cell quantification
#Joshua Hess

#Import necessary modules
import skimage.io
import h5py
import pandas as pd
import numpy as np
import os
import skimage.measure as measure
from pathlib import Path
import shutil
import tracemalloc
from joblib import Parallel, delayed
import argparse

#Import custom modules from regular quantification module
import SingleCellDataExtraction as quant


def ParallelMaskZstack(mask_loaded,image,channel_names_loaded,output,n_jobs):
    """This function will extract the stats for each cell mask through each channel
    in the input image
    mask: Tiff image mask that represents the cells in your image. Must end with the word mask!!
    z_stack: Multichannel z stack image"""

    def ParallelChannelMask(z):
        """Function for looping through z_stack and exporting single-cell measures
        Connector for parallel processing function and defined within this
        function so that objects defined inside ParallelMaskZstack can be used"""

        #Create a name for this channel
        exChan = Path(os.path.join(tmp_dir,im_name.stem.replace('.ome',''))+"_"+str(channel_names_loaded[z])+".csv")
        #Run the data Prep function
        image_loaded_z = quant.PrepareData(image,z)
        #Use the above information to mask z stack
        z_chan = pd.DataFrame(quant.MaskChannel(mask_loaded,image_loaded_z),columns=[str(channel_names_loaded[z])])
        #Write out the z_channel to a csv
        z_chan.to_csv(exChan,index=False)
        #Remove the channel from memory
        z_chan = 0
        #Remove the z channel from memory
        image_loaded_z = 0
        #Return the file name
        return exChan

    #Create a list to keep the filenames in
    fnames = []
    #Ensure the im_name is a pathlib object
    im_name = Path(image)
    #Create pathlib object for output
    output = Path(output)
    #Get home directory (Current)
    home_dir = os.getcwd()
    #Create a temporary directory to export files to
    tmp_dir = Path(os.path.join(str(output),im_name.stem.replace('.ome','')+"_tmp"))
    os.mkdir(str(tmp_dir))
    #Change to this directory
    os.chdir(str(tmp_dir))

    #Get the CellIDs for this dataset
    IDs = pd.DataFrame(quant.MaskIDs(mask_loaded))
    #Create a name for the cell ID file
    exIDs = Path(os.path.join(tmp_dir,im_name.stem)+"_cellIDs.csv")
    #Write out the cellIDs to csv
    IDs.to_csv(exIDs,index=False)
    #Update the list of filenames
    fnames.append(exIDs)

    #Remove the IDs from memory
    IDs = 0

    #Iterate through the z stack to extract intensities (Parallelized-joblib maintains marker order)
    scdata_fnames_z = Parallel(n_jobs=n_jobs)(delayed(ParallelChannelMask)(z) for z in range(len(channel_names_loaded)))
    #Update the return list of file names
    fnames = fnames + scdata_fnames_z

    #Change back to the home directory
    os.chdir(home_dir)
    #Return a list containing the file names that are individually exported
    return fnames, tmp_dir


def ParallelExtractSingleCells(mask,image,channel_names,output,n_jobs):
    """Function for extracting single cell information from input
    path containing single-cell masks, z_stack path, and channel_names path."""


    #Read the channels names
    channel_names_loaded = pd.read_csv(channel_names,header=None)
    #Add a column index for ease
    channel_names_loaded.columns = ["marker"]
    #Convert the channel names to a list
    channel_names_loaded = list(channel_names_loaded.marker.values)
    #Read the mask
    mask_loaded = skimage.io.imread(mask,plugin='tifffile')
    #Write single cell information for each channel individually
    scdata_fnames_z, tmp_dir = ParallelMaskZstack(mask_loaded,image,channel_names_loaded,output,n_jobs)

    #Remove the mask from memory
    mask_loaded = 0

    #Concatenate all of the channel files and cell IDs
    scdata_z = quant.ConcatenateData(scdata_fnames_z,tmp_dir)

    #Write the singe cell data to a csv file using the image name
    im_full_name = os.path.basename(image)
    im_name = im_full_name.split('.')[0]
    scdata_z.to_csv(str(Path(os.path.join(str(output),str(im_name+".csv")))),index=False)

    #Remove the single cell data from memory
    scdata_z = 0


def ParallelMultiExtractSingleCells(mask,image,channel_names,output,n_jobs):
    """Function for iterating over a list of z_stacks and output locations to
    export single-cell data from image masks"""

    print("Extracting single-cell data for "+str(image)+'...')
    #Start memory monitoring --- next line can be commented
    tracemalloc.start()
    # ---

    #Run the ExtractSingleCells function for this image
    ParallelExtractSingleCells(mask,image,channel_names,output,int(n_jobs))

    #Display memory monitoring --- next 3 lines can be commented
    current, peak = tracemalloc.get_traced_memory()
    print(f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB")
    tracemalloc.stop()
    # ---

    #Print update
    im_full_name = os.path.basename(image)
    im_name = im_full_name.split('.')[0]
    print("Finished "+str(im_name))


def ParseInputDataExtract():
   """Function for parsing command line arguments for input to single-cell
   data extraction parallel processing"""

   parser = argparse.ArgumentParser()
   parser.add_argument('--mask')
   parser.add_argument('--image')
   parser.add_argument('--channel_names')
   parser.add_argument('--output')
   parser.add_argument('--n_jobs')
   #parser.add_argument('--suffix')
   args = parser.parse_args()
   #Create a dictionary object to pass to the next function
   dict = {'mask': args.mask, 'image': args.image,\
    'channel_names': args.channel_names,'output':args.output,\
    'n_jobs':args.n_jobs}
   #Print the dictionary object
   print(dict)
   #Return the dictionary
   return dict


if __name__ == '__main__':
    #Parse the command line arguments
    args = ParseInputDataExtract()

    #Run the MultiExtractSingleCells function
    ParallelMultiExtractSingleCells(**args)
