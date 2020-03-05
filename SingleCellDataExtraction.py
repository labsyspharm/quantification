#Functions for reading in single cell imaging data
#Joshua Hess

#Import necessary modules
import skimage.io
import h5py
import pandas as pd
import numpy as np
import os
import skimage.measure as measure
from pathlib import Path


def MaskChannel(mask,channel):
    """Function for masking a single channel image and returning stats

    Returns a table with CellID according to the mask and the mean pixel intensity
    for the given channel for each cell"""
    #Perform the masking to get all measures
    dat = measure.regionprops(mask,channel)
    #Extract the stats
    intensity = []
    for i in range(len(dat)):
        intensity.append(dat[i].mean_intensity)
    #Return the dataframe object
    return intensity


def MaskIDs(mask):
    """This function will extract the CellIDs and the XY positions for each
    cell based on that cells centroid

    Returns a dictionary object"""
    #Get the CellIDs for this dataset
    dat = measure.regionprops(mask)
    #Extract the CellIDs
    labels = []
    xcoords = []
    ycoords = []
    area = []
    minor_axis_length=[]
    major_axis_length=[]
    eccentricity = []
    solidity = []
    extent=[]
    orientation=[]
    for i in range(len(dat)):
        #Get cell labels
        labels.append(dat[i].label)
        #Get x coordinate
        xcoords.append(dat[i].centroid[0])
        #Get y coordinate
        ycoords.append(dat[i].centroid[1])
        #Get Area
        area.append(dat[i].area)
        #Get the major_axis_length
        major_axis_length.append(dat[i].major_axis_length)
        #Get the minor_axis_length
        minor_axis_length.append(dat[i].minor_axis_length)
        #Get the eccentricity
        eccentricity.append(dat[i].eccentricity)
        #Get the solidity
        solidity.append(dat[i].solidity)
        #Get the extent
        extent.append(dat[i].extent)
        #Get the extent
        orientation.append(dat[i].orientation)

    #Form a dataframe from the lists
    IDs = {"CellID": labels, "X_position": xcoords, "Y_position": ycoords,"Area":area,\
        "MajorAxisLength":major_axis_length,"MinorAxisLength":minor_axis_length,\
        "Eccentricity":eccentricity,"Solidity":solidity,"Extent":extent,"Orientation":orientation}
    #Return the IDs object
    return IDs


def PrepareData(mask,z_stack,channel_names):
    """Function for preparing input for maskzstack function. Connecting function
    to use with mc micro ilastik pipeline"""

    #Get the image path
    im_name = Path(z_stack)
    #Create pathlib object for channel_names
    channel_names = Path(channel_names)

    #Check to see if the image is tiff
    if im_name.suffix == '.tiff' or im_name.suffix == '.tif':
        #Read the image
        z_stack = skimage.io.imread(im_name,plugin='tifffile')
        #Switch the axis order from cyx to yxc - consistent with reading single channel tif (mask)
        z_stack = np.swapaxes(z_stack,0,2)
        z_stack = np.swapaxes(z_stack,0,1)

    #Check to see if image is hdf5
    elif im_name.suffix == '.h5' or im_name.suffix == '.hdf5':
        #Read the image
        f = h5py.File(im_name,'r+')
        #Get the dataset name from the h5 file
        dat_name = list(f.keys())[0]
        ###If the hdf5 is exported from ilastik fiji plugin, the dat_name will be 'data'
        #Get the image data
        z_stack = np.array(f[dat_name])
        #Remove the first axis (ilastik convention)
        z_stack = z_stack.reshape((z_stack.shape[1],z_stack.shape[2],z_stack.shape[3]))
        ###If the hdf5 is exported from ilastik fiji plugin, the order will need to be
        ###switched as above --> z_stack = np.swapaxes(z_stack,0,2) --> z_stack = np.swapaxes(z_stack,0,1)

    #Read the mask
    mask = skimage.io.imread(mask,plugin='tifffile')
    #Read the channels names to pass to the
    channel_names = pd.read_csv(channel_names,header=None)
    #Add a column index for ease
    channel_names.columns = ["marker"]
    #Convert the channel names to a list
    channel_names = list(channel_names.marker.values)

    #Return the objects
    return mask, z_stack, channel_names


def MaskZstack(mask,z_stack,channel_names):
    """This function will extract the stats for each cell mask through each channel
    in the input image

    mask: Tiff image mask that represents the cells in your image. Must end with the word mask!!

    z_stack: Multichannel z stack image"""

    #Get the CellIDs for this dataset
    IDs = pd.DataFrame(MaskIDs(mask))
    #Iterate through the z stack to extract intensities
    list_of_chan = []
    for z in range(z_stack.shape[2]):
        #Get the z channel and the associated channel name from list of channel names
        list_of_chan.append(MaskChannel(mask,z_stack[:,:,z]))
    #Convert the channel names list and the list of intensity values to a dictionary and combine with CellIDs and XY
    dat = pd.concat([IDs,pd.DataFrame(dict(zip(channel_names,list_of_chan)))],axis=1)
    #Get the name of the columns in the dataframe so we can reorder to histoCAT convention
    cols = list(dat.columns.values)
    #Reorder the list (Move xy position to end with spatial information)
    cols.append(cols.pop(cols.index("X_position")))
    cols.append(cols.pop(cols.index("Y_position")))
    cols.append(cols.pop(cols.index("Area")))
    cols.append(cols.pop(cols.index("MajorAxisLength")))
    cols.append(cols.pop(cols.index("MinorAxisLength")))
    cols.append(cols.pop(cols.index("Eccentricity")))
    cols.append(cols.pop(cols.index("Solidity")))
    cols.append(cols.pop(cols.index("Extent")))
    cols.append(cols.pop(cols.index("Orientation")))
    #Reindex the dataframe with new order
    dat = dat.reindex(columns=cols)
    #Return the dataframe
    return dat


def ExtractSingleCells(mask,z_stack,channel_names,output):
    """Function for extracting single cell information from input
    path containing single-cell masks, z_stack path, and channel_names path."""

    #Get the name of the z_stack
    z_stack = Path(z_stack)
    im_name = z_stack.stem
    #Create pathlib object for output
    output = Path(output)
    #Get the name for the mask in the mask directory
    mask_name = Path(mask)
    #Run the data Prep function
    mask, z_stack, channel_names = PrepareData(mask_name,z_stack,channel_names)
    #Use the above information to mask z stack
    scdata = MaskZstack(mask,z_stack,channel_names)
    #Write the singe cell data to a csv file using the image name
    scdata.to_csv(str(Path(os.path.join(str(output),str(im_name+".csv")))),index=False)


def MultiExtractSingleCells(masks,z_stacks,channel_names,output):
    """Function for iterating over a list of z_stacks and output locations to
    export single-cell data from image masks"""

    #Iterate over each image in the list if only a single output
    if len(output) < 2:
        #Check to make sure the masks and image paths are equal in length
        if len(masks) != len(z_stacks):
            raise(ValueError("Number of masks not equal to number of z-stacks"))
        #Alternatively, iterate over masks and z stacks
        else:
            #Iterate through the images and export to the same location
            for i in range(len(z_stacks)):
                #Print update
                print("Extracting sinle-cell data for "+str(z_stacks[i])+'...')
                #Run the ExtractSingleCells function for this image
                ExtractSingleCells(masks[i],z_stacks[i],channel_names,output[0])
                #Print update
                print("Finished "+str(z_stacks[i]))
    #Alternatively, iterate over output directories
    else:
        #Check to make sure the output directories and image paths are equal in length
        if len(output) != len(z_stacks):
            raise(ValueError("Detected more than one output but not as many directories as images"))
        else:
            #Check to make sure the masks and image paths are equal in length
            if len(masks) != len(z_stacks):
                raise(ValueError("Number of masks not equal to number of z-stacks"))
            #Alternatively, iterate over masks and z stacks
            else:
                #Iterate through images and output directories
                for i in range(len(z_stacks)):
                    #Print update
                    print("Extracting sinle-cell data for "+str(z_stacks[i])+'...')
                    #Run the ExtractSingleCells function for this image and output directory
                    ExtractSingleCells(masks[i],z_stacks[i],channel_names,output[i])
                    #Print update
                    print("Finished "+str(z_stacks[i]))











#
