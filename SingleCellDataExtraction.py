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
import tracemalloc



def MaskChannel(mask_loaded,image_loaded_z):
    """Function for quantifying a single channel image

    Returns a table with CellID according to the mask and the mean pixel intensity
    for the given channel for each cell"""
    #Perform the masking to get all mean intensities
    dat = measure.regionprops(mask_loaded,image_loaded_z)
    # collect all intensities
    intensity_z = []
    for i in range(len(dat)):
        intensity_z.append(dat[i].mean_intensity)
    # return per channel
    return intensity_z


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


def PrepareData(image,z):
    """Function for preparing input for maskzstack function. Connecting function
    to use with mc micro ilastik pipeline"""

    image_path = Path(image)

    #Check to see if image tif(f)
    if image_path.suffix == '.tiff' or image_path.suffix == '.tif':
        #Check to see if the image is ome.tif(f)
        if  image.endswith(('.ome.tif','.ome.tiff')):
            #Read the image
            image_loaded_z = skimage.io.imread(image,img_num=z,plugin='tifffile')
            #print('OME TIF(F) found') 
        else:
            #Read the image
            image_loaded_z = skimage.io.imread(image,img_num=z,plugin='tifffile')
            #print('TIF(F) found')
            # Remove extra axis
            #image_loaded = image_loaded.reshape((image_loaded.shape[1],image_loaded.shape[3],image_loaded.shape[4]))

    #Check to see if image is hdf5
    elif image_path.suffix == '.h5' or image_path.suffix == '.hdf5':
        #Read the image
        f = h5py.File(image,'r+')
        #Get the dataset name from the h5 file
        dat_name = list(f.keys())[0]
        ###If the hdf5 is exported from ilastik fiji plugin, the dat_name will be 'data'
        #Get the image data
        image_loaded = np.array(f[dat_name])
        #Remove the first axis (ilastik convention)
        image_loaded = image_loaded.reshape((image_loaded.shape[1],image_loaded.shape[2],image_loaded.shape[3]))
        ###If the hdf5 is exported from ilastik fiji plugin, the order will need to be
        ###switched as above --> z_stack = np.swapaxes(z_stack,0,2) --> z_stack = np.swapaxes(z_stack,0,1)

    #Return the objects
    return image_loaded_z


def MaskZstack(mask_loaded,image,channel_names_loaded):
    """This function will extract the stats for each cell mask through each channel
    in the input image

    mask: Tiff image mask that represents the cells in your image. Must end with the word mask!!

    z_stack: Multichannel z stack image"""

    #Get the CellIDs for this dataset
    IDs = pd.DataFrame(MaskIDs(mask_loaded))
    #Iterate through the z stack to extract intensities
    list_of_chan = []
    #Get the z channel and the associated channel name from list of channel names
    for z in range(len(channel_names_loaded)):
        #Run the data Prep function
        image_loaded_z = PrepareData(image,z)
        #Use the above information to mask z stack
        list_of_chan.append(MaskChannel(mask_loaded,image_loaded_z))
    #Convert the channel names list and the list of intensity values to a dictionary and combine with CellIDs and XY
    dat = pd.concat([IDs,pd.DataFrame(dict(zip(channel_names_loaded,list_of_chan)))],axis=1)
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


def ExtractSingleCells(mask,image,channel_names,output):
    """Function for extracting single cell information from input
    path containing single-cell masks, z_stack path, and channel_names path."""

    #Create pathlib object for output
    output = Path(output)

    #Read the channels names
    channel_names_loaded = pd.read_csv(channel_names,header=None)
    #Add a column index for ease
    channel_names_loaded.columns = ["marker"]
    #Convert the channel names to a list
    channel_names_loaded = list(channel_names_loaded.marker.values)
    #Read the mask
    mask_loaded = skimage.io.imread(mask,plugin='tifffile')

    scdata_z = MaskZstack(mask_loaded,image,channel_names_loaded)
    #Write the singe cell data to a csv file using the image name

    im_full_name = os.path.basename(image)
    im_name = im_full_name.split('.')[0]
    scdata_z.to_csv(str(Path(os.path.join(str(output),str(im_name+".csv")))),index=False)


def MultiExtractSingleCells(mask,image,channel_names,output):
    """Function for iterating over a list of z_stacks and output locations to
    export single-cell data from image masks"""

    print("Extracting single-cell data for "+str(image)+'...')
    #Start memory monitoring --- next line can be commented
    tracemalloc.start()
    # ---

    #Run the ExtractSingleCells function for this image
    ExtractSingleCells(mask,image,channel_names,output)
    
    #Display memory monitoring --- next 3 lines can be commented
    current, peak = tracemalloc.get_traced_memory()
    print(f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB")
    tracemalloc.stop()
    # ---

    #Print update
    im_full_name = os.path.basename(image)
    im_name = im_full_name.split('.')[0]
    print("Finished "+str(im_name))


