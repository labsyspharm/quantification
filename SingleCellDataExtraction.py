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
import shutil
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
    #Remove excess region prop information from memory
    dat = 0
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

    #Remove excess region props from memory
    dat = 0

    #Form a dataframe from the lists
    IDs = {"CellID": labels, "X_position": xcoords, "Y_position": ycoords,"Area":area,\
        "MajorAxisLength":major_axis_length,"MinorAxisLength":minor_axis_length,\
        "Eccentricity":eccentricity,"Solidity":solidity,"Extent":extent,"Orientation":orientation}
    #Remove excess information from memory
    labels = 0
    xcoords = 0
    ycoords = 0
    area = 0
    minor_axis_length=0
    major_axis_length=0
    eccentricity = 0
    solidity = 0
    extent=0
    orientation=0
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


def MaskZstack(mask_loaded,image,channel_names_loaded,output):
    """This function will extract the stats for each cell mask through each channel
    in the input image
    mask: Tiff image mask that represents the cells in your image. Must end with the word mask!!
    z_stack: Multichannel z stack image"""

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
    IDs = pd.DataFrame(MaskIDs(mask_loaded))
    #Create a name for the cell ID file
    exIDs = Path(os.path.join(tmp_dir,im_name.stem)+"_cellIDs.csv")
    #Write out the cellIDs to csv
    IDs.to_csv(exIDs,index=False)
    #Remove the cell ID information from memory
    IDs = 0
    #Update the list of filenames
    fnames.append(exIDs)

    #Iterate through the z stack to extract intensities
    #list_of_chan = []
    tracemalloc.start()
    #Iterate through the z stack to extract intensities (Parallelize?)
    for z in range(len(channel_names_loaded)):
        #Create a name for this channel
        exChan = Path(os.path.join(tmp_dir,im_name.stem.replace('.ome',''))+"_"+str(channel_names_loaded[z])+".csv")
        #Run the data Prep function
        image_loaded_z = PrepareData(image,z)
        #Use the above information to mask z stack
        #list_of_chan.append(MaskChannel(mask_loaded,image_loaded_z))
        z_chan = pd.DataFrame(MaskChannel(mask_loaded,image_loaded_z),columns=[str(channel_names_loaded[z])])
        #Write out the z_channel to a csv
        z_chan.to_csv(exChan,index=False)
        #Remove the loaded image and csv frmo memory
        image_loaded_z = 0
        z_chan = 0
        #Update the list of filenames
        fnames.append(exChan)

        #Display memory monitoring --- next 3 lines can be commented
        current, peak = tracemalloc.get_traced_memory()
        print(f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB")
        print("Finished "+str(z))
    #Change back to the home directory
    os.chdir(home_dir)
    #Return a list containing the file names that are individually exported
    return fnames, tmp_dir


def ConcatenateData(fnames,tmp_dir):
    """Function for concatenating csv files that are exported individually"""

    #Read each file(CellIDs should be first in the list)
    IDs = pd.read_csv(fnames[0])
    #Read all of the channels in a list, should be in order if not parallel
    z_chans = pd.concat([pd.read_csv(f) for f in fnames[1:]],axis=1)
    #Concatenate the cellID and spatial information with mean intensity per channel
    dat = pd.concat([IDs,z_chans],axis=1)

    #Remove the excess data from memory
    IDs = 0
    z_chans = 0
    #Delete all of the temporary files that were created
    shutil.rmtree(tmp_dir)

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


    #Read the channels names
    channel_names_loaded = pd.read_csv(channel_names,header=None)
    #Add a column index for ease
    channel_names_loaded.columns = ["marker"]
    #Convert the channel names to a list
    channel_names_loaded = list(channel_names_loaded.marker.values)
    #Read the mask
    mask_loaded = skimage.io.imread(mask,plugin='tifffile')
    #Write single cell information for each channel individually
    scdata_fnames_z, tmp_dir = MaskZstack(mask_loaded,image,channel_names_loaded,output)

    #Concatenate all of the channel files and cell IDs
    scdata_z = ConcatenateData(scdata_fnames_z,tmp_dir)

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
