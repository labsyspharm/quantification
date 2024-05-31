#Functions for reading in single cell imaging data
#Joshua Hess

#Import necessary modules
import skimage.io
import h5py
import pandas as pd
import numpy as np
import os
import skimage.measure
import skimage.measure._regionprops
import tifffile

from pathlib import Path

#### Additional functions that can be specified by the user via intensity_props

## Function to calculate median intensity values per mask 
def intensity_median(mask, intensity):
    return np.median(intensity[mask])

## Function to sum intensity values (or in this case transcript counts)
def intensity_sum(mask, intensity):
    return np.sum(intensity[mask])

## Function to calculate the gini index: https://en.wikipedia.org/wiki/Gini_coefficient
def gini_index(mask, intensity):
    x = intensity[mask]
    sorted_x = np.sort(x)
    n = len(x)
    cumx = np.cumsum(sorted_x, dtype=float)
    return (n + 1 - 2 * np.sum(cumx) / cumx[-1]) / n

def MaskChannel(mask_loaded, image_loaded_z, intensity_props=["intensity_mean"]):
    """Function for quantifying a single channel image

    Returns a table with CellID according to the mask and the mean pixel intensity
    for the given channel for each cell"""
    standard_props = set(skimage.measure._regionprops.COL_DTYPES)
    # Look for regionprops in skimage
    builtin_props = set(intensity_props).intersection(standard_props)
    # Otherwise look for them in this module
    extra_props = set(intensity_props).difference(standard_props)
    dat = skimage.measure.regionprops_table(
        mask_loaded, image_loaded_z,
        properties = tuple(builtin_props),
        extra_properties = [globals()[n] for n in extra_props]
    )
    return dat


def MaskIDs(mask, mask_props=None):
    """This function will extract the CellIDs and the XY positions for each
    cell based on that cells centroid

    Returns a dictionary object"""

    all_mask_props = set(["label", "centroid", "area", "major_axis_length", "minor_axis_length", "eccentricity", "solidity", "extent", "orientation"])
    if mask_props is not None:
        all_mask_props = all_mask_props.union(mask_props)

    dat = skimage.measure.regionprops_table(
        mask,
        properties=all_mask_props
    )

    name_map = {
        "CellID": "label",
        "X_centroid": "centroid-1",
        "Y_centroid": "centroid-0",
        "Area": "area",
        "MajorAxisLength": "major_axis_length",
        "MinorAxisLength": "minor_axis_length",
        "Eccentricity": "eccentricity",
        "Solidity": "solidity",
        "Extent": "extent",
        "Orientation": "orientation",
    }
    for new_name, old_name in name_map.items():
        dat[new_name] = dat[old_name]
    for old_name in set(name_map.values()):
        del dat[old_name]

    return dat

def n_channels(image):
    """Returns the number of channel in the input image. Supports [OME]TIFF and HDF5."""

    image_path = Path(image)

    if image_path.suffix in ['.tiff', '.tif', '.btf']:
        s = tifffile.TiffFile(image).series[0]
        ndim = len(s.shape)
        if ndim == 2: return 1
        elif ndim == 3: return min(s.shape)
        else: raise Exception('mcquant supports only 2D/3D images.')

    elif image_path.suffix in ['.h5', '.hdf5']:
        f = h5py.File(image, 'r')
        dat_name = list(f.keys())[0]
        return f[dat_name].shape[3]

    else:
        raise Exception('mcquant currently supports [OME]TIFF and HDF5 formats only')

def PrepareData(image,z):
    """Function for preparing input for maskzstack function. Connecting function
    to use with mc micro ilastik pipeline"""

    image_path = Path(image)

    #Check to see if image tif(f)
    if image_path.suffix in ['.tiff', '.tif', '.btf']:
        image_loaded_z = tifffile.imread(image, key=z)

    #Check to see if image is hdf5
    elif image_path.suffix in ['.h5', '.hdf5']:
        #Read the image
        f = h5py.File(image,'r')
        #Get the dataset name from the h5 file
        dat_name = list(f.keys())[0]
        #Retrieve the z^th channel
        image_loaded_z = f[dat_name][0,:,:,z]

    else:
        raise Exception('mcquant currently supports [OME]TIFF and HDF5 formats only')

    #Return the objects
    return image_loaded_z


def MaskZstack(masks_loaded,image,channel_names_loaded, mask_props=None, intensity_props=["intensity_mean"]):
    """This function will extract the stats for each cell mask through each channel
    in the input image

    mask_loaded: dictionary containing Tiff masks that represents the cells in your image.

    z_stack: Multichannel z stack image"""

    #Get the names of the keys for the masks dictionary
    mask_names = list(masks_loaded.keys())

    #Create empty dictionary to store channel results per mask
    dict_of_chan = {m_name: [] for m_name in mask_names}
    #Get the z channel and the associated channel name from list of channel names
    for z in range(len(channel_names_loaded)):
        #Run the data Prep function
        image_loaded_z = PrepareData(image,z)

        #Iterate through number of masks to extract single cell data
        for nm in range(len(mask_names)):
            #Use the above information to mask z stack
            dict_of_chan[mask_names[nm]].append(
                MaskChannel(masks_loaded[mask_names[nm]],image_loaded_z, intensity_props=intensity_props)
            )
        #Print progress
        print("Finished "+str(z))

    # Column order according to histoCAT convention (Move xy position to end with spatial information)
    last_cols = (
        "X_centroid",
        "Y_centroid",
        "column_centroid",
        "row_centroid",
        "Area",
        "MajorAxisLength",
        "MinorAxisLength",
        "Eccentricity",
        "Solidity",
        "Extent",
        "Orientation",
    )
    def col_sort(x):
        if x == "CellID":
            return -2
        try:
            return last_cols.index(x)
        except ValueError:
            return -1

    #Iterate through the masks and format quantifications for each mask and property
    for nm in mask_names:
        mask_dict = {}
        # Mean intensity is default property, stored without suffix
        mask_dict.update(
            zip(channel_names_loaded, [x["intensity_mean"] for x in dict_of_chan[nm]])
        )
        # All other properties are suffixed with their names
        for prop_n in set(dict_of_chan[nm][0].keys()).difference(["intensity_mean"]):
            mask_dict.update(
                zip([f"{n}_{prop_n}" for n in channel_names_loaded], [x[prop_n] for x in dict_of_chan[nm]])
            )
        # Get the cell IDs and mask properties
        mask_properties = pd.DataFrame(MaskIDs(masks_loaded[nm], mask_props=mask_props))
        mask_dict.update(mask_properties)
        dict_of_chan[nm] = pd.DataFrame(mask_dict).reindex(columns=sorted(mask_dict.keys(), key=col_sort))

    # Return the dict of dataframes for each mask
    return dict_of_chan

def ExtractSingleCells(masks,image,channel_names,output, mask_props=None, intensity_props=["intensity_mean"]):
    """Function for extracting single cell information from input
    path containing single-cell masks, z_stack path, and channel_names path."""

    #Create pathlib object for output
    output = Path(output)

    #Read csv channel names
    channel_names_loaded = pd.read_csv(channel_names)
    #Check for the presence of `marker_name` column
    if 'marker_name' in channel_names_loaded:
        #Get the marker_name column if more than one column (CyCIF structure)
        channel_names_loaded_list = list(channel_names_loaded.marker_name)
    #Consider the old one-marker-per-line plain text format
    elif channel_names_loaded.shape[1] == 1:
        #re-read the csv file and add column name
        channel_names_loaded = pd.read_csv(channel_names, header = None)
        channel_names_loaded_list = list(channel_names_loaded.iloc[:,0])
    else:
        raise Exception('%s must contain the marker_name column'%channel_names)

    #Contrast against the number of markers in the image
    if len(channel_names_loaded_list) != n_channels(image):
        raise Exception("The number of channels in %s doesn't match the image"%channel_names)
    
    #Check for unique marker names -- create new list to store new names
    channel_names_loaded_checked = []
    for idx,val in enumerate(channel_names_loaded_list):
        #Check for unique value
        if channel_names_loaded_list.count(val) > 1:
            #If unique count greater than one, add suffix
            channel_names_loaded_checked.append(val + "_"+ str(channel_names_loaded_list[:idx].count(val) + 1))
        else:
            #Otherwise, leave channel name
            channel_names_loaded_checked.append(val)

    #Read the masks
    masks_loaded = {}
    #iterate through mask paths and read images to add to dictionary object
    for m in masks:
        m_full_name = os.path.basename(m)
        m_name = m_full_name.split('.')[0]
        masks_loaded.update({str(m_name):skimage.io.imread(m,plugin='tifffile')})

    scdata_z = MaskZstack(masks_loaded,image,channel_names_loaded_checked, mask_props=mask_props, intensity_props=intensity_props)
    #Write the singe cell data to a csv file using the image name

    # Determine the image name by cutting off its extension
    im_full_name = os.path.basename(image)
    im_tokens = im_full_name.split(os.extsep)
    if len(im_tokens) < 2: im_name = im_tokens[0]
    elif im_tokens[-2] == "ome": im_name = os.extsep.join(im_tokens[0:-2])
    else: im_name = os.extsep.join(im_tokens[0:-1])

    # iterate through each mask and export csv with mask name as suffix
    for k,v in scdata_z.items():
        # export the csv for this mask name
        scdata_z[k].to_csv(
                            str(Path(os.path.join(str(output),
                            str(im_name+"_{}"+".csv").format(k)))),
                            index=False
                            )


def MultiExtractSingleCells(masks,image,channel_names,output, mask_props=None, intensity_props=["intensity_mean"]):
    """Function for iterating over a list of z_stacks and output locations to
    export single-cell data from image masks"""

    print("Extracting single-cell data for "+str(image)+'...')

    #Run the ExtractSingleCells function for this image
    ExtractSingleCells(masks,image,channel_names,output, mask_props=mask_props, intensity_props=intensity_props)

    #Print update
    im_full_name = os.path.basename(image)
    im_name = im_full_name.split('.')[0]
    print("Finished "+str(im_name))
