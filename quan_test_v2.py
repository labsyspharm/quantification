#Functions for reading in single cell imaging data
#Joshua Hess

#Import necessary modules
import h5py
# TODO Create a reader for hdf5 images
import pandas as pd
import numpy as np
import skimage.measure as measure

import tifffile
import warnings
import textwrap
import pathlib


PROP_VALS = measure._regionprops.PROP_VALS

NAME_MAP = {
    'label': 'CellID',
    'centroid-1': 'X_centroid',
    'centroid-0': 'Y_centroid',
    'area': 'Area',
    'major_axis_length': 'MajorAxisLength',
    'minor_axis_length': 'MinorAxisLength',
    'eccentricity': 'Eccentricity',
    'solidity': 'Solidity',
    'extent': 'Extent',
    'orientation': 'Orientation'
}

def gini_index(mask, intensity):
    x = intensity[mask]
    sorted_x = np.sort(x)
    n = len(x)
    cumx = np.cumsum(sorted_x, dtype=float)
    return (n + 1 - 2 * np.sum(cumx) / cumx[-1]) / n

def median_intensity(mask, intensity):
    return np.median(intensity[mask])

EXTRA_PROPS = {
    'gini_index': gini_index,
    'median_intensity': median_intensity
}

def quantify_channel(
    mask, intensity_img, 
    intensity_props=None, channel_name=None,
    preprocess_func=None, preprocess_func_kwargs=None,
    return_img=False
):   
    if intensity_props is None:
        intensity_props = []
    
    intensity_props = intensity_props[:] + ['label', 'mean_intensity']
    # Look for regionprops in skimage
    builtin_props = set(intensity_props).intersection(PROP_VALS)
    # Otherwise look for them in this module
    extra_props = [
        EXTRA_PROPS[p] for p in
        set(intensity_props).difference(PROP_VALS)
    ]
    if preprocess_func is not None:
        if preprocess_func_kwargs is None:
            preprocess_func_kwargs = {}
        intensity_img = preprocess_func(intensity_img, **preprocess_func_kwargs)
    props_table = measure.regionprops_table(
        mask, intensity_img,
        properties = tuple(builtin_props),
        extra_properties = extra_props
    )
    if channel_name is not None:
        def format_name(prop_name):
            if len(extra_props) == 0: return channel_name
            else: return f"{channel_name}_{prop_name}"
        return {
            k if k in NAME_MAP.keys() else format_name(k): v
            for k, v in props_table.items()
        }
    if return_img:
        return props_table, intensity_img
    else:
        return props_table


def quantify_mask(mask, mask_props=None):   
    all_mask_props = set([
        "label", "centroid", "area",
        "major_axis_length", "minor_axis_length",
        "eccentricity", "solidity", "extent", "orientation"
    ])
    if mask_props is not None:
        all_mask_props = all_mask_props.union(mask_props)

    dat = measure.regionprops_table(
        mask,
        properties=all_mask_props
    )
    return dat


def load_marker_csv(csv_path):
    csv_path = pathlib.Path(csv_path)
    assert csv_path.suffix == '.csv', (
        f"{csv_path} is not a CSV file."
    )
    marker_name_df = pd.read_csv(csv_path)
    if 'marker_name' not in marker_name_df.columns:
        print(
            f"'marker_name' not in {csv_path.name} header\n"
            f"Assuming legacy format, first column is used as marker names"
        )
        marker_name_df = pd.read_csv(
            csv_path, header=None, 
            usecols=[0], names=['marker_name']
        )
    has_duplicates = marker_name_df.duplicated(keep=False)
    name_suffix = (
        marker_name_df.loc[has_duplicates]
            .groupby('marker_name')
            .cumcount()
            .map(lambda x: f"_{x + 1}")
    )
    marker_name_df.loc[has_duplicates, 'marker_name'] += name_suffix
    return marker_name_df['marker_name'].to_list()


def validate_masks(mask_paths):
    for p in mask_paths:
        assert pathlib.Path(p).exists(), (
            f"{p} does not exist"
        )
    mask_shapes = []
    for p in mask_paths:
        with tifffile.TiffFile(pathlib.Path(p)) as tiff:
            mask_shapes.append(tiff.series[0].shape)
    assert len(set(mask_shapes)) == 1, (
        f"Masks must be the same shape\n"
        '\n'.join([   f"{p} - {s}"
            for p, s in zip(mask_paths, mask_shapes)
        ])
    )
    shape = mask_shapes[0]
    ndim = len(shape)
    assert ndim == 2, (
        f"Only 2D masks are supported. Got a {ndim}D mask of shape {shape}"
    )
    return shape


def validate_img(img_path):
    suffixes = pathlib.Path(img_path).suffixes
    assert '.ome' in suffixes
    assert ('.tif' in suffixes) ^ ('.tiff' in suffixes)
    with tifffile.TiffFile(img_path) as tiff:
        img_shape = tiff.series[0].shape
    ndim = len(img_shape)
    assert (ndim == 2) ^ (ndim == 3), (
        f"Only 2D/3D images are supported. Got a {ndim}D image of shape {img_shape}"
    )
    if len(img_shape) == 2:
        img_shape = (1, *img_shape) 
    return img_shape


def format_mask_table(mask_table):
    mask_table = pd.DataFrame(mask_table)
    mask_table.rename(columns=NAME_MAP, inplace=True)
    mask_table.set_index('CellID', inplace=True)
    return mask_table


def write_table(
    df, output_dir, 
    mask_path, img_path=None, 
    prefix='', suffix='', flat=False
):
    if not img_path:
        img_path = ''
    mask_name, img_name = [
        pathlib.Path(p).stem.replace('.ome', '').replace('.', '_')
        for p in [mask_path, img_path]
    ]
    
    output_dir = pathlib.Path(output_dir)
    if not flat:
        output_dir = output_dir / img_name
    
    output_dir.mkdir(exist_ok=True, parents=True)
    
    output_filename = f"{prefix}{img_name}_{mask_name}{suffix}.csv"
    output_path = output_dir / output_filename
    df.to_csv(output_path)
    return output_path


def validate_props(props):
    if props is None:
        return
    for p in props:
        assert (p in PROP_VALS) ^ (p in EXTRA_PROPS.keys()), (
            f"{p} is not a valid property. Available properties are "
            f"[{', '.join(PROP_VALS)}] and [{', '.join(EXTRA_PROPS.keys())}]"
        )
    return
