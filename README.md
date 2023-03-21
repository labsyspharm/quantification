# Single cell quantification
Module for single-cell data extraction given a segmentation mask and multi-channel image. The CSV structure is aligned with histoCAT output.

**CommandSingleCellExtraction.py**:

* `--masks` Paths to where masks are stored (Ex: ./segmentation/cellMask.tif) -> If multiple masks are selected the first mask will be used for spatial feature extraction but all will be quantified

* `--image` Path to image(s) for quantification.  (Ex: ./registration/*.h5) -> works with .h(df)5 or .tif(f)

* `--output` Path to output directory. (Ex: ./feature_extraction)

* `--channel_names` csv file containing the channel names for the z-stack (Ex: ./my_channels.csv)

* `--mask_props` Space separated list of additional metrics to be calculated for every mask.
    This is intended for metrics that depend only on the cell mask. If the metric depends
    on signal intensity, use `--intensity-props` instead.
    See list at https://scikit-image.org/docs/dev/api/skimage.measure.html#regionprops

* `--intensity_props` Space separated list of additional metrics to be calculated for every marker separately.
    By default only mean intensity is calculated.
    If the metric doesn't depend on signal intensity, use `--mask-props` instead.
    See list at https://scikit-image.org/docs/dev/api/skimage.measure.html#regionprops

    Currently, the following additional properties can be specified:

    * `--intensity_props gini_index` : The Gini index calculates a single number
    between 0 and 1, representing how unequal the signal is distributed in each region.
    See https://en.wikipedia.org/wiki/Gini_coefficient for more information.
    * `--intensity_props intensity_median` : Will calculate the median of intensity values per labeled object in the mask.
    * `--intensity_props intensity_sum` : Will calculate the sum of intensity values per labelled object in the mask. This can be useful if you want to count RNA molecules from FISH based images for example.

# Run script
`python CommandSingleCellExtraction.py --masks ./segmentation/cellMask.tif ./segmentation/membraneMask.tif --image ./registration/Exemplar_001.h5  --output ./feature_extraction --channel_names ./my_channels.csv`

# Main developer
Denis Schapiro (https://github.com/DenisSch)

Joshua Hess (https://github.com/JoshuaHess12)

Jeremy Muhlich (https://github.com/jmuhlich)