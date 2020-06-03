# Single cell quantification
Module for single-cell data extraction given a segmentation mask and multi-channel image. The CSV structure is aligned with histoCAT output.

**CommandSingleCellExtraction.py**:

* mask: Path to directory where masks are stored (Ex: ./segmentation/cellMask.tif)

* image: Path to image(s) for quantification.  (Ex: ./registration/*.h5) -> works with .h(df)5 or .tif(f)

* output: Path to output directory. (Ex: ./feature_extraction)

* channel_names: csv file containing the channel names for the z-stack (Ex: ./my_channels.csv)

# Run script
python CommandSingleCellExtraction.py --mask ./segmentation/cellMask.tif --image ./registration/Exemplar_001.h5  --output ./feature_extraction --channel_names ./my_channels.csv

# Main developer
Denis Schapiro (https://github.com/DenisSch)

Joshua Hess (https://github.com/JoshuaHess12)

Jeremy Muhlich (https://github.com/jmuhlich)

