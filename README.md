# Single-cell Data Extraction Functions - histoCAT template
Module for single-cell data extraction given a mask and multi-channel z-stack image.

**CommandNoiseRemoval.py**:
Script for thresholding and masking noise channel after Ilastik training. To run an example, modify the ExampleNoiseRemovalCommand.txt document to suit your needs and copy paste to terminal.
Note: csv file with channel names must have the same format as the csv file in the folder
ExampleNoiseRemovalCommand.txt arguments:

* mask_dir: Path to directory where masks are stored (Ex: ./my_mask_dir)

* z_stacks: Path to images being used to mask. Either single directory or number of directories=to number of images (Ex: ./25546ON.hdf5 ./26531POST.hdf5)

* output: Path to output directory. (Ex: ./my_output_dir)

* channel_names: csv file containing the channel names for the z-stack (Ex: ./my_channels.csv)

* suffix: Series of suffixes applied to the oringal images upon exporting from ilastik/noise removal. Assumes that you only added suffixes to your original images upon masking and removing noise (Ex: orignal = my_image.hdf5 -> ilastik -> my_image_probabilities.tif -> noise removal -> my_image_probabilities_RemoveNoise.tif -> segmentation -> my_image_probabilities_RemoveNoise_mask.tif...suffix here is _probabilities_RemoveNoise_mask)
