#Script for parsing command line arguments and running single-cell
#data extraction functions
#Joshua Hess
import ParseInput
import SingleCellDataExtraction

#Parse the command line arguments
args = ParseInput.ParseInputDataExtract()

#Run the MultiExtractSingleCells function
SingleCellDataExtraction.MultiExtractSingleCells(**args)
