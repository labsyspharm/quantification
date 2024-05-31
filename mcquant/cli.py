#Script for parsing command line arguments and running single-cell
#data extraction functions
#Joshua Hess
from . import ParseInput
from . import SingleCellDataExtraction

def main():
    #Parse the command line arguments
    args = ParseInput.ParseInputDataExtract()
    #Run the MultiExtractSingleCells function
    SingleCellDataExtraction.MultiExtractSingleCells(**args)

if __name__ == '__main__':
    main()
