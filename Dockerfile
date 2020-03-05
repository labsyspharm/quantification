FROM python:3.6

RUN pip install scikit-image h5py pandas numpy pathlib

COPY CommandSingleCellExtraction.py ./app/CommandSingleCellExtraction.py
