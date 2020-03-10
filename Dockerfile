FROM python:3.6

RUN pip install scikit-image h5py pandas numpy pathlib

COPY . /app/
