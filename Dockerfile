FROM python:3.9

RUN pip install h5py pandas numpy pathlib scikit-image>=0.18.0 imagecodecs

COPY . /app/
