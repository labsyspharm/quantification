FROM python:3.11

RUN pip install --no-cache-dir h5py pandas numpy pathlib 'scikit-image>=0.23.2' imagecodecs 

COPY . /app/
