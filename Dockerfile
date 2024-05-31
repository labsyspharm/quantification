FROM python:3.11

COPY . src

RUN pip install --no-cache-dir build \
    && python -m build -w src \
    && pip install --no-cache-dir src/dist/mcquant*.whl \
    && rm -r src
