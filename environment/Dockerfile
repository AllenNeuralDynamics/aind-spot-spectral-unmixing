# hash:sha256:4c5228bd3ec9848ec8c6b2043b3308c53eb5c77d1a64253efeea75319683faa8
ARG REGISTRY_HOST
FROM $REGISTRY_HOST/codeocean/jupyterlab:3.0.9-miniconda4.9.2-python3.8-ubuntu20.04

ARG DEBIAN_FRONTEND=noninteractive

ARG GIT_ASKPASS
ARG GIT_ACCESS_TOKEN
COPY git-askpass /

RUN pip install -U --no-cache-dir \
    dask \
    joblib \
    matplotlib \
    numpy \
    pandas \
    scipy \
    seaborn \
    seaborn-image \
    tqdm \
    tqdm-joblib \
    zarr
