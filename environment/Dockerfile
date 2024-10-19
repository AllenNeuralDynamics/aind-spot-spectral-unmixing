# hash:sha256:13770f4295015504c4f5440afc1d78f9b575c4b3bcfca5bda77c15834db5c249
ARG REGISTRY_HOST
FROM $REGISTRY_HOST/codeocean/lightning-jupyterlab-py39:1.0

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
