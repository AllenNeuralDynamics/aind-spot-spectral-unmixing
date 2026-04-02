# hash:sha256:a9eb87555c5ad75743aafbc4596a7c8aca51cbb7aa08370a69f9e00e1716a77f
ARG REGISTRY_HOST
FROM $REGISTRY_HOST/codeocean/lightning-jupyterlab-py39:1.0

ARG DEBIAN_FRONTEND=noninteractive

ARG GIT_ASKPASS
ARG GIT_ACCESS_TOKEN
COPY git-askpass /

RUN pip install -U --no-cache-dir \
    dask==2023.5.0 \
    joblib==1.4.2 \
    matplotlib==3.7.5 \
    numpy==1.24.4 \
    pandas==2.0.3 \
    scipy==1.10.1 \
    seaborn==0.13.2 \
    seaborn-image==0.9.0 \
    tqdm==4.66.5 \
    tqdm-joblib==0.0.4 \
    zarr==2.16.1
