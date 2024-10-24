# hash:sha256:f387e3c6a3ff8099749e1a4639dee2f747f8972b39588e4c5f60bf2396050a16
ARG REGISTRY_HOST
FROM $REGISTRY_HOST/codeocean/lightning-jupyterlab-py39:1.0

ARG DEBIAN_FRONTEND=noninteractive

ARG GIT_ASKPASS
ARG GIT_ACCESS_TOKEN
COPY git-askpass /

RUN pip install -U --no-cache-dir \
    dask==2024.5.2 \
    joblib==1.4.2 \
    matplotlib==3.9.2 \
    numpy==1.26.4 \
    pandas==2.2.3 \
    scipy==1.13.1 \
    seaborn==0.13.2 \
    seaborn-image==0.10.0 \
    tqdm==4.66.5 \
    tqdm-joblib==0.0.4 \
    zarr==2.18.2
