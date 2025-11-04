# hash:sha256:eafcaad5be81a5bbbd32a661990cb9e322395c5c3412a78e79fa055dad1af153
ARG REGISTRY_HOST
FROM $REGISTRY_HOST/codeocean/lightning-jupyterlab-py39:1.0

ARG DEBIAN_FRONTEND=noninteractive
ARG GIT_ASKPASS
COPY git-askpass /

RUN pip install -U --no-cache-dir \
    boto3==1.40.65 \
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
    xmltodict==1.0.2 \
    zarr==2.18.2
