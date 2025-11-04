# hash:sha256:d4dbea6e4441c3294e05cfa87afd0506c780f7793e095d40e7af81501d5a12eb
ARG REGISTRY_HOST
FROM $REGISTRY_HOST/codeocean/lightning-jupyterlab-py39:1.0

ARG DEBIAN_FRONTEND=noninteractive
ARG GIT_ASKPASS
COPY git-askpass /

RUN pip install -U --no-cache-dir \
    boto3 \
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
    xmltodict \
    zarr==2.18.2
