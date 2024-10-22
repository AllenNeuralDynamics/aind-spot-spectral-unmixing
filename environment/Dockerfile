# hash:sha256:828d8c903c66d675a42be6e7c3b538d5d11e44a4241d3409f55d44d1b8ca8645
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
    numpy==2.0.2 \
    pandas==2.2.3 \
    scipy==1.13.1 \
    seaborn==0.13.2 \
    seaborn-image==0.10.0 \
    tqdm==4.66.5 \
    tqdm-joblib==0.0.4 \
    zarr==2.18.2

ADD "https://github.com/coder/code-server/releases/download/v4.9.0/code-server-4.9.0-linux-amd64.tar.gz" /.code-server/code-server.tar.gz
	
RUN cd /.code-server \
	&& tar -xvf code-server.tar.gz \
	&& rm code-server.tar.gz \
	&& ln -s /.code-server/code-server-4.9.0-linux-amd64/bin/code-server  /usr/bin/code-server
