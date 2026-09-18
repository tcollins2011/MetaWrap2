# MetaWrap2 container image.
#
# Ships the lightweight MetaWrap2 core (python3 + biopython + the metawrap2
# package). Heavy per-module tools live in the per-module conda environments
# under envs/ and can be layered on top as needed (see the build args below),
# which keeps the base image small for users who only need part of the pipeline.
#
# Build the core image:
#   docker build -t metawrap2:core -f containers/Dockerfile .
#
# Build with one or more module environments layered in, e.g. binning:
#   docker build --build-arg MODULE_ENVS="envs/binning.yaml" \
#     -t metawrap2:binning -f containers/Dockerfile .

FROM mambaorg/micromamba:1.5-jammy

# Optional space-separated list of env YAMLs to install into the base env.
ARG MODULE_ENVS=""

USER root
WORKDIR /opt/metawrap2

# Copy the repository into the image.
COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/metawrap2

# Core runtime: python + biopython + pip, then the metawrap2 package itself.
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        "python>=3.8" "biopython>=1.79" pip && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1
# Installs the metawrap2 package and its `metawrap2` console command on PATH.
RUN pip install --no-cache-dir /opt/metawrap2

# Optionally layer in per-module tool environments.
RUN if [ -n "$MODULE_ENVS" ]; then \
        for env_file in $MODULE_ENVS; do \
            echo "Installing $env_file into base env" && \
            micromamba install -y -n base -f "/opt/metawrap2/$env_file" ; \
        done && \
        micromamba clean --all --yes ; \
    fi

USER $MAMBA_USER
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["metawrap2", "--help"]
