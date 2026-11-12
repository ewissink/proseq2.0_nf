FROM mambaorg/micromamba:1.5.10

LABEL description="proseq2.0-nf pipeline tools (built from env/environment.yml)"

USER root
# procps provides `ps`, which Nextflow needs for task metrics
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps \
    && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER env/environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
    && micromamba clean -a -y

# Put the env on PATH so tools resolve without `micromamba run`
ENV PATH=/opt/conda/bin:$PATH
