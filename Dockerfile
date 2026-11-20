# syntax=docker/dockerfile:1
#
# LTR_peel — hermetic image.
#
# The detector normally fetches and builds five helper repos into --tools-dir on
# first use. That is convenient on a workstation and wrong in a container: it
# makes the image's behaviour depend on what GitHub served that morning. Here
# they are cloned at pinned commits during the build, and LTR_PEEL_TOOLS_DIR
# points every run at the result, so a container never reaches the network.
#
#   docker build -t ltr_peel:1.0.1 .
#   docker run --rm -v "$PWD:/data" -w /data ltr_peel:1.0.1 \
#       --genome genome.fa --proteins prot.fa --threads 8
#
# The ENTRYPOINT at the bottom is what lets the arguments be written bare like
# that, and it is also what Apptainer turns into the SIF's runscript, so a pulled
# image can be run as `./ltr_peel.sif --genome genome.fa`. See bin/entrypoint.sh.

# ---------------------------------------------------------------- build stage
FROM mambaorg/micromamba:1.5.10-jammy AS build

USER root
RUN apt-get update \
 && apt-get install -y --no-install-recommends build-essential git zlib1g-dev ca-certificates \
 && rm -rf /var/lib/apt/lists/*

# Pinned so a rebuild six months from now produces the same tools.
ARG KMER2LTR_REF=2c7e7f8d43697939dece93029e554418821d494a
ARG TESORTER2_REF=cc9b195ca25dbe63555d45af2b3e9380b697080b
ARG TRFMOD_REF=3e891db310124f7e5f7a630a1c006650be9d1f3a
ARG SDUST_REF=89c42cb41ba598e9cfa07c2ef99ae8c08f769b3e
ARG FASTAN_REF=843d0ed28b93b05ad428995d986312f02468afcd

# LTR_peel's adaptation of FastLTR, compiled against the pinned FASTAN sources
# below; ltr_peel.fastltr builds it the same way outside the container. The last
# line checks the build by the marker in its usage text. Run bare, FastLTR prints
# that text and exits 1, and this image's build shell sets pipefail, so the exit
# status is dropped (`|| true`) and grep alone decides.
COPY src/ltr_peel/fastltr/FastLTR.c /tmp/FastLTR.c

ENV TOOLS=/opt/ltr_peel/tools
RUN mkdir -p "$TOOLS" \
 && git clone https://github.com/cwb14/Kmer2LTR.git "$TOOLS/Kmer2LTR" \
 && git -C "$TOOLS/Kmer2LTR" checkout --quiet "$KMER2LTR_REF" \
 && git clone --branch feat/minimap2 https://github.com/cwb14/TEsorter2.git "$TOOLS/TEsorter2" \
 && git -C "$TOOLS/TEsorter2" checkout --quiet "$TESORTER2_REF" \
 && git clone https://github.com/lh3/TRF-mod "$TOOLS/TRF-mod" \
 && git -C "$TOOLS/TRF-mod" checkout --quiet "$TRFMOD_REF" \
 && make -C "$TOOLS/TRF-mod" -f compile.mak \
 && git clone https://github.com/lh3/sdust.git "$TOOLS/sdust" \
 && git -C "$TOOLS/sdust" checkout --quiet "$SDUST_REF" \
 && make -C "$TOOLS/sdust" \
 && git clone https://github.com/thegenemyers/FASTAN.git "$TOOLS/FASTAN" \
 && git -C "$TOOLS/FASTAN" checkout --quiet "$FASTAN_REF" \
 && cp /tmp/FastLTR.c "$TOOLS/FASTAN/FastLTR.c" \
 && cd "$TOOLS/FASTAN" \
 && gcc -O3 -fno-strict-aliasing -o FastLTR FastLTR.c alncode.c align.c ANO.c GDB.c gene_core.c ONElib.c -lm -lz -lpthread \
 && cd / \
 && find "$TOOLS" -name .git -type d -prune -exec rm -rf {} + \
 && test -x "$TOOLS/TRF-mod/trf-mod" && test -x "$TOOLS/sdust/sdust" \
 && { "$TOOLS/FASTAN/FastLTR" 2>&1 || true; } | grep -q "LTR_peel adaptation"

# --------------------------------------------------------------- runtime stage
FROM mambaorg/micromamba:1.5.10-jammy

LABEL org.opencontainers.image.title="LTR_peel" \
      org.opencontainers.image.description="Iterative detection of nested LTR retrotransposons" \
      org.opencontainers.image.source="https://github.com/cwb14/LTR_peel" \
      org.opencontainers.image.licenses="GPL-3.0-or-later" \
      org.opencontainers.image.version="1.0.1"

USER root
# procps: Nextflow's task tracer needs `ps`. awscli: AWS Batch stages S3 inputs
# with the CLI from *inside* the container (conf/awsbatch.config points at it).
RUN apt-get update \
 && apt-get install -y --no-install-recommends procps ca-certificates \
 && rm -rf /var/lib/apt/lists/*

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
 && micromamba install -y -n base -c conda-forge awscli \
 && micromamba clean --all --yes \
 && rm /tmp/environment.yml

ENV PATH=/opt/conda/bin:$PATH \
    LTR_PEEL_TOOLS_DIR=/opt/ltr_peel/tools \
    MPLBACKEND=Agg \
    PYTHONNOUSERSITE=1

COPY --from=build /opt/ltr_peel/tools /opt/ltr_peel/tools

COPY . /opt/ltr_peel/src
RUN pip install --no-cache-dir /opt/ltr_peel/src \
 && install -D -m 0755 /opt/ltr_peel/src/bin/entrypoint.sh \
        /opt/ltr_peel/bin/entrypoint \
 && install -D -m 0755 /opt/ltr_peel/src/bin/ltr_peel-container \
        /opt/ltr_peel/bin/ltr_peel-container \
 && rm -rf /opt/ltr_peel/src

# Fail the build, not the first user, if the image is wired up wrong.
RUN ltr_peel --help > /dev/null \
 && python -c "import ltr_peel, matplotlib, numpy; print('ltr_peel', ltr_peel.__version__)" \
 && for t in gt ltr_finder mmseqs hmmsearch blastn minimap2 miniprot mafft; do \
        command -v "$t" > /dev/null || { echo "missing: $t" >&2; exit 1; }; \
    done \
 && test -f "$LTR_PEEL_TOOLS_DIR/Kmer2LTR/src/kmer2ltr/cli.py" \
 && test -x "$LTR_PEEL_TOOLS_DIR/TRF-mod/trf-mod" \
 && test -x "$LTR_PEEL_TOOLS_DIR/FASTAN/FastLTR" \
 && test -x /opt/ltr_peel/bin/ltr_peel-container

# All three ways the entry point is reached, asserted here so a dispatch
# regression fails the build rather than the first person to run the image:
# bare flags, an explicit stage command, and an arbitrary command (which is how
# Nextflow drives the image).
RUN /opt/ltr_peel/bin/entrypoint --help > /dev/null \
 && /opt/ltr_peel/bin/entrypoint ltr_peel-gff3 --help > /dev/null \
 && /opt/ltr_peel/bin/entrypoint /bin/sh -c 'exit 0'

USER $MAMBA_USER
WORKDIR /data

# An empty ENTRYPOINT is what broke `./ltr_peel.sif --genome x.fa` in 1.0.1:
# Apptainer builds the runscript from ENTRYPOINT and CMD, and with only a CMD it
# *replaces* the command with the user's arguments, so the first flag was run as
# a program. With an ENTRYPOINT, Apptainer prepends it to the arguments instead.
# `docker run IMG ltr_peel ...` and `apptainer exec IMG ltr_peel ...` are
# unaffected: the dispatcher runs a first argument that names a command, and
# `exec` never reaches the entry point at all.
ENTRYPOINT ["/opt/ltr_peel/bin/entrypoint"]
CMD ["--help"]
