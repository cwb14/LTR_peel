# Installing LTR_peel

The short route is in the [README](../README.md). This page has the details: the
container launcher, every install route, and what to do on a cluster without
internet.

LTR_peel drives a stack of external tools (GenomeTools, LTR_finder, MMseqs2,
HMMER, miniprot, BLAST+, minimap2, and others). You do not have to install any of
them: the container has them all, and needs neither `sudo` nor a package manager.

## Apptainer / Singularity (recommended)

The right choice on a shared cluster, where you cannot install Docker and do not
have root. The whole tool is one file you own, and you run that file:

```bash
apptainer pull ltr_peel.sif docker://ghcr.io/cwb14/ltr_peel:latest
./ltr_peel.sif --help
```

`:latest` is the newest release; for a paper, replace `latest` with the number of
the release you used.

`./ltr_peel.sif` **is** the `ltr_peel` command: the flags are the same, and it
runs as you, in your own directory:

```bash
./ltr_peel.sif --genome your_genome.fa --proteins related.pep.fa --threads 20
```

The individual stages work the same way: `./ltr_peel.sif ltr_peel-gff3 --help`.
Every example in the docs is written as plain `ltr_peel ...`; with a container,
read it as `./ltr_peel.sif ...`.

On a cluster you may need `module load apptainer` (or `module load singularity`)
first. `singularity` works identically: swap the command name. If your home
directory is small, send the download cache to scratch before pulling:

```bash
export APPTAINER_CACHEDIR=/your/scratch/.apptainer
```

### When it says your genome is not found

A container sees only the directories that were *mounted* into it. Apptainer
mounts your home directory and your working directory, which covers most runs
and misses two common ones:

* **your working directory is reached through a symbolic link.** `/home/you/data
  -> /scratch/you` is the usual layout on a cluster. Apptainer sees a path inside
  your home, mounts nothing extra, and then cannot follow the link from inside.
  It gives up, leaves you in your home directory, and every relative path in the
  command is suddenly wrong.
* **an input lives somewhere else entirely**, as in `--genome
  /scratch/genomes/hg38.fa` run from your home directory.

Fix either one on the spot: resolve the link, or name the directory to mount.

```bash
cd -P . && ./ltr_peel.sif --genome hg38.fa --threads 20

apptainer exec -B /scratch ltr_peel.sif \
    ltr_peel --genome /scratch/genomes/hg38.fa --threads 20
```

Or fix it once, with the launcher that ships inside the image. It works the
mounts out from your command and hands the rest to the container unchanged:

```bash
apptainer exec ltr_peel.sif cat /opt/ltr_peel/bin/ltr_peel-container > ltr_peel
chmod +x ltr_peel

./ltr_peel --genome /scratch/genomes/hg38.fa --threads 20   # mounts /scratch/genomes
```

It is a short shell script with no dependencies: read it before you run it if you
like. Keep it beside `ltr_peel.sif`, or put it on your `PATH` and tell it where
the image lives with `export LTR_PEEL_SIF=/shared/containers/ltr_peel.sif`.

## conda / mamba

If you would rather have the tools on your `PATH` than in a container, clone the
repository and build the environment from it:

```bash
git clone https://github.com/cwb14/LTR_peel.git
cd LTR_peel

mamba env create -n ltr_peel -f environment.yml   # add --yes to overwrite an existing one
mamba activate ltr_peel
pip install .

ltr_peel --help
```

Run these one at a time rather than pasting the block. If an `ltr_peel`
environment already exists, `mamba env create` stops and asks whether to
overwrite it, and a question that arrives in the middle of a pasted block is
answered by whatever line you pasted after it, not by you. For a clean rebuild,
remove the old environment first:

```bash
mamba env remove -n ltr_peel
```

`environment.yml` carries `git`, `make` and a C compiler on purpose, mainly for
TRF-mod and FastLTR, which are cloned *and* compiled, so without a compiler the
first real run stops on `cc: not found`. Pass `--no-trf` if you would rather not
have TRF-mod at all. TEsorter2 is cloned the same way but never compiled.
Kmer2LTR needs neither step, normally: it is an ordinary Python package,
resolved from `PATH` if it is already installed there, and only cloned (never
built) as a fallback when it is not. BLAST+ and minimap2 are in `environment.yml`:
every run uses minimap2, re-boundarying (on by default) needs BLAST+, and strand
recovery (off by default) needs both.

## Docker

For your own machine, where you have root.

```bash
docker pull ghcr.io/cwb14/ltr_peel:latest
```

`:latest` is the newest release; for a paper, replace `latest` with the number of
the release you used.

Docker has to be told which directory to mount, which directory to work in and
which user to be, so its commands are longer than the Apptainer equivalents:

```bash
docker run --rm -v "$PWD:/data" -w /data -u "$(id -u):$(id -g)" \
  ghcr.io/cwb14/ltr_peel:latest --help
```

Naming the command explicitly (`... ghcr.io/cwb14/ltr_peel:latest ltr_peel
--help`) works too, and is how you reach the other stages.

## pip

Python parts only; external tools not included.

```bash
pip install ltr_peel
```

You then have to supply the external binaries yourself; see
[`environment.yml`](../environment.yml) for the list. TEsorter2, TRF-mod, sdust
and FASTAN are fetched into the tools directory on first use. TRF-mod and sdust
are then compiled, and FastLTR is compiled from FASTAN with LTR_peel's own
`FastLTR.c`; TEsorter2 is not compiled. Kmer2LTR is an ordinary Python package:
`pip install` it yourself, or let LTR_peel fetch the same clone into the tools
directory, unbuilt either way. The container prepares all five during its own
build and never reaches the network.

## Offline clusters

Compute nodes without internet cannot fetch those helper tools. Fetch them once,
on a machine with internet, into a folder you name, and give every later run the
same folder:

```bash
ltr_peel --genome small.fa --max-rounds 1 --tools-dir /shared/ltr_peel_tools   # once, online
ltr_peel --genome big.fa --proteins prot.fa --tools-dir /shared/ltr_peel_tools # offline
```

A folder you name with `--tools-dir` (or the environment variable
`LTR_PEEL_TOOLS_DIR`) is reused and never removed. Without it, the folder is
`<genome>_tools` in the run folder (`merged_tools` for several genomes, `NAME_tools`
with `--out_prefix NAME`), and a successful run removes it again only if
that run created it and it holds only helper tools; a folder that already existed
when the run started, or that holds anything else, is kept. A fetch that fails with no internet
says so and points at `--tools-dir`. The container needs none of this.

## Upgrading from 1.0.1

The 1.0.1 README suggested an alias, `alias ltr_peel='apptainer exec ...'`. It is
no longer needed, and if you later install with conda it will shadow the real
command and fail with `apptainer: command not found`. Check with `type -a
ltr_peel`, and drop it with `unalias ltr_peel` (and remove it from your shell
startup file, if it is there).
