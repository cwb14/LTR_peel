"""Single config module for the ltrkin_geno benchmark add-on (frozen values, paths, datasets)."""
import os
from pathlib import Path

BENCH = Path(os.environ["BENCH_ROOT"])
LG = BENCH / "ltrkin_geno"        # data root of the genotyping benchmark (calls/, pool/, truth/, runs/, report/)
# input data of the genotyping benchmark (data/Camara, bench/Camara/{sub,wga}); set LTRKIN_V2 to its location
V = Path(os.environ.get("LTRKIN_V2", str(LG / "v2")))

_CA = V / "data/Camara"
_SUB = V / "bench/Camara/sub"
_WGA = V / "bench/Camara/wga"
DATASETS = {
    "camara_hap2": dict(ref=_CA / "Camara_hap2.fa", other=_CA / "Camara_hap1.fa",
                        bam=_SUB / "hap2_d20_s1.bam", paf=_WGA / "hap2_on_hap1.paf.gz",
                        role="primary"),
    "camara_hap1": dict(ref=_CA / "Camara_hap1.fa", other=_CA / "Camara_hap2.fa",
                        bam=_SUB / "hap1_d20_s1.bam", paf=_WGA / "hap1_on_hap2.paf.gz",
                        role="replicate"),
}

TOOLS = ["ltr_peel", "ltrharvest", "ltr_finder", "ltr_retriever", "look4ltrs",
         "hirepeat", "hite", "ltrpred", "dante_ltr"]
SECONDARY = ["ltr_peel_hc"]
ALL_TOOLS = TOOLS + SECONDARY

PROT = "proteins/GCF_000004255.2_v.1.0_protein.faa.gz"
PY_LG = "python"          # interpreters resolved on PATH (conda envs ltrkin_geno and bench_tools)
PY_BENCH = "python"
SNAKEMAKE = "snakemake"

# Frozen values
NOM_ABS = 30
NOM_REL = 0.10
SPREAD = 0.20
SEG_FLANK = 1000
O_PAD = 200
MAX_LOCUS = 40_000
MIN_GAP = 50
NEAR_INDEL = 30
NEAR_BP = 100
MAX_DUP = 20
TSD_ASSUMED = 5
MIN_DUP = 4            # duplicated length < MIN_DUP: ends not pinned by a TSD (te_flag short_dup)
PSEUDO_LTR = 0.10
MIN_INTERNAL = 100     # true elements: median LTRs leaving < MIN_INTERNAL bp of internal region -> pseudo-LTRs
RO = 0.9
BIN_BP = (5, 20)
LEN_BINS = (2000, 5000, 10000)
N_BOOT = 1000
SEED = 1
BLOCK = 5_000_000
GQ_MIN = 20

# Run settings
QUEUE_BUDGET = 32
NICE = 10
LTRKIN_THREADS = 8
TRUTH_THREADS = 4      # hap_truth.py -t: only mappy's index build is threaded; per-element work runs on one core
