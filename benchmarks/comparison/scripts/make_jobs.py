#!/usr/bin/env python3
"""Write run directories (runs/<ds>/<tool>/cmd.sh) and a jobs.tsv for the queue.

Every tool runs with its documented default command (stated deviations: HiTE --te_type ltr; standalone LTR_Finder via
LTR_FINDER_parallel defaults), 16 threads, on bench/data/<ds>/genome.fa (normalized IDs). Each cmd.sh:
  - sets every *_NUM_THREADS to 16 (the server's /etc/environment caps them at 8 for all users),
  - puts temp files under bench/tmp/<ds>_<tool> (/ is nearly full),
  - prints 'STEP <name> start|end <epoch>' around each stage (component timing).
"""
import argparse
import csv
import os
import sys

B = os.environ["BENCH_ROOT"]
CONDA = 'source "$(conda info --base)/etc/profile.d/conda.sh"'
T = 16
TOOLS = ["ltr_peel", "ltrharvest", "ltr_finder", "ltr_retriever", "look4ltrs", "hirepeat", "hite", "ltrpred", "dante_ltr"]
# threads a tool actually uses (for the queue's budget); LTRharvest is single-threaded
USES = {"ltrharvest": 1}
PROT = {  # LTR_peel --proteins: closest RefSeq relative, never the target (sims: PrinTE's TAIR10.pep, the sims' genes)
    "rice": "proteins/GCF_000231095.2_ObraRS2_protein.faa.gz",
    "arab": "proteins/GCF_000004255.2_v.1.0_protein.faa.gz",
    "maize": "proteins/GCF_000003195.3_Sorghum_bicolor_NCBIv3_protein.faa.gz",
    "finch": "proteins/GCF_022539315.1_serCan2020_protein.faa.gz",
    "gorilla": "proteins/GCF_000001405.40_GRCh38.p14_protein.faa.gz",
    "asue": "proteins/GCF_000004255.2_v.1.0_protein.faa.gz",
    "pann": "proteins/GCF_019359855.2_Kyuss_2.0_protein.faa.gz",
    "smoke": "proteins/GCF_000004255.2_v.1.0_protein.faa.gz",
}
SIM_PROT = "proteins/TAIR10.pep.fa.gz"
ANIMAL = {"finch", "gorilla"}           # DANTE needs the Metazoa REXdb; everything else is plant
HARVEST_README = ("-minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 -similar 85 -vic 10 "
                  "-seed 20 -seqids yes")
LFP = f"perl {B}/build/LTR_FINDER_parallel/LTR_FINDER_parallel"
APPT = (f"conda activate apptainer; export APPTAINER_CACHEDIR={B}/tmp/apptainer_cache "
        f"APPTAINER_TMPDIR={B}/tmp/apptainer_tmp")


def step(name, body):
    return f'echo "STEP {name} start $(date +%s)"\n{body}\necho "STEP {name} end $(date +%s)"\n'


def tool_body(tool, ds, rundir):
    g = f"{ds}.fa"
    data = f"{B}/data/{ds}"
    if tool == "ltr_peel":
        prot = f"{B}/data/{SIM_PROT if ds.startswith('S') else PROT[ds]}"
        return step("ltr_peel", f"conda activate ltr_peel_bench\nexport LTR_PEEL_TOOLS_DIR={B}/build/ltr_peel_tools\n"
                    f"ltr_peel --genome {g} --proteins {prot} --threads {T} --detect-args \"ltr-timeout=300\" "
                    f"--run-sdust --no-plots --reboundary")
    if tool == "ltrharvest":
        return ("conda activate ltrr_bench\n"
                + step("suffixerator", f"gt suffixerator -db {g} -indexname {g} -tis -suf -lcp -des -ssp -sds -dna")
                + step("ltrharvest", f"gt ltrharvest -index {g} -gff3 ltrharvest.gff3 > ltrharvest.tabout"))
    if tool == "ltr_finder":
        return "conda activate ltrr_bench\n" + step("ltr_finder_parallel", f"{LFP} -seq {g} -threads {T}")
    if tool == "ltr_retriever":
        return ("conda activate ltrr_bench\n"
                + step("suffixerator", f"gt suffixerator -db {g} -indexname {g} -tis -suf -lcp -des -ssp -sds -dna")
                + step("ltrharvest", f"gt ltrharvest -index {g} {HARVEST_README} > {g}.harvest.scn")
                + step("ltr_finder_parallel", f"{LFP} -seq {g} -threads {T} -harvest_out -size 1000000 -time 300")
                + f"cat {g}.harvest.scn {g}.finder.combine.scn > {g}.rawLTR.scn\n"
                + step("ltr_retriever", f"{B}/build/LTR_retriever/LTR_retriever -genome {g} -inharvest {g}.rawLTR.scn "
                       f"-threads {T}"))
    if tool == "look4ltrs":
        return step("look4ltrs", f"mkdir -p {rundir}/out\n{B}/build/Look4LTRs/bin/look4ltrs --fasta {data}/split/ "
                    f"--out {rundir}/out/ --parallel {T}")
    if tool == "hirepeat":
        return step("hirepeat", f"conda activate hirepeat_bench\nexport PYTHONPATH={B}/build/HiRepEat\n"
                    f"python {B}/build/HiRepEat/engine/detect.py --genome {g} --threads {T} --output {rundir}")
    if tool == "hite":  # README's recommended conda route (the 3.3.3 image lacks the LtrDetector binary)
        return step("hite", f"set +u; conda activate hite_bench; set -u\n"
                    f"export PATH={B}/build/HiTE/module:{B}/build/HiTE/tools:$PATH\nmkdir -p {rundir}/out {rundir}/work\n"
                    f"python {B}/build/HiTE/main.py --genome {rundir}/{g} --out_dir {rundir}/out --thread {T} "
                    f"--work_dir {rundir}/work --te_type ltr")
    if tool == "ltrpred":
        return step("ltrpred", f"{APPT}\napptainer exec --cleanenv --writable-tmpfs -B {B} --pwd {rundir} {B}/containers/ltrpred.sif "
                    f"Rscript -e 'LTRpred::LTRpred(genome.file = \"{g}\", cores = {T})'")
    if tool == "dante_ltr":
        db = " -D Metazoa_v3.1" if ds in ANIMAL else ""
        return ("conda activate dante_bench\n"
                + step("dante", f"dante -q {g} -o dante.gff3 -c {T}{db}")
                + "conda activate dantel_bench\n"
                + step("dante_ltr", f"dante_ltr -g dante.gff3 -s {g} -o dante_ltr -c {T}")
                + step("dante_ltr_to_library", f"dante_ltr_to_library -g dante_ltr.gff3 -s {g} -o dante_ltr_lib -c {T}"))
    raise ValueError(tool)


def write_job(ds, tool, force=False):
    rd = f"{B}/runs/{ds}/{tool}"
    os.makedirs(rd, exist_ok=True)
    g = f"{B}/data/{ds}/genome.fa"
    if not os.path.exists(g):
        sys.exit(f"error: {g} missing (run prep_genome.py first)")
    for suf in ("", ".fai"):
        link = f"{rd}/{ds}.fa{suf}"
        if not os.path.lexists(link):
            os.symlink(g + suf, link)
    tmp = f"{B}/tmp/{ds}_{tool}"
    nt = " ".join(f"{v}={T}" for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                                       "NUMEXPR_NUM_THREADS", "NUMBA_NUM_THREADS", "BLIS_NUM_THREADS",
                                       "RAYON_NUM_THREADS"))
    cmd = (f"#!/usr/bin/env bash\n# {ds} / {tool}: generated by make_jobs.py\nset -euo pipefail\n{CONDA}\n"
           f"export {nt}\nexport TMPDIR={tmp}\nmkdir -p $TMPDIR\ncd {rd}\n" + tool_body(tool, ds, rd))
    path = f"{rd}/cmd.sh"
    if force or not os.path.exists(path):
        with open(path, "w") as f:
            f.write(cmd)
    return dict(job_id=f"{ds}__{tool}", ds=ds, tool=tool, threads=USES.get(tool, T), outdir=rd,
                size=os.path.getsize(g))


def main(argv=None):
    ap = argparse.ArgumentParser(description="Write runs/<ds>/<tool>/cmd.sh and a jobs.tsv.")
    ap.add_argument("--ds", nargs="+", required=True, help="datasets (bench/data/<ds>)")
    ap.add_argument("--tools", nargs="+", default=TOOLS, choices=TOOLS, help="tools (default: all 9)")
    ap.add_argument("--out", required=True, help="jobs.tsv to write")
    ap.add_argument("--force", action="store_true", help="overwrite existing cmd.sh")
    a = ap.parse_args(argv)
    rows = [write_job(ds, t, a.force) for ds in a.ds for t in a.tools]
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"make_jobs: {len(rows)} jobs -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
