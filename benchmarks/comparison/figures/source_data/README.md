# Figure source data

Each folder holds the values one figure draws: every element, call, curve point and score on it, and nothing else.
With these tables the figure redraws exactly (same pixels, same `_numbers.tsv`) without the benchmark dataset.

```bash
cd benchmarks/comparison
python figures/<stem>.py --from-source-data      # redraw from source_data/<stem>/ only
BENCH_ROOT=/path/to/bench python figures/<stem>.py -v   # rebuild from the full dataset; rewrites source_data/<stem>/
```

Run in the `bench_tools` environment (`envs/`); `figS4_ltrkin.py` uses `longreads`. A full run writes the tables,
reads them back and draws only from what it read, so both routes give the same figure. The README figures
(`benchmark_summary`, `nesting_depth`) write their PNG to `docs/images/`.

**Full dataset.** Every run, call set, mask, library and score behind these tables (about 640 GB) is on the Ou lab
server (asc-oulab) at `/data2/chris/poa_synthesis/benchmarking_ltrrts/bench`; the genotyping inputs of test 8 (the
*Cardamine amara* haplotypes, reads and alignments) are at `/data2/chris/poa_synthesis/discretionary_analysis3/v2`
(set `LTRKIN_V2` to it).

**Format.** Tab-separated, one header row; tables over 1 MB are gzipped. Tool keys are those of the benchmark
(`ltr_peel`, `ltr_retriever`, ...). Genome coordinates are 0-based half-open (BED style); `nesting_depth` gives
positions within each FASTA record, 1-based closed. File names start with the panel they feed (`A_`, `B_`, ...);
`scalars.tsv` holds single values and `legend_*.tsv` values printed only in the legend.

| Figure | Tables |
|---|---|
| `fig_benchmark` (Fig. 1) | A: true elements, each tool's calls, their TSDs and domains, which element each tool finds; B: scorecard; C, D: recall curves and recall-halved points (coding, protein-less); E: masking scores; F: curated families recovered; G: curated coverage along every *A. thaliana* locus and its moved span, age classes; H: genotyping by caller and class |
| `figS1_nonauto` | A: the rice example and its ends; B: every simulated element and which tools found it; C: recall by length; D: scorecard and basis; E: support of protein-less calls in real genomes |
| `figS2_masking` | A: genes, curated hits, calls and masked intervals in the rice window; B: masked bases; C: masking scores; missing cells |
| `figS3_library_age` | A: the two example copies, their curated hits and domains; B: curated coverage in 50 bins along every locus of four genomes and its moved span (380,330 rows, the largest table at 8.4 MB); legend shares |
| `figS4_ltrkin` | A: example locus, flank alignment blocks, spanning reads, calls and their features; B: junction sequences; C, D: genotyping by caller and class |
| `figS5_resources` | every run's wall time, CPU time and memory; genome sizes |
| `benchmark_summary` (README) | A: true elements, calls, TSDs, domains; B: scorecard; C: coding recall curves; D: masking scores |
| `nesting_depth` (README) | A, B: the three nested elements, their masked runs and domains; C: the letters drawn at the join |
