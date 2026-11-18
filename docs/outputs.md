# LTRquest output reference

What a run does by default, what it keeps, and the full format of every file it
writes. For what the tool is and how to run it, see the [README](../README.md);
to install it, [install.md](install.md); for the Nextflow pipeline,
[nextflow.md](nextflow.md); for changing the code, [developers.md](developers.md).

1. [What runs by default and how to turn it off](#1-what-runs-by-default-and-how-to-turn-it-off)
2. [Run folder clean-up and --keep-intermediates](#2-run-folder-clean-up-and---keep-intermediates)
3. [Capped chunks and FastLTR](#3-capped-chunks-and-fastltr)
4. [Primary outputs](#4-primary-outputs-depthn_clean_ltrtsvfa)
5. [Multiple genomes](#5-multiple-genomes-shared-family-names)
6. [GFF3](#6-gff3-annotation)
7. [Strand recovery](#7-strand-recovery---strand-recovery)
8. [Re-boundarying](#8-re-boundarying)
9. [Plots](#9-plots-prefix_plots)
10. [Benchmarks](#10-benchmarks)
11. [Settings not in the short help](#11-settings-not-in-the-short-help)

Every flag below is listed by `ltrquest --help-all`.

---

## 1. What runs by default and how to turn it off

Three things run unless you switch them off. Each was chosen because it removes
false calls or recovers real ones on most genomes; the switches are for a genome
where one of them hurts, or a machine that lacks a tool.

| Default | What it does | Turn off |
|---|---|---|
| Per-chunk caller cap, 300 s | Each caller may spend 300 s on one 500 kb chunk, so a satellite or rDNA array cannot stall the run. See [section 3](#3-capped-chunks-and-fastltr). | `--ltr-timeout 0` (no cap) or another number of seconds |
| sdust low-complexity filter | Drops candidates made mostly of low-complexity sequence (poly-A, simple repeats) early. Together with `--max-dust-frac` (default 0.55; raise toward 0.80 for AT-rich or fungal genomes) and `--sdust-args` (default `-w 64 -t 15`). | `--no-sdust` |
| Re-boundarying | Moves the ends of truncated calls onto the ends that carry a target-site duplication. [Section 8](#8-re-boundarying). Needs BLAST+. | `--no-reboundary` |

Strand recovery is off by default: on large genomes it takes about as long as
detection and reaches few elements ([section 7](#7-strand-recovery---strand-recovery));
`--strand-recovery sensitive` turns it on.

On a machine without BLAST+ the run stops at startup; `--no-reboundary` runs
without it.

The cap is wall-clock time, not work done. A borderline chunk that finishes in
280 s on an idle machine can be capped on a loaded one, so two runs of the same
genome can differ in which chunks were capped. FastLTR searches capped chunks
(section 3), which is why this rarely changes the final set.

Changing a default changes what is found. Each detected genome records its
settings (the cap, the sdust switch, strand recovery and others), and a later run
in the same folder with different ones stops rather than mixing them; pass
`--redetect` to detect again ([section 5](#5-multiple-genomes-shared-family-names)).
A run folder made by an earlier version, before these defaults, therefore needs
`--redetect` too.

Other on-by-default behavior has its own switch: `--no-tandem-filter` (next
section), `--no-trf` (TRF tandem-repeat masking), `--no-plots`, and
`--no-fastltr` (same as `--fastltr off`).

### Tandem arrays

Complete (telomere-to-telomere) assemblies include the rDNA arrays and
centromeric satellites, and an LTR caller pairs any two copies of a repeat unit
as "LTRs". The `_clean_` tables ([section 4](#4-primary-outputs-depthn_clean_ltrtsvfa))
drop such a call along with the false-positive families when its repeat carries
on past **both** of its ends: at least half of the 15-mers in the 1 kb outside
each LTR recur in the 2 kb inside the element next to the other LTR. A genuine
element's flanks are unrelated target-site sequence and score near 0. Two
elements sharing a middle LTR each continue on one side only, so they are kept;
three or more in head-to-tail tandem are treated as an array. The unfiltered
`depth<N>_ltr.tsv` tables keep these calls, `<prefix>_fpcheck.tandem.tsv` lists
every call's two scores and verdict (both kept only with `--keep-intermediates`),
and they never count toward the FP fraction that triggers masking. Detection is
untouched, so a reused genome needs no `--redetect` to gain or lose the filter.
`--no-tandem-filter` keeps these calls in the `_clean_` outputs.

## 2. Run folder clean-up and --keep-intermediates

A run builds a lot on the way. When it succeeds it keeps only the results:

- `<prefix>_all_depth_LTR_cleaned.gff3`
- `<prefix>_depth<N>_clean_ltr.tsv` and `.fa`
- `<prefix>_plots/`

and removes the rest: the round tables and work directories, the unfiltered depth
tables and FASTAs, the rejoined, strand and re-boundarying tables
(`_rejoined.tsv`, `_strand_recovery.tsv`, `_reboundary.tsv`), the false-positive
check files, the pooled cluster files, the protein GFF3
(`_all_depth_protein_LTR_cleaned.gff3`), the reuse record (`.detect.json`), and a
helper-tools folder this run created. Only names the run itself gives its working
files, directly inside the run folder, are touched; your inputs never are.

`--keep-intermediates` keeps everything. You need it for:

- **adding a genome to this run later** ([section 5](#5-multiple-genomes-shared-family-names)),
  which reuses the detection records of the genomes already done;
- **re-plotting** with `ltrquest-plots`, whose per-family summary pages and TEGV
  gene track read the cluster files and gene tracks in the work directories;
- **`ltrquest-reboundary --posthoc`** on a finished run;
- **inspecting a step**: the per-round work directories hold every caller's
  candidates and the tables of what each filter dropped.

Everything is also kept automatically, with a message, when plotting, strand
recovery or re-boundarying failed, so the step can be redone. Nothing is removed
for a genome that has no GFF3 or `_clean_` tables (when the false-positive step
never ran on degenerate input), since its unfiltered tables are then all it has.

The helper tools (TEsorter2, TRF-mod, sdust, FASTAN/FastLTR, Kmer2LTR) are fetched
into a tools folder. One you name with `--tools-dir DIR` (or `LTRQUEST_TOOLS_DIR`)
is reused by later runs and never removed. The default, a `*_tools` folder in the
run folder (`<genome>_tools`, as in `Athal_tair10_chr2_tools`; `merged_tools` for
several genomes; `NAME_tools` with `--out_prefix NAME`), is removed at the end only
if this run created it and it holds only helper tools. A folder that already
existed when the run started, or that holds anything else, is kept. For machines
without internet, fill a named folder once online
([install.md](install.md#offline-clusters)).

Nothing is written beside your genome: LTRquest's own indexes and scratch files live
in the run folder and go with the clean-up. Helper tools may still put small
temporary files in the system temp folder (`$TMPDIR`, usually `/tmp`).

## 3. Capped chunks and FastLTR

Each round cuts the genome into 500 kb chunks and runs LTRharvest and LTR_FINDER
on every chunk, plus, where `--fastltr` says so, a third caller: FastLTR
(LTRquest's adaptation of Gene Myers' FastLTR, from
[FASTAN](https://github.com/thegenemyers/FASTAN), the detector behind
[HiRepEat](https://github.com/AnantMaheshwari/HiRepEat)). Their candidates are
pooled before Kmer2LTR, and every later stage treats them alike.

On a satellite, an rDNA array or a long telomere, LTRharvest and LTR_FINDER can
run for tens of minutes on one chunk. `--ltr-timeout N` (default 300; 0 = no cap)
stops each such call after N seconds, but neither tool writes anything until it
finishes, so a chunk it is cut short on loses **all** of that caller's
candidates, not only the array's. The cap is wall-clock time: a borderline chunk
can be capped on a loaded machine and not on an idle one.

FastLTR crosses an array in seconds, so it still searches those chunks. It counts
8-mer matches along the diagonals of 32 kb windows and aligns only the promising
ones; LTRquest's version marks each tandem alignment so it crosses an array once
instead of once per spacing, and it searches the same LTR spacings as LTRharvest
(100-20,000 bp). `--fastltr rescue` (the default) runs it only on chunks a cap
fired on; `all` runs it on every chunk (more elements, lower precision); `off`
never (`--no-fastltr` is the same). `--fastltr-args` (default `-d100 -D20000 -w8`)
passes options to it. Cap and FastLTR together:

```bash
ltrquest --genome t2t.fa --proteins prot.fa --threads 20 --ltr-timeout 30
```

The detector log names each capped chunk (`[WARN] ltrharvest timeout ...`) and
counts them per caller at the end of step 5. In `--detect-args`, `ltr-timeout=N`
overrides `--ltr-timeout`.

FastLTR reports every non-overlapping direct-repeat pair, where LTRharvest keeps
only the best of overlapping ones, so inside a degenerate tandem array it would
pair each unit with the next. A FastLTR pair is dropped when either copy recurs
one spacing further out (an LTR-RT has two copies, an array more), and step 5
counts those too. That test reaches 20 kb; an array whose unit is longer (human
45S rDNA, 43 kb) is recognized neither by it nor by the tandem filter
([section 1](#tandem-arrays)), and in human rDNA rescue calls copies of a
structure inside the unit (`LTR/unknown`, low confidence), as LTRharvest and
LTR_FINDER do when they finish such a chunk. Each round's work directory keeps
every caller's candidates in LTRharvest's tabular layout, chromosome last:
`<prefix>_r<N>.work/<prefix>_r<N>.{ltrharvest,ltrfinder,fastltr}.stitched.scn`
(kept with `--keep-intermediates`).

Against uncapped runs, on the chunks a 30 s cap cut both callers short on:

| | *A. thaliana* T2T | rice T2T-NIP |
|---|---|---|
| curated-supported elements recovered, cap alone | 0 of 97 | 1 of 76 |
| ... cap + `--fastltr rescue` | **66 of 97** | **65 of 76** |
| elements kept where no cap fired | 97.7% | 99.6% |
| wall time: cap / cap + rescue / uncapped | 0.33 / 0.37 / 1.48 h | 0.81 / 0.70 / 1.26 h |

Inside rice's CENH3 centromeres the cap alone leaves 97 elements, rescue 150 and
the uncapped run 155. Of 401 LTR-RTs planted into real CEN180 and rDNA arrays, a
30 s cap finds 26%, with rescue 88%, uncapped 92%. Benchmarks and how to rerun
them: [benchmarks/fastltr/](../benchmarks/fastltr/README.md).

## 4. Primary outputs: `depth<N>_clean_ltr.{tsv,fa}`

Elements are bucketed by nesting depth. `N` =
how many layers of LTR-RT are nested *inside* the element:

```
depth0:  [===============]                     no LTR-RT inside ("un-nested")
depth1:  [=====[depth0]=====]                  one LTR-RT nested inside
depth2:  [===[depth1 [depth0] ]===]            two layers nested inside
...
```

For the test run you'd get, e.g.:

```
Athal_tair10_chr2_LTRs_depth0_clean_ltr.tsv   ← table of un-nested elements
Athal_tair10_chr2_LTRs_depth0_clean_ltr.fa    ← their sequences
Athal_tair10_chr2_LTRs_depth1_clean_ltr.tsv   ← single-nested elements (ie, one LTR-RT nested inside)
Athal_tair10_chr2_LTRs_depth1_clean_ltr.fa
...
```

### 4.1 `depth<N>_clean_ltr.tsv` format

Tab-separated, one row per element, 34 named columns (`#`-prefixed header
line). Columns 1–29 are Kmer2LTR's own output: Kmer2LTR runs ahead of
classification, so its boundary calls, divergence figures and TSDs are the
round's single source of truth, with one correction (see *Canonical ends*
below). LTRquest appends the last five — `strand`, `family`, `confidence`,
`domains`, `nest_status` — once classification and clustering have run.

| # | Column | Meaning |
|---|---|---|
| 1 | `seq_id` | Element id: `chrom:start-end#Order/Superfamily/Clade`, e.g. `Chr2:102001-110500#LTR/Gypsy/Reina` |
| 2 | `seq_len` | Length of the record as supplied |
| 3 | `status` | `pass`, `weak_pair`, `k2p_undefined`, `no_pair`, `too_short`, or `all_ambiguous`. Only `pass` rows reach LTRquest's output |
| 4 | `ltr5_start` | 1-based inclusive, relative to the record. Always `1` here — see below |
| 5 | `ltr5_end` | 1-based inclusive, relative to the record: last bp of the 5′ LTR |
| 6 | `ltr3_start` | 1-based inclusive, relative to the record: first bp of the 3′ LTR |
| 7 | `ltr3_end` | 1-based inclusive, relative to the record. Always equals `seq_len` here — see below |
| 8 | `ltr5_len` | Length of the 5′ LTR |
| 9 | `ltr3_len` | Length of the 3′ LTR |
| 10 | `flank5_len` | bp of over-extension the detector added on the 5′ side. Always `0` here, by construction — see below |
| 11 | `flank3_len` | bp of over-extension the detector added on the 3′ side. Always `0` here, by construction — see below |
| 12 | `aln_len` | Length of the LTR-pair alignment, **including** gap columns |
| 13 | `n_sites` | Ungapped, unambiguous columns — the denominator for `identity`, `p_dist` and `k2p` |
| 14 | `n_ts` | Transitions |
| 15 | `n_tv` | Transversions |
| 16 | `n_gapcols` | Alignment columns with a gap on either side |
| 17 | `identity` | `n_match / n_sites` |
| 18 | `p_dist` | `(n_ts + n_tv) / n_sites` |
| 19 | `k2p` | Kimura 2-parameter distance. `NA` when saturated, never clamped |
| 20 | `k2p_se` | Standard error of `k2p` |
| 21 | `bitscore` | Global alignment score of the two LTRs under a fixed generic model |
| 22 | `flank_margin_bits` | How decisively the boundary was called. `NA` when the discovery core already reached both termini |
| 23 | `cigar` | Extended CIGAR of the LTR-pair alignment; query = 5′ LTR |
| 24 | `motif` | The two terminal dinucleotides, lowercased, e.g. `tg...ca`. Read off the called boundary, never used to find it |
| 25 | `k2p_time` | `round(k2p / (2 × mutation_rate))`, years since insertion; set by `--mutation-rate` (default `3e-8`) |
| 26 | `orientation` | Whether the record is stored reverse-complemented relative to its own header interval, as written to `{prefix}_ltr.fa`: `-` means the library holds the reverse complement, `+` means it holds the interval as read off the genome. A storage fact, not biology, and always known — see `strand` below |
| 27 | `tsd` | Target-site duplication at the called boundary. `.` = searched and absent; `NA` = could not be searched. On TG..CA ends it is the exact 4–6 bp duplication at those ends |
| 28 | `tsd_offset` | `d5,d3` — how far each boundary had to move for `tsd` to appear, positive meaning *into* the element. `0,0` on canonical ends. `NA` whenever `tsd` carries no sequence, whether that is `.` or `NA` |
| 29 | `tsd_input` | The same measurement at the record's termini as originally supplied, before re-bounding |
| 30 | `strand` | The element's biological orientation, resolved by LTRquest from classification, domain order and pass-2 homology, and from strand recovery when it is on: `+`, `-`, or `.` if none of those resolved it |
| 31 | `family` | Family id from the pooled consensus-LTR clustering, `<prefix>_fam00001` … — see [section 5](#5-multiple-genomes-shared-family-names) |
| 32 | `confidence` | Hallmarks of an insertion at the reported ends: `high` = TG..CA termini and an exact TSD (`tsd` with `tsd_offset` `0,0`), `medium` = one of the two, `low` = neither — see below |
| 33 | `domains` | Protein domains with genomic coordinates: `GENE\|clade@start-end;...`, or `.` |
| 34 | `nest_status` | `nest-outer:chrom:s-e` (that element is inside me) and/or `nest-inner:chrom:s-e` (I am inside that element), `;`-joined, or `.` if un-nested |

**Canonical ends.** An LTR-RT is integrated into a staggered cut, which leaves
the same few target bases on both sides (the TSD; 5 bp in most elements, 4–6 in
general), and integrase leaves TG at its 5′ end and CA at its 3′ end. Kmer2LTR
reports a TSD it finds 1 bp off the called ends but does not move them, and its
first hit is not always the element's ends (`X1..X5|TG…CA|X1..X5` called 1 bp
wide carries an exact 4-bp TSD at the wrong ends). So right after Kmer2LTR, each
round moves every element whose ends are not TG..CA with an exact TSD onto the
nearest ends that are, if any lie within 2 bp of each end (fewest bases moved
first), and Kmer2LTR re-measures it on the pair it was called with. Two calls
that reach the same ends are one element; the second is removed.

**Overlapping calls.** Real insertions are disjoint or nested: a nested element
lies strictly inside its host with its own LTR pair. Any other overlap — two
calls crossing, one inside another but sharing an LTR, the same span twice —
means at most one of them is right, and only the call with the stronger
structural evidence is kept: first the hallmarks at its own ends (both, then a
TSD alone, then TG..CA alone, then neither), then the LTR pair's `bitscore`, then
the number of distinct protein domains. The rule is applied within each round and
again when the rounds are pooled. Each round's dropped calls, and what they lost to,
are in `<prefix>_r<N>.work/<prefix>_r<N>.conflicts.tsv` (kept with `--keep-intermediates`).

**Hosts called with their tenant.** A round keeps a call only if TEsorter2
classifies it. A host found in a later round than its tenant is classified with
the tenant masked, but one called in the same round is classified with the tenant
still inside, where the tenant's protein domains can classify it. So a nest-outer
every one of whose HMM hits lies inside its tenants, and that has no TSD to admit
it, is withdrawn, and the round's overlapping calls are resolved again without it.
A later round, with the tenant masked, can still call it and judge it on its own
sequence. Each round's withdrawn hosts, with their tenants, are in
`<prefix>_r<N>.work/<prefix>_r<N>.withdrawn_outers.tsv` (kept with `--keep-intermediates`).

**`confidence`.** A real insertion is born with both hallmarks. Its TSD copies and
termini then mutate like the rest of the element, so an old element can lose
them, but a young one keeps them: in a young element (`identity` ≥ 0.99, say)
`low` marks a call that is very likely wrong, while in an old one it is
expected. Species whose elements do not end TG..CA reach at most `medium`; filter
on `tsd` there instead.

Three things this table will burn you on if you skim it:

- **`.` and `NA` are not the same value.** `tsd` (and `tsd_input`, the same
  measurement taken at the record's original termini) write `.` when the
  search ran and found no duplication — a real negative result — and `NA`
  when the search had no chance to run at all, for instance because the
  element's contig could not be matched back to the reference genome. Treat
  the two as one "missing" bucket and you will misjudge every TSD call:
  presence-vs-absence and ran-vs-could-not-run are different findings, and
  only one of them says anything about the boundary.
- **`strand` and `orientation` answer different questions.** `orientation`
  is a storage fact: whether the record in `{prefix}_ltr.fa` is the reverse
  complement of its own header interval (`-`) or reads that interval forward
  off the genome (`+`) — the call the library-writing step makes so that
  family clustering sees one orientation per family. `strand` is a biological
  claim: which way the element is transcribed, decided later by LTRquest from
  classification, domain order and pass-2 homology. The two are independent —
  an element stored forward (`orientation=+`) can carry either strand, and
  vice versa. `strand` can land on `.` when none of those sources resolve
  it; `orientation` never does, because it records a decision the writer
  always makes. With strand recovery the two do line up on the elements
  recovery calls, and only those: in the `_clean_` FASTAs and tables it flips
  their records to coding sense and restates `orientation` to match, exactly as
  the classification step already does for the elements TEsorter2 strands. Every
  other row is untouched.
- **Every row satisfies `ltr5_start == 1` and `ltr3_end == seq_len`.** That is
  what "Kmer2LTR-bounded" means: the record has already been cut down to the
  called LTR pair, so both flanks report `0` and there is nothing upstream or
  downstream of the element left in the record. This is not merely intended —
  `assert_bounded` checks it against every row at the end of each round and
  raises before the round can finish if a single one fails.

### 4.2 `depth<N>_clean_ltr.fa` format — nested regions are masked

Each record is the **full-length LTR-RT**, with every
nested inner element **hard-masked** by an IUPAC letter chosen by the *inner
element's own depth*:

| Inner element's depth | Mask letter |
|---|---|
| 0 | `N` |
| 1 | `R` |
| 2 | `D` |
| 3+ | `Y`, `S`, `W`, `K`, `M`, `B`, `H` |

So a `depth1` record looks like `ACGT...NNNN...ACGT` (its depth0 insert masked
with `N`), and a `depth2` record looks like `ACGT...RRR..NNN..RRR...ACGT`
(its depth1 child masked `R`, whose own depth0 child is masked `N`). 

**To remove the masked inner sequence** and keep only the outer element's sequence:

```bash
sed '/^>/!s/[^ATCGatcg]//g' Athal_tair10_chr2_LTRs_depth1_clean_ltr.fa > depth1_removed_nested.fa
```

Line lengths become uneven and some lines empty, which `samtools faidx` refuses;
rewrap with `seqkit seq -w 60` (seqkit is on bioconda). Alternatively, without
seqkit, this awk command writes each sequence on one line, which `samtools faidx`
indexes as it is:

```bash
awk '/^>/{if(s)print s; s=""; print; next}{gsub(/[^ACGTacgt]/,""); s=s $0}END{if(s)print s}' \
  Athal_tair10_chr2_LTRs_depth1_clean_ltr.fa > depth1_removed_nested.fa
```

`depth0` records contain no masked nest-ins, so this is only needed for depth ≥ 1.

### 4.3 Hosts rebuilt around a tenant in one LTR: `<prefix>_rejoined.tsv`

An element that inserts into another element's LTR splits it: `[1..p] TENANT
[p+1..L]`. LTRharvest and LTR_FINDER pair the two near-side fragments, so on
their own they report the host cut short at the tenant, with no TSD, and the
tenant as a free element beside it. When two calls abut (within 20 bp),
`ltrquest.reconcile` asks whether the bases beyond one continue the other's own
LTR from the junction on: identity ≥ 0.70, E ≤ 1e-5, the alignment starting at
both junctions. Where one side and only one does, the host is rebuilt to the end
of that continuation, and the tenant becomes its `nest-inner`. Such a host's row
differs from others in four ways:

- the LTR interval holding the tenant (`ltr5_start..ltr5_end` or
  `ltr3_start..ltr3_end`) spans it, and the record masks it like any other inner;
- `ltr5_len` / `ltr3_len` count LTR bases only, not the tenant inside;
- `tsd`, `tsd_offset`, `tsd_input` and `motif` are read at the rebuilt ends;
- the divergence columns (`aln_len` … `k2p_time`) are Kmer2LTR's measurement on
  the fragment pair it was called with: the same LTR pair, fewer sites.

A pair is left as called when both sides continue, when the continuation would
run through another call, when an element would be in two rejoins, or when a
rebuilt LTR would be longer than 7 kb, the longest LTR the callers report. Past
that length the "continuation" is the next unit of a tandem array. Such a host
is never a re-boundarying target. The sidecar (kept only with `--keep-intermediates`) has one row per rebuilt host, and
the strand stages use it to carry the call made on the detected fragment to the
new key:

| Column | Meaning |
|---|---|
| `old_key` | `chrom:start-end` as detected |
| `new_key` | `chrom:start-end` as rebuilt |
| `tenant` | the element inside the split LTR |
| `ltr` | which LTR holds it, in genome orientation: `left` or `right` |
| `remainder_bp` | bases of the split LTR found beyond the tenant |
| `identity`, `evalue` | of that continuation |
| `tsd` | at the rebuilt ends (`.` none, `NA` not searchable) |

Reconciling without `--genome` rebuilds nothing and writes the header alone.

## 5. Multiple genomes: shared family names

Pass several genomes to give closely-related species **one family vocabulary**,
which is what makes between-species family comparisons meaningful:

```bash
ltrquest \
  --genome   ../tests/data/Athal_tair10_chr2.fa.gz Alyrata_chr2.fa \
  --proteins ../tests/data/Athal.pep.gz \
  --threads  40
```

Each genome is detected on its own, then **all** detected elements are pooled
for a single clustering pass, so `family` means the same thing everywhere:

```
Athal_tair10_chr2_LTRs_depth0_clean_ltr.tsv     family = merged_fam00001
Alyrata_chr2_LTRs_depth0_clean_ltr.tsv          family = merged_fam00001
Athal_tair10_chr2_LTRs_all_depth_LTR_cleaned.gff3
Alyrata_chr2_LTRs_all_depth_LTR_cleaned.gff3
merged_all_ltr.consensus_id0.75_cluster.tsv     ← the one shared table
```

False-positive families are also called over the pooled set rather than one
species at a time, so a repeat that looks convincing in isolation but wrong
across species gets caught.

Naming: every genome keeps its own `<basename>_LTRs` prefix. `--out_prefix`
names only the shared pool — `--out_prefix Arabidopsis` gives
`Arabidopsis_fam00001` and `Arabidopsis_all_ltr.*`, default `merged`. With a
single genome the two are the same thing, so nothing changes.

One `--proteins` file serves every genome.

**Adding a genome later.** The earlier run must have used `--keep-intermediates`:
the reuse needs the detection records and tables a default run removes
([section 2](#2-run-folder-clean-up-and---keep-intermediates)). Re-run in the same
directory with the longer list (`--genome A.fa B.fa C.fa`). Genomes already
detected there are reused (a single-genome re-run too); only `C` is detected.
Families, `_clean_` tables and GFF3s are then rebuilt for all three, so family
ids can be renumbered, but the families are the same ones a run of all three from
scratch gives: reused genomes enter the pool exactly as detection wrote them.
Each detected genome leaves a `<prefix>.detect.json`; the run stops if a reused
genome's FASTA, `--proteins` or detection settings no longer match it. Delete
that file to redo one genome, or pass `--redetect` to redo all. A run folder made
by an earlier version with strand recovery on holds re-oriented unfiltered FASTAs;
add a genome to one of those with `--redetect`.

> **Sequence IDs must be unique across genomes.** Two files both calling a
> chromosome `Chr2` are rejected before any work starts: pooled clustering keys
> elements on `chrom:start-end`, so a shared ID would cross-assign families
> *and* cross-purge real elements between species. Rename first:
>
> ```bash
> awk '/^>/{sub(/^>/,">Aly_")}1' Alyrata.fa > Alyrata.renamed.fa
> ```

## 6. GFF3 annotation

Two files, pooled across all depths and built from the FP- and
tandem-array-purged `_clean_ltr.tsv` set (falling back to the raw set, with a warning, if the FP
stage never ran):

| Output | Contents |
|---|---|
| `<prefix>_all_depth_LTR_cleaned.gff3` | The LTR-RTs. Kept |
| `<prefix>_all_depth_protein_LTR_cleaned.gff3` | The same, plus every miniprot protein alignment. Only written when the run had `--proteins`, and removed at the end unless `--keep-intermediates` |

Each element is one block:

```
chr  LTRquest  LTR_retrotransposon   13031  17307  .  -  .  ID=…_LTRRT_00001;Name=chr:13031-17307;
                                                            classification=LTR/Gypsy/Tekay;superfamily=Gypsy;clade=Tekay;
                                                            family=…_fam00001;family_size=5;depth=0;
                                                            K2P_d=0.026305;K2P_T=438418;strand_source=tesorter
chr  LTRquest  long_terminal_repeat  13031  14048  .  -  .  ID=…_LTRRT_00001.lLTR;Parent=…_LTRRT_00001
chr  LTRquest  long_terminal_repeat  16290  17307  .  -  .  ID=…_LTRRT_00001.rLTR;Parent=…_LTRRT_00001
chr  LTRquest  protein_match         14482  14712  .  -  .  ID=…_LTRRT_00001.CHD.1;Parent=…;Name=CHD;clade=Tekay
###
```

Notes:

- Nested elements are **flat top-level features**, not children of their host:
  GFF3 `Parent` means part-of, and a nested LTR-RT is not a part of the element
  it landed in. Use the `depth` and `nest_status` attributes instead.
- Blocks are coordinate-sorted, and lines within a block are too, but
  **overlapping blocks are not interleaved** — so a nested element's block
  follows its host's in full.
- GFF3 attribute names are their own vocabulary, independent of the element
  table's column names: `K2P_d`/`K2P_T` above read from `k2p`/`k2p_time` in
  section 4.1. Keeping the attribute names stable regardless of what the
  table's own columns are called is what lets an existing GFF3 consumer keep
  working unmodified.
- **Column 9 is meant to be read, not just parsed**, so it carries no
  percent-escapes. GFF3 reserves `;` `=` `&` `,` inside an attribute value,
  and the two LTRquest attributes that are *lists* — `tsd_offset` (`d5,d3`)
  and `nest_status` — would otherwise arrive as `0%2C0` and
  `nest-outer:…%3Bnest-inner:…`. Both are instead written with `,`, which is
  GFF3's own separator for multiple values within one attribute: `tsd_offset=0,0`,
  `nest_status=nest-outer:chr1:200-300,nest-inner:chr1:50-900`. Spec-legal and
  legible. Note this differs from the element table, where `nest_status` joins
  on `;` (section 4.1).
- Escaping is still applied as a backstop to text LTRquest does not control —
  a sequence named `scaf;1` cannot be written literally without corrupting the
  file. If that ever happens the run says so on stderr rather than letting you
  discover a `%3B` by eye.
- `strand_source` names which tier of the cascade decided column 7:
  `tesorter`, `domain_order`, `pass2`, `homology` or `ppt` from strand
  recovery, `none` where the strand is `.`, or `table` if the strand in the
  element table is not one the cascade would reach — a hand-edited table, say.
- `confidence` (`high`, `medium`, `low`) is the table's column of the same name,
  computed from the same `tsd`, `tsd_offset` and `motif` (section 4.1).

## 7. Strand recovery (`--strand-recovery`)

Off by default; `--strand-recovery sensitive` (or `balanced`, `conservative`) turns it on, with the
polypurine-tract fallback. It is slow on large genomes and reaches few elements: on three grass genomes it
took about as long as detection and gave a strand to 10-22% of the elements still unstranded, because most
of them sit in families with no stranded member to copy from. A meaningful fraction of elements come out of the strand
cascade at `.`, most often because TEsorter2 could not classify them and they
carry too few protein domains to read an order off. Recovery transfers
orientation onto those from the elements that already carry a strand: the
unstranded element, in genome-forward orientation, is aligned against the
stranded ones held in coding sense, and the alignment's orientation is the
query's own genomic strand.

Four independent views run — dc-megablast and minimap2, each on the whole
element and on the internal region between the LTRs. A view calls only when
every one of its own alignments agrees, and any disagreement between views
vetoes the locus outright. The preset sets how much agreement a surviving locus
needs:

| Preset | Views | Agreement | Recall | Disagreement with the cascade |
|---|---|---|---|---|
| `conservative` | the two dc-megablast views | 1 | 0.58 | 0.19% |
| `balanced` | all four | 2 | 0.79 | 0.70% |
| `sensitive` | all four | 1, unless another contradicts | 0.84 | 1.06% |

Measured by leave-one-out against LTRquest's own tesorter strand, pooled over
five Brassicaceae genomes and 3616 labeled elements. Disagreement tracks how
well LTRquest could label the element in the first place: 0.37% where five or
more protein domains support the label, 2.54% where only one does.

The polypurine-tract fallback, on whenever strand recovery runs (`--no-strand-recovery-ppt` turns it
off), covers loci homology cannot reach. It scores both strand hypotheses on the
two windows flanking the internal region and calls only when one wins by a margin
and clears an absolute floor. Measured at 93% on real genomes, which is lower than
the homology calls, so its calls are labeled `strand_source=ppt` and can be
filtered out.

What changes in the outputs:

| Output | Change |
|---|---|
| `<prefix>_strand_recovery.tsv` | New. `locus`, `strand`, `source`, `evidence`, one row per element recovery was asked about — declines included, so comparing two presets is a diff. Kept only with `--keep-intermediates` |
| `depth<N>_clean_ltr.tsv` | `strand` filled on the recovered rows; `orientation` restated on those it flipped |
| `depth<N>_clean_ltr.fa` | Recovered minus elements stored in coding sense, like TEsorter2-called ones |
| both GFF3s | `strand` on the recovered elements, and `strand_source=homology` or `ppt` |

It re-orients the `_clean_` depth FASTAs and tables and nothing else. The
unfiltered `<prefix>_depth<N>_ltr.{tsv,fa}`, `<prefix>_r<N>_ltr.fa` (the
per-round library) and `<prefix>_all_ltr.fa` (the pooled clustering input) all
describe the state before annotation ran and keep detection's orientation. That
is why adding a genome later gives the same families as a run from scratch, and
why re-clustering them after a recovered run reproduces the original answer
rather than the recovered one.

Two further things it deliberately does not do. Family clustering runs before
annotation, so it sees the pre-flip sequences: recovery does not re-cluster and
does not change any family label. And the structure PDFs re-infer strand from
domain order rather than reading the `strand` column, so they will not show a
recovered call.

Needs `blast+` on `PATH` for every preset, and `minimap2` as well for
`balanced` and `sensitive`. Both are in `environment.yml` and the image, and the
wrapper checks for them while parsing its arguments, so a missing one fails the
run immediately rather than after every expensive stage has already succeeded.

## 8. Re-boundarying

LTRharvest and LTR_FINDER extend an LTR pair outward from a seed and stop at the
first obstacle near an LTR end (an indel or a mutation-dense patch), and Kmer2LTR can
only trim. Such calls end short of the element and have no TSD (target-site
duplication, the short direct repeat an insertion leaves on both sides) at their
called ends. A few calls run a few bases long, with the same result.

Re-boundarying is on by default (`--no-reboundary` turns it off; `--reboundary` is still accepted). It, or `ltrquest-reboundary --posthoc` on a finished run, moves the ends
of every call whose `tsd` is `.` or `NA`. A call with any TSD, including one Kmer2LTR
found only 1 bp off (`tsd_offset` other than `0,0`), is never moved.

How the new ends are chosen:

- **Templates** are calls whose TSD sits exactly at the called ends (`tsd_offset`
  `0,0`): independent evidence that both ends are right. They are pooled over every
  genome, family, strand and nesting depth. A call is not a template if either LTR
  contains `N` or hosts a nested element. At most 100 per family are searched (the
  youngest half by K2P, plus a draw seeded by the family name).
- **Search.** One BLASTN (E ≤ 1e-10, both strands) of each call's two LTRs finds its
  5 best templates (`--templates`).
- **Placement.** On each side, the template's LTR is aligned end to end around the
  call's LTR, using the template end that the template's own TSD verifies. If the
  template's outer 30 bp match below 80% identity (`--min-identity`), usually because
  a large indel separates the true end from the rest of the LTR, those 30 bp
  (`--anchor-len`) are searched for up to 5,000 bp (`--max-indel`) beyond the call:
  an *anchor*. An anchor is refused where a complete copy of the template LTR sits
  (a solo LTR or a neighboring copy, not this element's end). `--no-anchor` turns
  anchors off.
- **Vote.** Per side, placements whose outer 30 bp are ≥ 80% identical and whose
  whole LTR is ≥ 60% identical vote. The median wins when at least 2 of them
  (`--min-support`) agree with it within 2 bp. The end then moves outward by at
  least 1 bp (`--min-ext`), or inward by at most 10 bp (`--max-trim`; `0` never
  trims). A move that would add `N` bases (an assembly gap) is refused. The new ends
  come only from the templates: nothing searches for a TSD or TG..CA (the `motif`
  column above records this as `tg...ca`).
- **Split calls.** A detector sometimes calls one element twice, staggered, each call
  holding one of its true ends. If a move reaches, within 5 bp, the outer end of
  another call that lies inside the new span, has the same strand, depth and
  nesting, hosts nothing and has no TSD, the two are one element. The other call is
  removed (its row, FASTA record and `nest_status` references), provided the
  surviving call passes Kmer2LTR. Family labels are not compared, because a split
  call's family was clustered from its partly wrong LTR pair.
- **Conflicts.** A move is rejected if it would take the element outside its host,
  give it its host's exact span or the exact span of a call inside it, make it
  contain or overlap another element, or leave a nested element outside a trimmed
  call. When two moves claim the same bases, or the same ends (two calls of one
  element, one inside the other, one extended and one trimmed to the same span),
  the first in genome order wins.
- **Re-measurement.** Kmer2LTR re-measures each moved element on the LTR pair it was
  called with, taking the templates' outer ends as given. Each LTR's inner end moves
  only as far as real homology with the other LTR reaches, and bases with no partner
  count as gaps rather than substitutions, so `k2p` and `k2p_time` are not inflated.
  The move is kept if Kmer2LTR's status is `pass`, both LTRs are ≥ 100 bp, the
  alignment is ≥ 90 bp and the element passes detection's size floor
  (end − start ≥ 300). It is refused if the re-measured element would break the
  overlap rule (*Overlapping calls*, section 4.1) with a call it kept to before the
  move: an LTR overlapping one of a call nested in it or hosting it, or an end shared
  with one. Only the re-measurement shows this, because the conflicts above are
  checked on spans. A pair already in conflict before the move (a tenant left in its
  host's LTR, section 4.3) is not counted against the move, and when two moves conflict
  only with each other, the first in genome order wins. Kmer2LTR then supplies the
  element's columns 1–29 except `orientation`, which stays as stored; `nest_status`
  follows its neighbors' new names.

The stage runs once, by design. On rice, a second pass over its own output changed
112 more elements, and a curated LTR library agreed with only 88% of those ends,
against 99.5% for the first pass. `ltrquest` uses the defaults; to try
other settings, run `ltrquest-reboundary --posthoc` (see `--help`) on a run made
with `--no-reboundary --keep-intermediates`.

**What changes:** the `_clean_` depth tables and FASTAs (a host's masked region is
repainted over added bases and restored over trimmed ones), the GFF3s, the plots, and
the new `<prefix>_reboundary.tsv` (kept only with `--keep-intermediates`). The unfiltered `<prefix>_depth<N>_ltr.{tsv,fa}` tables
(also kept only with it)
keep the calls as detected. Family labels do not change; the pooled `<run>_all_ltr.*`
clustering files describe the calls before re-boundarying. A genome with no
`_depth<N>_clean_ltr.tsv` tables is skipped with a warning.

**Requirements:** BLAST+ (`blastn` from `PATH` or `--blastn`, with `makeblastdb`
beside it or on `PATH`) and a Kmer2LTR with the homology-paired credit. An older
Kmer2LTR is refused before anything changes: `ltrquest-reboundary` exits 3 and
`ltrquest` stops.

> When run through the Nextflow pipeline (on by default there too; `--reboundary false`
> skips it), the task that writes the re-bounded `_clean_` depth tables and
> `<prefix>_reboundary.tsv` is pooled, because an element's templates are drawn
> from every sample, but each file is still published to its own sample's
> `${params.outdir}/<sample>/`. See [nextflow.md](nextflow.md#what-the-pipeline-does).

### 8.1 `<prefix>_reboundary.tsv`

One row per candidate (a call whose templates would move an end under the rules
above, moved or not), plus one row per removed split call. Every genome it ran on gets
one; it holds only the header when there is nothing to report.

| Column | Meaning |
|---|---|
| `old_seq_id`, `new_seq_id` | the call before and after; `new_seq_id` is `.` when rejected, and on a `merged` row the surviving call's new id |
| `family` | the call's family |
| `method`, `templates` | always `templates`; the templates placed on the call, `prefix:seq_id`, comma-separated |
| `decision`, `reason` | `moved`; `merged`, with reason `split_call`; or `rejected`, with a reason from the table below |
| `ext5`, `ext3` | signed bp at the element's biological 5′ / 3′ end: positive added, negative trimmed. On a rejected row, what the placements proposed |
| `end_source5`, `end_source3` | `core` (template aligned end to end), `anchor` (template's outer bases found past a large indel), `.` (end did not move) |
| `id_outer5`, `id_outer3` | identity of the templates' outer 30 bp at that end; `.` where it did not move |
| `support5`, `support3` | template placements that agree on that end |
| `k2l_status` | Kmer2LTR's status for the new record |
| `tsd_called`, `tsd_new`, `tsd_null` | TSD at the called ends (always `.` or `NA`: only calls without one are moved), at the new ends, and at the new ends with the right flank displaced 1 kb (the null). The last two are `.` where no end moved |
| `k2p_called`, `k2p_new` | LTR-pair K2P before and after |
| `unpaired5`, `unpaired3` | bases at the new pair's biological 5′ / 3′ outer end with no partner in the other LTR: gap columns in the new alignment |
| `obstacle5`, `obstacle3` | what sits between the two LTRs at an old end that moved out: `gap:<bp>`, `mm:<rate>`, `none`, `.` |
| `merged_into` | on a `merged` row, the survivor's new id |

5′ and 3′ follow the element's `strand`; for an element with none, the strand its
templates imply (a majority vote), else genomic left and right.

| `reason` | Why the call was not moved |
|---|---|
| `gate_identity` | no placement passed the identity gates |
| `no_support` | fewer than `--min-support` placements agree |
| `gap_in_added` | the move would add `N` bases |
| `host_exceeded`, `duplicates_host`, `duplicates_element`, `engulfs_element`, `overlaps_element`, `nest_broken` | a conflict, as above |
| `duplicates_move` | another call, earlier in genome order, moves to the same span; this one keeps its call |
| `merge_partner_failed` | this call was to be merged into another whose move then failed; both keep their calls |
| `record_missing` | the call has no FASTA record |
| `record_mismatch` | the stored FASTA record is not the genome at the row's span |
| `kmer2ltr_not_pass` | Kmer2LTR's status for the new record is not `pass` |
| `kmer2ltr_moved_ends` | Kmer2LTR moved an end it was given (an internal consistency check) |
| `length_filter` | after re-measuring: an LTR under 100 bp, the alignment under 90 bp, or the element under detection's size floor (end − start < 300) |
| `ltr_conflict` | after re-measuring: an LTR would overlap one of a call nested in this one or hosting it, or the two would share an end (with a call as it stands, or with a call moved earlier in genome order); this one keeps its call |

### 8.2 GFF3

Moved elements carry `boundary_source=templates` and `boundary_shift=<ext5>,<ext3>`
(signed bp at the biological 5′ / 3′ end; negative = trimmed). The GFF3 writer reads
the sidecar (`--reboundary-map`) so family and strand-provenance attributes still find
renamed elements. Sidecars from the pre-release design (decision `extended`) still
load.

### 8.3 Finished runs

`--posthoc` needs the intermediates, so the run must have used `--keep-intermediates`.

```bash
ltrquest-reboundary --posthoc --indir RUN --prefix P1 P2 ... --genome G1 G2 ... --threads 32
ltrquest-reboundary --restore --indir RUN --prefix P1 P2 ...     # undo
```

`--posthoc` moves the originals (clean tables and FASTAs, GFF3s, plots) into
`<prefix>_pre_reboundary/` once, re-bounds, then rewrites the GFF3s and plots
(`--no-plots` skips the plots). A second `--posthoc` run starts again from that
backup, so trying other settings never stacks on an earlier result. It refuses a run
that `ltrquest` already re-boundaried in place (a sidecar and no backup):
that would be a second pass. Benchmarks: `benchmarks/reboundary/README.md`.

## 9. Plots (`<prefix>_plots/`)

| Output | What it shows |
|---|---|
| `struct/<clade>_average.pdf` | Average element structure per classification clade with bootstrap 95% CI whiskers on feature positions |
| `struct/<clade>_individual.pdf` | Every element of the clade drawn individually, all domains shown |
| `struct/all_elements.pdf` | All elements on one page |
| `<prefix>_summary.pdf` | Multi-page summary |
| `<prefix>_TEGV.html` | Self-contained interactive genome browser (open in web browser) |

Re-plotting a finished run with `ltrquest-plots` needs `--keep-intermediates`: the summary PDF's per-family pages and TEGV's gene track read files a default run removes.

## 10. Benchmarks

How accurate LTRquest is, how long it takes, and how it compares with other
callers: [benchmarks/README.md](../benchmarks/README.md).

## 11. Settings not in the short help

`ltrquest --help` lists the few flags most runs need; `--help-all` lists the rest.
Detail on the ones that need explaining:

- **`--mutation-rate R`** (default `3e-8`, a plant-generic rate in substitutions per
  site per year) only feeds `k2p_time`, the insertion age in years. Nothing else
  depends on it. Against a mammalian genome the run still finishes and still
  reports ages, computed from the wrong clock, with nothing in the output to
  flag it. Set a rate that describes your organism; `7e-9` is the usual choice
  for *Arabidopsis thaliana*.
- **`--max-rounds N`** (default up to 10) and **`--terminate_count N`** (default
  100). Rounds continue until one finds fewer than N new elements. `--max-rounds 1`
  finds un-nested elements only; 2 adds one level of nesting.
- **`--fp-mask-threshold F`** (default 0.10). If more than this fraction of the
  elements are in false-positive families, those are masked and everything is
  re-run on the masked genome.
- **`--keep-weak-hmm-pass2-matches`** and **`--pass2-aligner minimap2|blast`**
  (default `minimap2`). The second TEsorter2 pass assigns a clade to elements the
  first pass left unclassified; by default a pass-2 match whose target has protein
  hits but no clade is dropped.
- **`--detect-args "KEY=VALUE ..."`** passes options to `ltrquest.detect` in every
  round, for example `--detect-args "tesorter-rule=70-70-80"`;
  **`--detect-args-from-round N "..."`** does the same from round N on. Detection
  settings are recorded per genome, so changing them needs `--redetect` on an
  existing run folder.
- **`--out_prefix NAME`**: one genome, the prefix of every output (default
  `<genome>_LTRs`), which also names the families; several genomes, the name of
  the shared pool (default `merged`). See [section 5](#5-multiple-genomes-shared-family-names).
- **`--redetect`** detects every genome again, ignoring records from an earlier run.
- **`--no-trf`** skips TRF-mod tandem-repeat masking; **`--no-plots`** skips the
  plots.
