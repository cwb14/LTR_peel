# Benchmarks

Each directory holds one benchmark: a page with its results and how they were made, and the code to rerun it.

**[comparison/](comparison/README.md)** compares LTR_peel 1.0.1 (commit 5042a1e) with eight other LTR retrotransposon
callers (LTR_retriever, LTRharvest, LTR_Finder, Look4LTRs, HiRepEat, HiTE, LTRpred and DANTE_LTR) on simulated genomes
with exact truth, five real genomes scored against curated repeat libraries, two allopolyploids and real insertion
polymorphisms genotyped from long reads. LTR_peel found the most nested and protein-less elements and, with
DANTE_LTR, the most old ones, while libraries built from the calls of HiTE, LTR_retriever and DANTE_LTR reproduced the
curated annotation best. The page lists every command, dataset, unfinished run and conflict of interest (LTR_peel's
author ran it).

**[fastltr/](fastltr/README.md)** asks whether LTR_peel's FastLTR caller recovers the elements that `--ltr-timeout`
throws away when it stops LTRharvest and LTR_FINDER on slow chunks such as centromeric and rDNA arrays, and what that
costs elsewhere. In chunks where a 30 s cap stopped both callers, the cap alone kept 0 of 97 (*A. thaliana*) and 1 of
76 (rice) curated-supported elements, and the cap with `--fastltr rescue` kept 66 of 97 and 65 of 76. It also scores
elements planted into real arrays, a simulated genome and the human T2T assembly.

**[reboundary/](reboundary/README.md)** measures how well `--reboundary` moves the ends of truncated LTR-RT calls to
their true positions, using two signals the stage never uses: target-site duplications and a curated library. On rice
it changed 1,726 elements, and at 99.5% of the 1,965 new ends that the curated library could judge, the library's
terminus sat exactly there.
