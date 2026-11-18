"""True element ends from a segment pair (tool-independent).

r_seq carries flank + TSD + element + TSD + flank; o_seq carries flank + one TSD + flank.
Global alignment (not semi-global): every r_seq base must be explained, so the element's long
r-only gap is forced and placed where it maximizes the flank alignment; a semi-global aligner
with free end gaps would simply skip the flank rather than pay for a long gap.
Align them, find the long r_seq-only gap of the element, left/right-normalize it, and return the
innermost element consistent with the alignment plus the duplicated sequence.
Which gap is the element's: without a hint, the only r_seq-only gap >= MIN_GAP in the segment;
with hint = the candidate locus span ([start, end) offsets in r_seq), the only one overlapping the
hint by >= 1 bp. Other long gaps (unrelated indels in the flanks) are then ignored unless the near
rule catches them: any indel >= NEAR_INDEL within NEAR_BP of the normalized gap still means None.
The alignment is run on both orientations (forward and reversed strings), because the aligner
places a gap leftmost and a SNP at the junction can then shift the ends; the valid result with
the larger duplicated length is kept (tie: forward). The reversed alignment's indels are mapped
back to forward offsets before the checks, so the hint is applied to both orientations in the
same frame (equivalent to mirroring it to (n - end, n - start) for the reversed strings).
"""
from dataclasses import dataclass

import parasail

from lg_config import MAX_DUP, MAX_LOCUS, MIN_GAP, NEAR_BP, NEAR_INDEL

# Affine gap open 26 / extend 1 (minimap2 asm20 long-gap costs).
_MATRIX = parasail.matrix_create("ACGTN", 1, -4)
GAP_OPEN, GAP_EXTEND = 26, 1


@dataclass
class TrueEnds:
    start: int      # 0-based offset in r_seq; element is r_seq[start:end]
    end: int
    dup_len: int
    dup_seq: str


def _indels(o_seq, r_seq):
    """Internal indels as (kind, r_start, r_end, length) in r_seq coordinates.

    CIGAR convention (o_seq = s1/query, r_seq = s2/reference): 'D' = bases only in r_seq,
    'I' = bases only in o_seq. Leading/trailing I or D ops are end gaps (o_seq padding or
    r_seq overhang) and are dropped; leading D still advances the r_seq position.
    """
    res = parasail.nw_trace_scan_32(o_seq, r_seq, GAP_OPEN, GAP_EXTEND, _MATRIX)
    cg = res.cigar
    ops = [(parasail.Cigar.decode_op(x).decode(), parasail.Cigar.decode_len(x)) for x in cg.seq]
    rp = cg.beg_ref
    while ops and ops[0][0] in "ID":        # leading end gap: advance r_seq position past it
        op, n = ops.pop(0)
        rp += n if op == "D" else 0
    while ops and ops[-1][0] in "ID":
        ops.pop()
    out = []
    for op, n in ops:
        if op in "=XM":
            rp += n
        elif op == "D":
            out.append(("D", rp, rp + n, n))
            rp += n
        elif op == "I":
            out.append(("I", rp, rp, n))
        else:
            raise ValueError(f"unexpected CIGAR op {op!r}")
    return out


def _normalize(r: str, g0: int, G: int):
    """Return (a_L, a_R): left- and right-normalized start of the G-long gap starting at g0."""
    n = len(r)
    a_L = g0
    while a_L > 0 and r[a_L - 1] == r[a_L - 1 + G]:
        a_L -= 1
    a_R = g0
    while a_R + G < n and r[a_R] == r[a_R + G]:
        a_R += 1
    return a_L, a_R


def _one(r: str, ind, hint=None):
    """Clean-gap checks and normalization on forward-coordinate indels; TrueEnds or None.

    The element gap is the only r-only gap >= MIN_GAP (hint None), or the only one overlapping
    hint = [start, end) by >= 1 bp; the other long gaps stay in `ind` for the near rule."""
    big = [x for x in ind if x[0] == "D" and x[3] >= MIN_GAP]
    if hint is not None:
        big = [x for x in big if max(x[1], hint[0]) < min(x[2], hint[1])]
    if len(big) != 1:
        return None
    _, g0, _, G = big[0]
    a_L, a_R = _normalize(r, g0, G)
    dup_len = a_R - a_L
    if not 0 <= dup_len <= MAX_DUP:
        return None
    for x in ind:
        if x is big[0] or x[3] < NEAR_INDEL:
            continue
        for pt in (a_L, a_R + G):
            if max(x[1] - pt, pt - x[2], 0) <= NEAR_BP:
                return None
    return TrueEnds(start=a_R, end=a_L + G, dup_len=dup_len, dup_seq=r[a_L:a_R])


def true_ends(r_seq: str, o_seq: str, hint: "tuple[int, int] | None" = None) -> "TrueEnds | None":
    """Return TrueEnds for the clean element gap, or None if not clean.

    hint: the candidate locus span as [start, end) offsets in r_seq. None: the element gap must be
    the only r-only gap >= MIN_GAP in the segment. Given: it must be the only one overlapping the
    hint (none or two or more -> None); other long gaps count only through the near rule."""
    if not r_seq or not o_seq:
        raise ValueError("true_ends: r_seq and o_seq must be non-empty")
    if hint is not None and not hint[0] < hint[1]:
        raise ValueError(f"true_ends: hint must be (start, end) with start < end, got {hint!r}")
    r, o = r_seq.upper(), o_seq.upper()
    if len(o) > MAX_LOCUS:
        return None
    n = len(r)
    fwd = _one(r, _indels(o, r), hint)
    # reversed strings: indels mapped back to forward offsets, so the same (forward) hint applies
    rev = _one(r, [(k, n - e, n - b, ln) for k, b, e, ln in _indels(o[::-1], r[::-1])], hint)
    if fwd is None or (rev is not None and rev.dup_len > fwd.dup_len):
        return rev
    return fwd
