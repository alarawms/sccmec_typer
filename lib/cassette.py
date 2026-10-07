"""Putative SCC cassette boundaries from orfX direct repeats (assembly mode only).

SCC elements integrate at the 3' end of orfX (rlmH); each inserted element is followed by
a direct repeat (DR) of the orfX 3'-end core. The core is taken from the genome itself
(last CORE_LEN bp of its orfX), and copies are searched downstream of orfX on the same
contig. DR copies with <= HIGH_CONF_MM mismatches split the region into segments; every
mec/ccr locus is assigned to a segment, which shows e.g. whether an extra ccr lies in the
SCCmec or in an adjacent SCC element.

This is reported, not used for typing: proper attL/attR calls also need the inverted
repeats and the empty-site junction, and DR copies can diverge (candidates up to
CANDIDATE_MM mismatches are listed with their mismatch counts).
"""
import os
import tempfile

from lib.aligner import run_minimap2

ORFX_REF = os.path.join(os.path.dirname(__file__), "data", "orfX_NCTC8325.fasta")
CORE_LEN = 18
SEARCH_BP = 120_000
HIGH_CONF_MM = 2
CANDIDATE_MM = 4
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _revcomp(s):
    return s.translate(COMPLEMENT)[::-1]


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line.upper())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def locate_orfx(assembly, threads=4):
    """(contig, end3, strand) of orfX's 3' end in the assembly, or None."""
    with tempfile.TemporaryDirectory() as tmp:
        paf = os.path.join(tmp, "orfx.paf")
        run_minimap2(assembly, ORFX_REF, paf, threads, preset="asm20", target_is_db=False)
        best = None
        with open(paf) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) < 12:
                    continue
                qlen, qs, qe = int(p[1]), int(p[2]), int(p[3])
                strand, contig, ts, te, matches, blen = p[4], p[5], int(p[7]), int(p[8]), int(p[9]), int(p[10])
                score = matches / max(1, blen) * (qe - qs) / qlen
                if best is None or score > best[0]:
                    best = (score, qlen, qs, qe, strand, contig, ts, te)
    if not best or best[0] < 0.6:
        return None
    _, qlen, qs, qe, strand, contig, ts, te = best
    end3 = te + (qlen - qe) if strand == "+" else ts - (qlen - qe)
    return contig, end3, strand


def dr_candidates(seq, end3, strand, core_len=CORE_LEN, max_mm=CANDIDATE_MM, search=SEARCH_BP):
    """[(offset_from_orfX_3prime_end, mismatches)] of the orfX core downstream of orfX.

    Offsets are in orfX orientation; offset 0 is the core at the end of orfX itself (attL).
    """
    if strand == "+":
        core = seq[end3 - core_len:end3]
        window = seq[end3 - core_len:end3 + search]
    else:
        core = _revcomp(seq[end3:end3 + core_len])
        window = _revcomp(seq[max(0, end3 - search):end3 + core_len])
    if len(core) < core_len:
        return []
    out = []
    for i in range(len(window) - core_len + 1):
        mm = sum(1 for a, b in zip(window[i:i + core_len], core) if a != b)
        if mm <= max_mm:
            out.append((i, mm))
    return out


def segment_loci(candidates, loci, high_conf_mm=HIGH_CONF_MM):
    """Split by high-confidence DRs; assign (gene, offset) loci to segments.

    Returns [{"start": kb, "end": kb or None (open), "genes": [...]}].
    """
    bounds = sorted({0} | {off for off, mm in candidates if mm <= high_conf_mm})
    segs = [{"start": b, "end": (bounds[i + 1] if i + 1 < len(bounds) else None), "genes": []}
            for i, b in enumerate(bounds)]
    for gene, off in sorted(loci, key=lambda x: x[1]):
        if off < 0:
            continue
        for s in segs:
            if off >= s["start"] and (s["end"] is None or off < s["end"]):
                s["genes"].append(gene)
                break
    return segs


def cassette_map(assembly, loci_hits, threads=4):
    """Cassette map dict for the JSON/TSV, or None when orfX is not found.

    loci_hits: de-duplicated mec/ccr locus hits (gene, contig, start, end).
    """
    orfx = locate_orfx(assembly, threads)
    if not orfx:
        return None
    contig, end3, strand = orfx
    seq = read_fasta(assembly).get(contig, "")
    cands = dr_candidates(seq, end3, strand)
    loci = []
    for h in loci_hits:
        if h.get("contig") != contig or not (h["gene"].startswith("ccr") or h["gene"].startswith("mec")):
            continue
        mid = (h["start"] + h["end"]) // 2
        off = mid - end3 if strand == "+" else end3 - mid
        if 0 <= off <= SEARCH_BP:
            loci.append((h["gene"], off))
    segs = segment_loci(cands, loci)
    return {
        "orfx_contig": contig, "orfx_3prime_end": end3, "orfx_strand": strand,
        "dr_candidates": [{"offset": o, "mismatches": m} for o, m in cands],
        "segments": segs,
        "summary": summarize(segs),
        "note": "putative DR boundaries from the genome's own orfX 3'-end core; not used for typing",
    }


def summarize(segs):
    parts = []
    for i, s in enumerate(segs, 1):
        if not s["genes"]:
            continue
        end = f"{s['end'] / 1000:.1f}" if s["end"] is not None else "?"
        parts.append(f"seg{i}({s['start'] / 1000:.1f}-{end}kb):{','.join(s['genes'])}")
    return "|".join(parts)


def interpret(segs):
    """Short note when mec and an extra ccr sit in different segments (tandem SCCs)."""
    with_mec = [i for i, s in enumerate(segs) if any(g.startswith("mec") for g in s["genes"])]
    ccr_only = [i for i, s in enumerate(segs)
                if s["genes"] and all(g.startswith("ccr") for g in s["genes"])]
    if with_mec and ccr_only:
        return ("ccr-only segment(s) adjacent to the SCCmec segment: extra ccr likely in a separate "
                "(tandem) SCC element, not inside SCCmec")
    return ""
