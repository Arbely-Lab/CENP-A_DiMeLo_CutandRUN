#USAGE

# Auto mode (uses pileup since regions provided)
#python mod_extract.py sample.bam --mod C+m --bed regions.bed \
#  --threshold 0.5 --region-fraction-mode [weighted,unweighted,both] --perpos-output per_pos.csv --region-output per_region.csv

# Single position but still get a region summary for that one site
#python mod_extract.py sample.bam --mod C+m --positions chr1:123456 \
#  --perpos-output per_pos.csv --region-output per_region.csv

# Force read mode + per-read table (streamed)
#python mod_extract.py sample.bam --mod A+m --positions chr1:1000-1020 chrX:2001 \
#  --mode read --perpos-output per_pos.csv --per-read-output per_read.csv


#srun -t 12:00:00 --cluster smp --partition smp --mem 50G --pty bash
#cd /ix1/yarbely/mam835/Dimelo
#module load python/3.12.0
#source activate Dimelo
#python mod_extract_mem_v2.py ./Ragini/U2OS_2/U2OS_2_CENPA_mods_clean.bam --mod C+m --positions NC_060925.1:1-5000 --threshold 0.5 --region-fraction-mode both --perpos-output U2OS_per_pos_chr1_tv3.csv --region-output U2OS_per_region_chr1_tv3.csv


#!/usr/bin/env python3
# Memory-efficient mod extractor with pileup (position-first) and read-iteration (read-first) paths.
# Fraction = (# reads with prob >= threshold at pos) / (# reads covering pos with canonical base).

import pysam
import argparse
import re
import csv
import os
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Tuple, Iterable, Optional, Set

# ---------------------- position helpers ----------------------

def parse_positions_cli(position_args: Optional[List[str]]) -> List[Tuple[str, int, int]]:
    """
    Accept any chromosome name up to the first colon (so dots, dashes, etc. are fine).
    Examples: 'NC_060925.1:1-5000', 'chr1:1000', 'chrUn_KI270442v1:50-200'
    """
    out = []
    if not position_args:
        return out
    for token in position_args:
        m = re.match(r"^([^:]+):(\d+)-(\d+)$", token)
        if m:
            chrom, s, e = m.groups()
            out.append((chrom, int(s), int(e)))
            continue
        m = re.match(r"^([^:]+):(\d+)$", token)
        if m:
            chrom, p = m.groups()
            p = int(p)
            out.append((chrom, p, p))
            continue
        raise ValueError(f"Could not parse position/range: '{token}'. Expected 'chr:pos' or 'chr:start-end'.")
    return out

def parse_bed(path: str) -> List[Tuple[str, int, int]]:
    # BED is 0-based half-open; convert to 1-based inclusive
    out = []
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            chrom = cols[0]
            start0 = int(cols[1]); end0 = int(cols[2])
            if end0 <= start0:
                continue
            out.append((chrom, start0 + 1, end0))
    return out

def expand_positions(ranges: Iterable[Tuple[str, int, int]]) -> List[Tuple[str, int]]:
    s: Set[Tuple[str, int]] = set()
    for chrom, s1, e1 in ranges:
        for p in range(s1, e1 + 1):
            s.add((chrom, p))
    return sorted(s)

def merge_intervals(ranges):
    """Merge overlapping/adjacent 1-based inclusive intervals per chromosome."""
    by_chr = defaultdict(list)
    for chrom, s, e in ranges:
        if s > e: s, e = e, s
        by_chr[chrom].append((s, e))
    merged = {}
    for chrom, segs in by_chr.items():
        segs.sort()
        out = []
        cs, ce = segs[0]
        for s, e in segs[1:]:
            if s <= ce + 1:
                ce = max(ce, e)
            else:
                out.append((cs, ce))
                cs, ce = s, e
        out.append((cs, ce))
        merged[chrom] = out
    return merged  # {chrom: [(start,end), ...]}

# ---------------------- mod parsing helpers ----------------------

def parse_mod_string(mod_str: str) -> Tuple[str, str]:
    m = re.match(r"^([ACGTN])\+([A-Za-z])$", mod_str)
    if not m:
        raise ValueError(f"Invalid --mod '{mod_str}'. Expected like 'C+m' or 'A+m'.")
    return m.group(1), m.group(2)

def build_mappings(read: pysam.AlignedSegment):
    q2r: Dict[int, int] = {}
    r2q: Dict[int, int] = {}
    for qpos, rpos in read.get_aligned_pairs(matches_only=True):
        if qpos is not None and rpos is not None:
            ref1 = rpos + 1
            q2r[qpos] = ref1
            r2q[ref1] = qpos
    return q2r, r2q

def get_mod_sites_from_read(read: pysam.AlignedSegment, base: str, mod_code: str) -> List[Tuple[int, float]]:
    if not hasattr(read, "modified_bases") or read.modified_bases is None:
        return []
    q2r, _ = build_mappings(read)
    out: List[Tuple[int, float]] = []
    for (canon, strand, mod), entries in read.modified_bases.items():
        if canon != base or mod != mod_code:
            continue
        for entry in entries:
            if isinstance(entry, tuple) and len(entry) >= 2:
                qpos, prob_byte = entry[0], entry[1]
            else:
                continue
            if qpos in q2r and prob_byte is not None:
                prob = float(prob_byte) / 255.0
                ref_pos = q2r[qpos]
                out.append((ref_pos, prob))
    return out

def get_base_coverage_positions_for_read(read: pysam.AlignedSegment, base: str) -> Set[int]:
    if read.query_sequence is None:
        return set()
    _, r2q = build_mappings(read)
    cov: Set[int] = set()
    for ref1, qpos in r2q.items():
        if 0 <= qpos < len(read.query_sequence):
            b = read.query_sequence[qpos]
            if b and b.upper() == base:
                cov.add(ref1)
    return cov

def get_prob_at_qpos_for_mod(read, base: str, mod_code: str, qpos: int) -> Optional[float]:
    """Return normalized prob (0..1) for chosen (base, mod_code) at this query position, else None."""
    if not hasattr(read, "modified_bases") or read.modified_bases is None:
        return None
    for (canon, strand, mod), entries in read.modified_bases.items():
        if canon != base or mod != mod_code:
            continue
        for entry in entries:
            if isinstance(entry, tuple) and len(entry) >= 2 and entry[0] == qpos:
                prob_byte = entry[1]
                if prob_byte is not None:
                    return float(prob_byte) / 255.0
    return None

# ---------------------- one-pass aggregators ----------------------

@dataclass
class PosAgg:
    covered: int = 0        # reads with canonical base covering pos
    modified: int = 0       # reads with prob >= threshold
    n_probs: int = 0        # reads that reported a probability
    mean: float = 0.0       # Welford running mean of probs
    M2: float = 0.0         # Welford running sum of squares

    def add_prob(self, x: float):
        self.n_probs += 1
        delta = x - self.mean
        self.mean += delta / self.n_probs
        delta2 = x - self.mean
        self.M2 += delta2 * delta2

    def variance(self) -> float:
        return (self.M2 / (self.n_probs - 1)) if self.n_probs > 1 else 0.0

# ---- read-iteration path (can stream per-read CSV) ----

def parse_bam_streaming(
    bam_path: str,
    ranges: List[Tuple[str, int, int]],
    mod_str: str,
    threshold: float,
    per_read_writer: Optional[csv.writer],
    per_read_positions: Optional[List[Tuple[str, int]]],
    progress_every: int = 100000
) -> Dict[Tuple[str, int], PosAgg]:
    import sys, time

    base, mod_code = parse_mod_string(mod_str)
    agg: Dict[Tuple[str, int], PosAgg] = defaultdict(PosAgg)

    merged = merge_intervals(ranges) if ranges else None
    pos_index: Dict[Tuple[str, int], int] = {}
    if per_read_writer and per_read_positions:
        pos_index = {k: i for i, k in enumerate(per_read_positions)}

    start = time.time()
    processed = 0

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        def process_read(read: pysam.AlignedSegment):
            nonlocal processed
            if read.is_unmapped:
                return
            processed += 1
            chrom = read.reference_name

            mod_sites = get_mod_sites_from_read(read, base, mod_code)
            cov_sites = get_base_coverage_positions_for_read(read, base)

            # coverage increments
            for ref_pos in cov_sites:
                if merged and not any(s <= ref_pos <= e for s, e in merged.get(chrom, [])):
                    continue
                agg[(chrom, ref_pos)].covered += 1

            # probabilities / modified increments
            for ref_pos, prob in mod_sites:
                if merged and not any(s <= ref_pos <= e for s, e in merged.get(chrom, [])):
                    continue
                pa = agg[(chrom, ref_pos)]
                pa.add_prob(prob)
                if prob >= threshold:
                    pa.modified += 1

            # per-read wide row (streamed)
            if per_read_writer and per_read_positions:
                row = [""] * (1 + len(per_read_positions))
                row[0] = read.query_name
                if pos_index:
                    for ref_pos, prob in mod_sites:
                        k = (chrom, ref_pos)
                        idx = pos_index.get(k)
                        if idx is not None:
                            row[1 + idx] = f"{prob:.3f}"
                per_read_writer.writerow(row)

            if processed % progress_every == 0:
                dt = time.time() - start
                rate = processed / dt if dt > 0 else 0.0
                print(f"[progress] processed {processed:,} reads ({rate:,.0f} reads/s)",
                      file=sys.stderr, flush=True)

        if merged:
            # region-restricted iteration (fast); fetch uses 0-based, half-open
            for chrom, intervals in merged.items():
                for s1, e1 in intervals:
                    for read in bam.fetch(chrom, s1 - 1, e1):
                        process_read(read)
        else:
            for read in bam:
                process_read(read)

    return agg

# ---- pileup path (position-first; lowest memory) ----

def parse_bam_streaming_pileup(
    bam_path: str,
    ranges: List[Tuple[str, int, int]],
    mod_str: str,
    threshold: float,
    progress_every_cols: int = 10000,
) -> Dict[Tuple[str, int], PosAgg]:
    import sys, time

    if not ranges:
        raise ValueError("Pileup mode requires --positions and/or --bed.")

    base, mod_code = parse_mod_string(mod_str)
    agg: Dict[Tuple[str, int], PosAgg] = defaultdict(PosAgg)
    merged = merge_intervals(ranges)

    start = time.time()
    cols_seen = 0

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for chrom, intervals in merged.items():
            for s1, e1 in intervals:
                for col in bam.pileup(
                    chrom, s1 - 1, e1,
                    truncate=True,
                    stepper="all",
                    ignore_overlaps=False,
                    ignore_orphans=True,
                ):
                    pos1 = col.reference_pos + 1
                    key = (chrom, pos1)

                    for pr in col.pileups:
                        if pr.is_del or pr.is_refskip:
                            continue
                        read = pr.alignment
                        qpos = pr.query_position
                        if qpos is None or read.query_sequence is None:
                            continue
                        b = read.query_sequence[qpos].upper()
                        if b != base:
                            continue

                        pa = agg[key]
                        pa.covered += 1

                        prob = get_prob_at_qpos_for_mod(read, base, mod_code, qpos)
                        if prob is not None:
                            pa.add_prob(prob)
                            if prob >= threshold:
                                pa.modified += 1

                    cols_seen += 1
                    if cols_seen % progress_every_cols == 0:
                        dt = time.time() - start
                        rate = cols_seen / dt if dt > 0 else 0.0
                        print(f"[progress] processed {cols_seen:,} columns "
                              f"({rate:,.0f} cols/s) in {chrom}:{s1}-{e1}",
                              file=sys.stderr, flush=True)

    return agg

# ---------------------- region summarization ----------------------

@dataclass(frozen=True)
class Region:
    chrom: str
    start: int  # 1-based inclusive
    end: int    # 1-based inclusive

def to_regions(ranges):
    return [Region(c, s, e) for c, s, e in ranges]

def summarize_regions(
    regions: List[Region],
    agg: Dict[Tuple[str, int], PosAgg],
    mod_str: str,
    threshold: float,
    mode: str = "weighted",  # "weighted" | "unweighted" | "both"
):
    """
    Builds region rows. Always returns Covered_reads_sum, Covered_reads_avg, Modified_reads_sum, Num_prob_observations,
    plus region fractions per 'mode':
      - weighted: Modified_reads_sum / Covered_reads_sum
      - unweighted: mean over sites of (Modified_reads / Covered_reads), excluding sites with Covered_reads == 0
    """
    rows = []
    mod_label = mod_str.replace("+", "plus")

    by_chr: Dict[str, Set[int]] = defaultdict(set)
    for (chrom, pos) in agg.keys():
        by_chr[chrom].add(pos)

    for r in regions:
        cov_sum = 0
        mod_sum = 0
        prob_weighted_sum = 0.0
        prob_count = 0

        site_fracs: List[float] = []

        for pos in by_chr.get(r.chrom, set()):
            if r.start <= pos <= r.end:
                pa = agg.get((r.chrom, pos))
                if not pa:
                    continue
                cov_sum += pa.covered
                mod_sum += pa.modified
                if pa.n_probs > 0:
                    prob_weighted_sum += pa.mean * pa.n_probs
                    prob_count += pa.n_probs
                if pa.covered > 0:
                    site_fracs.append(pa.modified / pa.covered)

        region_fraction_weighted = (mod_sum / cov_sum) if cov_sum > 0 else 0.0
        region_fraction_unweighted = (sum(site_fracs) / len(site_fracs)) if site_fracs else 0.0
        region_avg_prob = (prob_weighted_sum / prob_count) if prob_count > 0 else 0.0

        region_len = (r.end - r.start + 1) if (r.end >= r.start) else 0
        covered_avg = (cov_sum / region_len) if region_len > 0 else 0.0

        base = {
            "Region": f"{r.chrom}:{r.start}-{r.end}",
            f"Region_AverageProb_{mod_label}": f"{region_avg_prob:.3f}",
            "Covered_reads_sum": cov_sum,
            "Covered_reads_avg": f"{covered_avg:.3f}",
            "Modified_reads_sum": mod_sum,
            "Num_prob_observations": prob_count,
        }
        if mode in ("weighted", "both"):
            base[f"Region_Fraction_{mod_label}(thresh={threshold})"] = f"{region_fraction_weighted:.3f}"
        if mode in ("unweighted", "both"):
            base[f"Region_Fraction_unweighted_{mod_label}(thresh={threshold})"] = f"{region_fraction_unweighted:.3f}"

        rows.append(base)

    return rows

# ---------------------- sanity checks ----------------------

def ensure_bam_index_exists(bam_path: str):
    candidates = [bam_path + ".bai", os.path.splitext(bam_path)[0] + ".bai"]
    if not any(os.path.exists(p) for p in candidates):
        raise FileNotFoundError(
            f"No BAM index (.bai) found for '{bam_path}'. Create one with:\n  samtools index {bam_path}"
        )

# ---------------------- CLI / main ----------------------

def main():
    parser = argparse.ArgumentParser(
        description="Extract a chosen base modification (e.g., C+m, A+m) from a BAM with low memory.\n"
                    "Fraction = (# reads with prob >= threshold) / (# reads covering pos with canonical base)."
    )
    parser.add_argument("bam", help="Input BAM (indexed).")
    parser.add_argument("--mod", required=True, help="Modification to extract (e.g., C+m for 5mC, A+m for m6A).")
    parser.add_argument("--positions", nargs="+",
                        help="Genomic positions/ranges like chr1:100-120 chrX:2001 (1-based, inclusive).")
    parser.add_argument("--bed", help="BED file of regions (0-based, half-open).")
    parser.add_argument("--threshold", type=float, default=0.5,
                        help="Probability threshold (0..1) to call a read 'modified' (default: 0.5).")
    parser.add_argument("--perpos-output", required=True,
                        help="Output CSV path for per-position aggregated metrics.")
    parser.add_argument("--region-output",
                        help="Write region-level summary CSV (uses --bed and/or --positions as regions).")
    parser.add_argument("--per-read-output",
                        help="Optional: stream a per-read wide table (only in read mode AND with regions).")
    parser.add_argument("--mode", choices=["auto", "pileup", "read"], default="auto",
                        help="auto=pileup if regions else read; or force one.")
    parser.add_argument(
        "--region-fraction-mode",
        choices=["weighted", "unweighted", "both"],
        default="weighted",
        help="Aggregate per-site fractions within regions as coverage-weighted (modified_sum/covered_sum), "
             "unweighted mean of per-site fractions, or both."
    )
    args = parser.parse_args()

    if not (0.0 <= args.threshold <= 1.0):
        raise ValueError("--threshold must be between 0 and 1.")

    ensure_bam_index_exists(args.bam)

    # Collect regions from both sources
    ranges: List[Tuple[str, int, int]] = []
    ranges.extend(parse_positions_cli(args.positions) if args.positions else [])
    if args.bed:
        ranges.extend(parse_bed(args.bed))

    # Sanity: if user supplied regions but parser yielded zero, stop loudly
    if (args.positions or args.bed) and not ranges:
        raise ValueError("No valid regions parsed from --positions/--bed. "
                         "Check contig names and formatting (e.g., 'NC_060925.1:1-5000').")

    # Optional: log what we parsed
    import sys
    if ranges:
        print(f"[info] Parsed {len(ranges)} region(s); first 3: {ranges[:3]}", file=sys.stderr, flush=True)
    else:
        print("[info] No regions provided; falling back to full-BAM iteration.", file=sys.stderr, flush=True)

    # Mode selection
    if args.mode == "pileup":
        if not ranges:
            raise ValueError("pileup mode requires --positions and/or --bed.")
        use_pileup = True
    elif args.mode == "read":
        use_pileup = False
    else:  # auto
        use_pileup = bool(ranges)

    # per-read output constraints
    per_read_writer = None
    per_read_positions: Optional[List[Tuple[str, int]]] = None
    pr_file = None
    if args.per_read_output:
        if use_pileup:
            raise ValueError("--per-read-output is only supported in read mode.")
        if not ranges:
            raise ValueError("--per-read-output requires --positions and/or --bed to constrain columns.")
        per_read_positions = expand_positions(ranges)
        # write header
        with open(args.per_read_output, "w", newline="") as prf:
            prw = csv.writer(prf)
            header = ["Read Name"] + [f"{c}:{p}|{args.mod}" for (c, p) in per_read_positions]
            prw.writerow(header)
        # reopen for appending during streaming
        pr_file = open(args.per_read_output, "a", newline="")
        per_read_writer = csv.writer(pr_file)

    try:
        if use_pileup:
            agg = parse_bam_streaming_pileup(
                args.bam, ranges, args.mod, args.threshold
            )
        else:
            agg = parse_bam_streaming(
                args.bam, ranges, args.mod, args.threshold,
                per_read_writer=per_read_writer,
                per_read_positions=per_read_positions
            )
    finally:
        if pr_file:
            pr_file.close()

    # Per-position CSV (only positions we actually saw)
    all_positions = sorted(agg.keys())
    mod_label = args.mod.replace("+", "plus")
    with open(args.perpos_output, "w", newline="") as pf:
        pw = csv.writer(pf)
        pw.writerow([
            "Chrom", "Pos",
            "Covered_reads",
            f"Modified_reads(thresh={args.threshold})",
            f"Fraction_{mod_label}",
            f"AverageProb_{mod_label}",
            f"Variance_{mod_label}",
            "Num_prob_observations"
        ])
        for (chrom, pos) in all_positions:
            pa = agg[(chrom, pos)]
            frac = (pa.modified / pa.covered) if pa.covered > 0 else 0.0
            pw.writerow([
                chrom, pos,
                pa.covered,
                pa.modified,
                f"{frac:.3f}",
                f"{pa.mean:.3f}",
                f"{pa.variance():.3f}",
                pa.n_probs
            ])

    # ----- Region CSV (BED and/or POSITIONS) -----
    if args.region_output:
        region_ranges: List[Tuple[str, int, int]] = []
        if args.bed:
            region_ranges.extend(parse_bed(args.bed))
        if args.positions:
            region_ranges.extend(parse_positions_cli(args.positions))

        if not region_ranges:
            raise ValueError("--region-output requires at least one region from --bed or --positions.")

        regions = to_regions(region_ranges)
        region_rows = summarize_regions(regions, agg, args.mod, args.threshold, mode=args.region_fraction_mode)
        with open(args.region_output, "w", newline="") as rf:
            # field order: include whichever fraction(s) were requested
            sample_row = region_rows[0] if region_rows else {}
            fieldnames = list(sample_row.keys()) if sample_row else [
                "Region",
                f"Region_AverageProb_{mod_label}",
                "Covered_reads_sum", "Covered_reads_avg",
                "Modified_reads_sum", "Num_prob_observations",
                f"Region_Fraction_{mod_label}(thresh={args.threshold})"
            ]
            rw = csv.DictWriter(rf, fieldnames=fieldnames)
            rw.writeheader()
            for row in region_rows:
                rw.writerow(row)

if __name__ == "__main__":
    main()
