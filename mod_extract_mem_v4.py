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



#!/usr/bin/env python3
import argparse, csv, re, sys, math
from collections import defaultdict

import pysam
import pandas as pd

# ---------------------------
# Region parsing / utilities
# ---------------------------

def parse_positions(position_args):
    """
    Accepts:
      - chr:start-end
      - chr:pos
    Returns list of (chrom, start, end) 1-based inclusive
    """
    out = []
    for pos in position_args or []:
        m = re.match(r"^([\w\.\-]+):(\d+)-(\d+)$", pos.replace(",", ""))
        if m:
            c, s, e = m.groups()
            s, e = int(s), int(e)
            if e < s: s, e = e, s
            out.append((c, s, e))
            continue
        m = re.match(r"^([\w\.\-]+):(\d+)$", pos.replace(",", ""))
        if m:
            c, p = m.groups()
            p = int(p)
            out.append((c, p, p))
            continue
        raise SystemExit(f"Bad --positions entry: {pos}")
    return out

def load_bed_regions(path):
    """
    Simple BED loader (0-based, half-open). Convert to 1-based inclusive.
    Columns: chrom start end [name...]
    """
    regs = []
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            regs.append((chrom, start + 1, end))
    return regs

def merge_intervals(regions):
    """
    Merge overlapping/adjacent intervals per chromosome.
    Input/Output: list of (chrom, start, end), 1-based inclusive.
    """
    from itertools import groupby
    out = []
    regions = sorted(regions, key=lambda x: (x[0], x[1], x[2]))
    for chrom, grp in groupby(regions, key=lambda x: x[0]):
        cur = []
        for _, s, e in grp:
            if not cur:
                cur = [chrom, s, e]
            else:
                if s <= cur[2] + 1:
                    cur[2] = max(cur[2], e)
                else:
                    out.append(tuple(cur))
                    cur = [chrom, s, e]
        if cur:
            out.append(tuple(cur))
    return out

def region_overlaps_read(chrom, read_start_1based, read_end_1based, include_regions):
    # quick overlap test for a read against any region with same chrom
    for c, s, e in include_regions:
        if c != chrom: continue
        if read_end_1based >= s and read_start_1based <= e:
            return True
    return False

# ---------------------------
# Basemod helpers
# ---------------------------

def parse_basemod_arg(basemod_str):
    s = (basemod_str or "CG").replace(" ", "").upper()
    if s in ("C+M", "C", "CG", "MCCG"):
        return {"CG"}
    if s in ("A+A", "A", "MA", "M6A"):
        return {"A"}
    if s in ("A+CG", "CG+A", "A+C", "C+A"):
        return {"A", "CG"}
    raise SystemExit(f"Unsupported --basemod '{basemod_str}'. Use 'CG', 'A', or 'A+CG'.")

# ---------------------------
# CpG index for DiMeLo mode
# ---------------------------

def load_cpg_sites_from_bed(bed_path, include_regions):
    """
    Return dict: chrom -> sorted list (unique) of CpG C positions (1-based) within include_regions.
    Expect CpG BED intervals of width 2bp (C,G) in 0-based coordinates.
    """
    keep = defaultdict(list)
    with open(bed_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"): continue
            chrom, start, end = line.split()[:3]
            start, end = int(start), int(end)
            if end - start != 2:
                continue
            cpos = start + 1  # 1-based C position
            if any(chrom == c and s <= cpos <= e for c, s, e in include_regions):
                keep[chrom].append(cpos)
    for k in keep:
        keep[k] = sorted(set(keep[k]))
    return keep

def build_cpg_index_from_fasta(fasta_path, include_regions):
    """
    Scan FASTA in include_regions; return dict chrom -> sorted list of 1-based CpG C positions.
    """
    fa = pysam.FastaFile(fasta_path)
    idx = defaultdict(list)
    for chrom, s, e in include_regions:
        seq = fa.fetch(chrom, s - 1, e + 1).upper()  # pysam: 0-based, end-exclusive
        for i in range(0, len(seq) - 1):
            if seq[i] == "C" and seq[i + 1] == "G":
                idx[chrom].append(s + i)
    for k in idx:
        idx[k] = sorted(set(idx[k]))
    return idx

# ---------------------------
# Robust modified_bases normalization
# ---------------------------

def _collect_qprobs_from_modified_bases(mb, c_codes=None, a_codes=None):
    """
    Normalize pysam read.modified_bases -> two dicts:
      qprobC[qpos] = max prob among accepted C-mod entries at that qpos
      qprobA[qpos] = max prob among accepted A-mod entries at that qpos

    Accepts many key encodings:
      - b"C+m", b"A+a"
      - (b"C", b"m"), (b"A", b"a")
      - str variants like "C+m"
    If c_codes is None -> accept ANY C mod; same for a_codes.
    """
    qprobA, qprobC = {}, {}

    if not mb:
        return qprobC, qprobA

    for k, entries in mb.items():
        # Normalize key to (base, code_str)
        base = None
        code = ""

        if isinstance(k, (bytes, bytearray)):
            ks = k.decode()
            if "+" in ks:
                base, code = ks.split("+", 1)
            else:
                base, code = ks, ""
        elif isinstance(k, tuple) and len(k) >= 1:
            kb = k[0]
            base = kb.decode() if isinstance(kb, (bytes, bytearray)) else str(kb)
            # some pysam builds use (base, strand_flag, code)
            if len(k) >= 3:
                kc = k[2]
            elif len(k) >= 2:
                kc = k[1]
            else:
                kc = ""
            code = kc.decode() if isinstance(kc, (bytes, bytearray)) else str(kc)

        else:
            ks = str(k)
            if "+" in ks:
                base, code = ks.split("+", 1)
            else:
                base, code = ks, ""

        base = (base or "").upper()
        code_l = (code or "").lower()

        # Decide whether to accept this key
        if base == "C":
            if (c_codes is None) or (code_l in c_codes) or (code_l in ("5mc", "5mc".lower())):
                for tup in entries:
                    if not tup: 
                        continue
                    qpos = tup[0]
                    prob = tup[1] if len(tup) > 1 else None
                    if prob is None:
                        continue
                    if (qpos not in qprobC) or (prob > qprobC[qpos]):
                        qprobC[qpos] = prob

        elif base == "A":
            if (a_codes is None) or (code_l in a_codes) or (code_l in ("6ma", "m6a")):
                for tup in entries:
                    if not tup:
                        continue
                    qpos = tup[0]
                    prob = tup[1] if len(tup) > 1 else None
                    if prob is None:
                        continue
                    if (qpos not in qprobA) or (prob > qprobA[qpos]):
                        qprobA[qpos] = prob

        else:
            continue

    return qprobC, qprobA

# ---------------------------
# DiMeLo-style counting (multi-mod, no MD tag needed)
# ---------------------------

def dimelo_count_positions(bam_path, include_regions, basemods, cpg_index,
                           threshA_raw, threshC_raw, c_codes=None, a_codes=None, probe_keys=0):
    """
    DiMeLo-like pooled counts per position (no MD tag needed).
      - For 'CG': only CpG C positions using cpg_index (coverage+calls from C-mods).
      - For 'A' : any aligned read base == 'A' (coverage) with A-mod calls.
    c_codes: set of accepted C mod codes (e.g., {'m'}). None => accept ANY C mod
    a_codes: set of accepted A mod codes (e.g., {'a'}). None => accept ANY A mod
    probe_keys: if >0, print up to N examples of modified_bases keys observed (debugging).
    """
    counts = defaultdict(lambda: [0, 0])  # (meth, total)
    cpg_sets = {chrom: set(poslist) for chrom, poslist in (cpg_index or {}).items()}

    seen_keys = set()
    printed = 0

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if not (read.has_tag("Mm") and read.has_tag("Ml")):
                continue

            chrom = bam.get_reference_name(read.reference_id)
            rs1 = read.reference_start + 1
            re1 = read.reference_end
            if not region_overlaps_read(chrom, rs1, re1, include_regions):
                continue

            # Probe keys if requested
            if probe_keys and printed < probe_keys:
                try:
                    mb_probe = read.modified_bases
                    if mb_probe:
                        for k in mb_probe.keys():
                            seen_keys.add(k)
                        printed += 1
                except AttributeError:
                    pass

            # Collect qprobs
            try:
                mb = read.modified_bases
            except AttributeError:
                mb = None

            qprobC, qprobA = _collect_qprobs_from_modified_bases(mb, c_codes=c_codes, a_codes=a_codes)

            # Iterate aligned pairs (no ref seq needed)
            for qpos, rpos in read.get_aligned_pairs(matches_only=True):
                if qpos is None or rpos is None:
                    continue
                gpos = rpos + 1
                if not any(chrom == c and s <= gpos <= e for c, s, e in include_regions):
                    continue

                # CG path
                if "CG" in basemods and cpg_sets and chrom in cpg_sets and gpos in cpg_sets[chrom]:
                    counts[("CG", chrom, gpos)][1] += 1
                    prob = qprobC.get(qpos)
                    if prob is not None and prob >= threshC_raw:
                        counts[("CG", chrom, gpos)][0] += 1

                # A path
                if "A" in basemods and read.query_sequence:
                    rb = read.query_sequence[qpos]
                    if rb is not None and rb.upper() == "A":
                        counts[("A", chrom, gpos)][1] += 1
                        prob = qprobA.get(qpos)
                        if prob is not None and prob >= threshA_raw:
                            counts[("A", chrom, gpos)][0] += 1

    if probe_keys and seen_keys:
        sys.stderr.write("[probe] modified_bases keys observed (unique up to 50):\n")
        for k in sorted(seen_keys, key=lambda x: str(x))[:50]:
            sys.stderr.write(f"  {repr(k)}\n")

    return counts

# ---------------------------
# Dimelo outputs (dtype-safe smoothing)
# ---------------------------

def write_dimelo_outputs_multi(counts_dict, out_prefix, write_bed=False, smooth=1000, min_periods=100):
    """
    counts_dict: {(mod, chrom, pos): (meth, total)}  with mod in {'CG','A'}
    Writes per-mod:
      - {out_prefix}_{mod}_perpos.csv
      - {out_prefix}_{mod}_smoothed.csv
      - {out_prefix}_{mod}.bed (optional)
    """
    mods = sorted({k[0] for k in counts_dict.keys()})
    for mod in mods:
        rows = []
        for (m, chrom, pos), (meth, tot) in counts_dict.items():
            if m != mod:
                continue
            frac = (meth / tot) if tot > 0 else float("nan")
            rows.append((chrom, pos, meth, tot, frac))

        if not rows:
            continue

        df = pd.DataFrame(
            rows, columns=["chrom", "pos", "methylated_bases", "total_bases", "frac"]
        )

        # enforce numeric types
        df["pos"] = pd.to_numeric(df["pos"], errors="coerce").astype("Int64")
        df["methylated_bases"] = pd.to_numeric(df["methylated_bases"], errors="coerce")
        df["total_bases"] = pd.to_numeric(df["total_bases"], errors="coerce")
        df["frac"] = pd.to_numeric(df["frac"], errors="coerce")

        df = df.dropna(subset=["pos"]).copy()
        df["pos"] = df["pos"].astype(int)

        df = df.sort_values(["chrom", "pos"])
        df.to_csv(f"{out_prefix}_{mod}_perpos.csv", index=False)

        # rolling smoothing per chromosome
        sm_list = []
        for chrom, sub in df.groupby("chrom", sort=False):
            tmp = sub[["pos", "frac"]].copy()
            if smooth <= 1 or len(tmp) == 0:
                tmp["frac_smoothed"] = tmp["frac"]
            else:
                tmp["frac_smoothed"] = (
                    tmp.rolling(window=smooth, min_periods=min_periods, on="pos")
                       .mean()["frac"]
                )
            tmp.insert(0, "chrom", chrom)
            sm_list.append(tmp[["chrom", "pos", "frac_smoothed"]])

        sm = pd.concat(sm_list, ignore_index=True) if sm_list else pd.DataFrame(columns=["chrom","pos","frac_smoothed"])
        sm = sm.sort_values(["chrom", "pos"])
        sm.to_csv(f"{out_prefix}_{mod}_smoothed.csv", index=False)

        if write_bed:
            bed = df.copy()
            bed["start"] = bed["pos"] - 1
            bed["end"]   = bed["pos"]
            bed[["chrom", "start", "end", "methylated_bases", "total_bases"]].to_csv(
                f"{out_prefix}_{mod}.bed", sep="\t", header=False, index=False
            )

# ---------------------------
# Region summaries
# ---------------------------

def compute_region_stats_from_counts(pos_counts, include_regions):
    """
    pos_counts: dict (chrom,pos) -> (modified, covered)
    Returns list of dicts per region for CSV.
    """
    rows = []
    for chrom, start, end in include_regions:
        perpos_fracs = []
        cov_sum = 0
        mod_sum = 0
        n_cov_pos = 0
        for pos in range(start, end + 1):
            rec = pos_counts.get((chrom, pos))
            if not rec:
                continue
            mod, cov = rec
            cov_sum += cov
            mod_sum += mod
            if cov > 0:
                perpos_fracs.append(mod / cov)
                n_cov_pos += 1
        frac_unw = (sum(perpos_fracs) / len(perpos_fracs)) if perpos_fracs else 0.0
        var_unw = (pd.Series(perpos_fracs).var(ddof=1) if len(perpos_fracs) > 1 else 0.0)
        cov_avg = (cov_sum / n_cov_pos) if n_cov_pos > 0 else 0.0
        rows.append({
            "Region": f"{chrom}:{start}-{end}",
            "Covered_reads_sum": cov_sum,
            "Covered_reads_avg": cov_avg,
            "Modified_reads_sum": mod_sum,
            "Region_Fraction_unweighted": frac_unw,
            "Region_Variance_unweighted": var_unw,
            "Num_fraction_positions": n_cov_pos
        })
    return rows

# ---------------------------
# Main
# ---------------------------

def main():
    ap = argparse.ArgumentParser(description="Memory-lean modified-base extractor with DiMeLo-compatible mode.")
    ap.add_argument("bam", help="Input BAM (indexed).")

    # Unified basemod selector (both modes)
    ap.add_argument("--basemod", type=str, default="CG",
                    help="Which modification(s) to count. Use 'CG', 'A', or 'A+CG' (aliases: 'C','C+m','A','A+a').")

    # Regions
    g = ap.add_argument_group("Regions")
    g.add_argument("--positions", nargs="+", help="One or more chr:start-end (or chr:pos).")
    g.add_argument("--bed", help="BED of regions (0-based, half-open).")
    g.add_argument("--merge", action="store_true", help="Merge overlapping/adjacent regions.")

    # DiMeLo-compatible mode
    dm = ap.add_argument_group("DiMeLo-compatible mode")
    dm.add_argument("--dimelo-mode", action="store_true",
                    help="DiMeLo-like counting (CpG for 'CG' and any A for 'A'), raw Ml thresholds, pooled per-position counts; optional rolling smoothing & BED.")
    dm.add_argument("--threshC-raw", type=int, default=129, help="Raw Ml threshold (0–255) for calling mC (DiMeLo default 129 ≈ 0.506).")
    dm.add_argument("--threshA-raw", type=int, default=129, help="Raw Ml threshold (0–255) for calling mA (DiMeLo default 129 ≈ 0.506).")
    dm.add_argument("--smooth", type=int, default=1000, help="Rolling window (bp) for smoothed fraction output (DiMeLo default 1000).")
    dm.add_argument("--min-periods", type=int, default=100, help="Min positions in window (DiMeLo default 100).")
    dm.add_argument("--fasta", type=str, help="Reference FASTA (for CpG detection if no --cpg-bed). Required if basemod includes 'CG' and no --cpg-bed.")
    dm.add_argument("--cpg-bed", type=str, help="BED of CpG sites (2bp windows). Faster than scanning FASTA.")
    dm.add_argument("--dimelo-out-prefix", default="dimelo_like", help="Prefix for DiMeLo-like outputs.")
    dm.add_argument("--write-bed", action="store_true", help="Write DiMeLo-like BED with methylated_bases and total_bases.")
    dm.add_argument("--c-mod-codes", default="m",
                    help="Comma-separated C modification codes to count (e.g., 'm' for 5mC, 'm,h' to include 5hmC). Use '*' to accept ANY C mods. Default: m")
    dm.add_argument("--a-mod-codes", default="a",
                    help="Comma-separated A modification codes to count (e.g., 'a' for m6A). Use '*' to accept ANY A mods. Default: a")
    dm.add_argument("--probe-mod-keys", type=int, default=0,
                    help="If >0, print up to N examples of modified_bases keys observed (for debugging key formats).")

    # Outputs (normal mode or optional extras in DiMeLo mode)
    ap.add_argument("--perpos-output", help="Write per-position CSV (mod,chrom,pos,modified,covered,fraction).")
    ap.add_argument("--region-output", help="Write per-region CSV (summary).")

    args = ap.parse_args()

    # collect regions
    regs = []
    if args.positions:
        regs.extend(parse_positions(args.positions))
    if args.bed:
        regs.extend(load_bed_regions(args.bed))
    if not regs:
        raise SystemExit("Provide --positions and/or --bed.")
    if args.merge:
        regs = merge_intervals(regs)

    basemods = parse_basemod_arg(args.basemod)

    # -------------------
    # DiMeLo-compatible path
    # -------------------
    if args.dimelo_mode:
        cpg_index = None
        if "CG" in basemods:
            if not args.cpg_bed and not args.fasta:
                raise SystemExit("In --dimelo-mode for 'CG', provide either --cpg-bed or --fasta.")
            cpg_index = load_cpg_sites_from_bed(args.cpg_bed, regs) if args.cpg_bed else build_cpg_index_from_fasta(args.fasta, regs)

        c_codes = None if args.c_mod_codes.strip() == "*" else set(x.strip().lower() for x in args.c_mod_codes.split(",") if x.strip())
        a_codes = None if args.a_mod_codes.strip() == "*" else set(x.strip().lower() for x in args.a_mod_codes.split(",") if x.strip())

        counts = dimelo_count_positions(
            bam_path=args.bam,
            include_regions=regs,
            basemods=basemods,
            cpg_index=cpg_index,
            threshA_raw=args.threshA_raw,
            threshC_raw=args.threshC_raw,
            c_codes=c_codes,
            a_codes=a_codes,
            probe_keys=args.probe_mod_keys
        )

        write_dimelo_outputs_multi(
            counts_dict=counts,
            out_prefix=args.dimelo_out_prefix,
            smooth=args.smooth,
            min_periods=args.min_periods,
            write_bed=args.write_bed
        )

        # Optional region summary from these counts (unweighted mean per-position fractions) — per mod
        if args.region_output:
            for mod in basemods:
                flat = {(chrom,pos):(meth,tot) for (m,chrom,pos),(meth,tot) in counts.items() if m==mod}
                rows = compute_region_stats_from_counts(flat, regs)
                tag = "CG" if mod=="CG" else "A"
                pd.DataFrame(rows).to_csv(f"{args.region_output.rsplit('.csv',1)[0]}_{tag}.csv", index=False)

        # Optional per-position combined CSV
        if args.perpos_output:
            with open(args.perpos_output, "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["mod","chrom","pos","modified","covered","fraction"])
                for (m, chrom, pos), (meth, tot) in sorted(counts.items()):
                    frac = (meth/tot) if tot>0 else ""
                    w.writerow([m, chrom, pos, meth, tot, f"{frac:.6f}" if frac!="" else ""])
        print("[dimelo-mode] Done.")
        return

    # -------------------
    # Normal memory-lean path (per-position/per-region) honoring basemod (threshold=129 default)
    # -------------------
    pos_counts_by_mod = {k: defaultdict(lambda: [0,0]) for k in ("CG","A")}
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if not (read.has_tag("Mm") and read.has_tag("Ml")):
                continue
            chrom = bam.get_reference_name(read.reference_id)
            if not region_overlaps_read(chrom, read.reference_start + 1, read.reference_end, regs):
                continue

            try:
                mb = read.modified_bases
            except AttributeError:
                mb = None

            # For normal mode, accept C+m and A+a only (can be extended similarly if you want)
            qprobC, qprobA = _collect_qprobs_from_modified_bases(mb, c_codes={"m"}, a_codes={"a"})

            for qpos, rpos in read.get_aligned_pairs(matches_only=True):
                if qpos is None or rpos is None:
                    continue
                gpos = rpos + 1
                if not any(chrom == c and s <= gpos <= e for c, s, e in regs):
                    continue

                if "CG" in basemods:
                    pos_counts_by_mod["CG"][(chrom, gpos)][1] += 1
                    prob = qprobC.get(qpos)
                    if prob is not None and prob >= 129:
                        pos_counts_by_mod["CG"][(chrom, gpos)][0] += 1

                if "A" in basemods and read.query_sequence:
                    rb = read.query_sequence[qpos]
                    if rb is not None and rb.upper() == "A":
                        pos_counts_by_mod["A"][(chrom, gpos)][1] += 1
                        prob = qprobA.get(qpos)
                        if prob is not None and prob >= 129:
                            pos_counts_by_mod["A"][(chrom, gpos)][0] += 1

    # Outputs
    if args.perpos_output:
        with open(args.perpos_output, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["mod","chrom","pos","modified","covered","fraction"])
            for mod in [m for m in ("CG","A") if m in basemods]:
                for (chrom, pos), (modc, cov) in sorted(pos_counts_by_mod[mod].items()):
                    frac = (modc/cov) if cov>0 else ""
                    w.writerow([mod, chrom, pos, modc, cov, f"{frac:.6f}" if frac!="" else ""])

    if args.region_output:
        for mod in [m for m in ("CG","A") if m in basemods]:
            rows = compute_region_stats_from_counts(pos_counts_by_mod[mod], regs)
            tag = "CG" if mod=="CG" else "A"
            pd.DataFrame(rows).to_csv(f"{args.region_output.rsplit('.csv',1)[0]}_{tag}.csv", index=False)

    print("[standard mode] Done.")

if __name__ == "__main__":
    main()

