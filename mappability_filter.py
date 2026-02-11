#!/usr/bin/env python3
"""
mappability_filter.py

Filter genes and transcripts based on exonic mappability scores using
UCSC Hoffman lab mappability tracks.

Computes the fraction of exonic bases that are unmappable for each gene
and transcript, then outputs lists of IDs to remove from downstream analysis.

Supports:
  - Automatic download of UCSC mappability files for k-mer lengths: 24, 36, 50, 100
  - bigBed (.bb) — Umap Unique Mappability (intervals = uniquely mappable regions)
  - bigWig (.bw) — continuous per-base mappability scores

Requirements:
    pip install pyBigWig

Usage:
    # Auto-download k36 mappability and filter:
    python mappability_filter.py --kmer 36 --gtf annotations.gtf

    # Use local bigBed file:
    python mappability_filter.py --bigbed k36.Unique.Mappability.bb --gtf annotations.gtf

    # Use local bigWig file:
    python mappability_filter.py --bigwig mappability.bw --gtf annotations.gtf
"""

import argparse
import gzip
import os
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path

# UCSC Hoffman Mappability base URL
UCSC_MAPPABILITY_URL = "https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability"

# GENCODE GTF base URL
GENCODE_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human"

# Available k-mer lengths for unique mappability
AVAILABLE_KMERS = [24, 36, 50, 100]

# Available GENCODE versions (recent releases)
AVAILABLE_GENCODE_VERSIONS = list(range(38, 48))  # v38 to v47


def get_mappability_url(kmer: int, track_type: str = "Unique.Mappability") -> str:
    """Get the UCSC download URL for a mappability file."""
    return f"{UCSC_MAPPABILITY_URL}/k{kmer}.{track_type}.bb"


def get_gencode_url(version: int, annotation_type: str = "basic") -> str:
    """
    Get the GENCODE download URL for a GTF file.

    annotation_type: "basic" or "comprehensive"
    """
    if annotation_type == "basic":
        return f"{GENCODE_URL}/release_{version}/gencode.v{version}.basic.annotation.gtf.gz"
    else:
        return f"{GENCODE_URL}/release_{version}/gencode.v{version}.annotation.gtf.gz"


def download_file(url: str, output_path: Path, description: str) -> Path:
    """Download a file from URL if not already present."""
    if output_path.exists():
        print(f"[INFO] Using cached: {output_path}", file=sys.stderr)
        return output_path

    print(f"[INFO] Downloading {description}: {url}", file=sys.stderr)
    print(f"[INFO] Destination: {output_path}", file=sys.stderr)

    try:
        urllib.request.urlretrieve(url, output_path)
        print(f"[INFO] Download complete: {output_path}", file=sys.stderr)
    except Exception as e:
        print(f"[ERROR] Failed to download {url}: {e}", file=sys.stderr)
        sys.exit(1)

    return output_path


def download_mappability(kmer: int, output_dir: Path, track_type: str = "Unique.Mappability") -> Path:
    """Download a mappability file from UCSC if not already present."""
    filename = f"k{kmer}.{track_type}.bb"
    output_path = output_dir / filename
    url = get_mappability_url(kmer, track_type)
    return download_file(url, output_path, f"k{kmer} mappability")


def download_gencode(version: int, output_dir: Path, annotation_type: str = "basic") -> Path:
    """Download a GENCODE GTF file if not already present."""
    if annotation_type == "basic":
        filename = f"gencode.v{version}.basic.annotation.gtf.gz"
    else:
        filename = f"gencode.v{version}.annotation.gtf.gz"
    output_path = output_dir / filename
    url = get_gencode_url(version, annotation_type)
    return download_file(url, output_path, f"GENCODE v{version}")


# ---------------------------------------------------------------------------
# GTF parsing
# ---------------------------------------------------------------------------

def parse_gtf_exons(gtf_path: str):
    """
    Parse exon records from a GTF file (Ensembl or GENCODE).

    Returns:
        gene_exons:  dict  gene_id  -> list of (chrom, start, end)  [0-based half-open]
        tx_exons:    dict  tx_id    -> list of (chrom, start, end)
        tx_to_gene:  dict  tx_id    -> gene_id
    """
    gene_exons = defaultdict(list)
    tx_exons = defaultdict(list)
    tx_to_gene = {}

    opener = gzip.open if gtf_path.endswith(".gz") else open

    print(f"[INFO] Parsing GTF: {gtf_path}", file=sys.stderr)
    n = 0
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "exon":
                continue

            chrom = cols[0]
            start = int(cols[3]) - 1  # GTF 1-based -> 0-based
            end = int(cols[4])
            attrs = cols[8]

            gene_id = _extract_attr(attrs, "gene_id")
            tx_id = _extract_attr(attrs, "transcript_id")
            if gene_id is None or tx_id is None:
                continue

            # Strip version numbers (ENSG00000123456.1 -> ENSG00000123456)
            gene_id = gene_id.split(".")[0]
            tx_id = tx_id.split(".")[0]

            gene_exons[gene_id].append((chrom, start, end))
            tx_exons[tx_id].append((chrom, start, end))
            tx_to_gene[tx_id] = gene_id

            n += 1
            if n % 500_000 == 0:
                print(f"  ... {n:,} exon records", file=sys.stderr)

    print(f"[INFO] {n:,} exons — {len(gene_exons):,} genes, {len(tx_exons):,} transcripts",
          file=sys.stderr)
    return gene_exons, tx_exons, tx_to_gene


def _extract_attr(attr_str: str, key: str) -> str | None:
    """Extract an attribute value from a GTF attributes string."""
    needle = f'{key} "'
    idx = attr_str.find(needle)
    if idx == -1:
        return None
    start = idx + len(needle)
    end = attr_str.index('"', start)
    return attr_str[start:end]


# ---------------------------------------------------------------------------
# Interval utilities
# ---------------------------------------------------------------------------

def merge_intervals(intervals):
    """Merge list of (chrom, start, end) -> dict chrom -> [(start, end), ...]"""
    by_chrom = defaultdict(list)
    for chrom, s, e in intervals:
        by_chrom[chrom].append((s, e))
    merged = {}
    for chrom, ivs in by_chrom.items():
        ivs.sort()
        out = [ivs[0]]
        for s, e in ivs[1:]:
            if s <= out[-1][1]:
                out[-1] = (out[-1][0], max(out[-1][1], e))
            else:
                out.append((s, e))
        merged[chrom] = out
    return merged


def _resolve_chrom(chrom: str, available_chroms: dict) -> str | None:
    """Handle chr prefix mismatch between Ensembl GTF and UCSC files."""
    if chrom in available_chroms:
        return chrom
    if f"chr{chrom}" in available_chroms:
        return f"chr{chrom}"
    if chrom.startswith("chr") and chrom[3:] in available_chroms:
        return chrom[3:]
    return None


# ---------------------------------------------------------------------------
# bigBed scoring — Umap Unique Mappability
# ---------------------------------------------------------------------------

def _covered_bases_in_interval(mappable_ivs, query_start: int, query_end: int) -> int:
    """
    Given a list of (start, end, ...) tuples from bb.entries() that overlap
    [query_start, query_end), compute how many bases within the query are
    covered (i.e. mappable).
    """
    if mappable_ivs is None:
        return 0
    covered = 0
    for entry in mappable_ivs:
        m_start = entry[0]
        m_end = entry[1]
        ov_start = max(m_start, query_start)
        ov_end = min(m_end, query_end)
        if ov_start < ov_end:
            covered += ov_end - ov_start
    return covered


def score_intervals_bigbed(bb, merged_intervals: dict) -> tuple[int, int]:
    """
    Score mappability using a bigBed of uniquely-mappable regions.
    Returns (total_bases, unmappable_bases).
    """
    total = 0
    unmappable = 0
    bb_chroms = bb.chroms()

    for chrom, ivs in merged_intervals.items():
        query_chrom = _resolve_chrom(chrom, bb_chroms)

        if query_chrom is None:
            for s, e in ivs:
                total += e - s
                unmappable += e - s
            continue

        chrom_len = bb_chroms[query_chrom]

        for s, e in ivs:
            s_clip = max(0, s)
            e_clip = min(e, chrom_len)
            if s_clip >= e_clip:
                continue

            length = e_clip - s_clip
            total += length

            try:
                entries = bb.entries(query_chrom, s_clip, e_clip)
            except Exception:
                unmappable += length
                continue

            covered = _covered_bases_in_interval(entries, s_clip, e_clip)
            unmappable += length - covered

    return total, unmappable


# ---------------------------------------------------------------------------
# bigWig scoring — per-base mappability values
# ---------------------------------------------------------------------------

def score_intervals_bigwig(bw, merged_intervals: dict, mappability_threshold: float) -> tuple[int, int]:
    """
    Score mappability using a bigWig with per-base float scores.
    A base is unmappable if score < mappability_threshold (or NaN/missing).
    Returns (total_bases, unmappable_bases).
    """
    total = 0
    unmappable = 0
    bw_chroms = bw.chroms()

    for chrom, ivs in merged_intervals.items():
        query_chrom = _resolve_chrom(chrom, bw_chroms)

        if query_chrom is None:
            for s, e in ivs:
                total += e - s
                unmappable += e - s
            continue

        chrom_len = bw_chroms[query_chrom]

        for s, e in ivs:
            s_clip = max(0, s)
            e_clip = min(e, chrom_len)
            if s_clip >= e_clip:
                continue

            length = e_clip - s_clip
            total += length

            try:
                vals = bw.values(query_chrom, s_clip, e_clip)
            except Exception:
                unmappable += length
                continue

            for v in vals:
                if v is None or v != v:  # NaN check
                    unmappable += 1
                elif v < mappability_threshold:
                    unmappable += 1

    return total, unmappable


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Filter genes/transcripts with low exonic mappability.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"""
Available k-mer lengths: {', '.join(map(str, AVAILABLE_KMERS))}
Available GENCODE versions: {min(AVAILABLE_GENCODE_VERSIONS)}-{max(AVAILABLE_GENCODE_VERSIONS)}

Examples:
  # Auto-download both mappability and GTF:
  python mappability_filter.py --kmer 36 --gencode 47

  # Use local GTF with auto-downloaded mappability:
  python mappability_filter.py --kmer 36 --gtf annotations.gtf

  # Use all local files:
  python mappability_filter.py --bigbed k50.Unique.Mappability.bb --gtf annotations.gtf

  # Stricter filtering (flag if >=30% unmappable):
  python mappability_filter.py --kmer 36 --gencode 47 --frac-unmappable-cutoff 0.3
        """,
    )

    # Mappability input options (mutually exclusive)
    inp = parser.add_mutually_exclusive_group(required=True)
    inp.add_argument(
        "--kmer", "-k", type=int, choices=AVAILABLE_KMERS,
        help=f"K-mer length for auto-download from UCSC. Options: {AVAILABLE_KMERS}",
    )
    inp.add_argument(
        "--bigbed", "-bb",
        help="Local bigBed (.bb) file. Intervals = uniquely mappable regions.",
    )
    inp.add_argument(
        "--bigwig", "-bw",
        help="Local bigWig (.bw) file with per-base mappability scores.",
    )

    # GTF input options (mutually exclusive)
    gtf_group = parser.add_mutually_exclusive_group(required=True)
    gtf_group.add_argument(
        "--gtf", "-g",
        help="Local GTF annotation file (Ensembl/GENCODE, plain or .gz).",
    )
    gtf_group.add_argument(
        "--gencode", type=int, metavar="VERSION",
        help=f"GENCODE version to auto-download (v{min(AVAILABLE_GENCODE_VERSIONS)}-v{max(AVAILABLE_GENCODE_VERSIONS)}).",
    )
    parser.add_argument(
        "--gencode-comprehensive", action="store_true",
        help="Download comprehensive annotation instead of basic (larger file).",
    )

    # Output options
    parser.add_argument(
        "--out-dir", "-o", default=".",
        help="Output directory (default: current directory)",
    )
    parser.add_argument(
        "--out-prefix", "-p", default="",
        help="Prefix for output files (default: none)",
    )

    # Thresholds
    parser.add_argument(
        "--mappability-threshold", type=float, default=0.5,
        help="[bigWig only] Per-base score below this = unmappable (default: 0.5)",
    )
    parser.add_argument(
        "--frac-unmappable-cutoff", "-c", type=float, default=0.5,
        help="Flag gene/tx if >= this fraction of exonic bases are unmappable (default: 0.5)",
    )

    # Download options
    parser.add_argument(
        "--cache-dir", default=None,
        help="Directory to cache downloaded mappability files (default: output directory)",
    )

    args = parser.parse_args()

    # ---- Setup paths ----
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cache_dir = Path(args.cache_dir) if args.cache_dir else out_dir
    cache_dir.mkdir(parents=True, exist_ok=True)

    prefix = f"{args.out_prefix}_" if args.out_prefix else ""

    out_genes = out_dir / f"{prefix}genes_to_remove.txt"
    out_tx = out_dir / f"{prefix}transcripts_to_remove.txt"
    out_full = out_dir / f"{prefix}mappability_scores.tsv"

    # ---- Validate GENCODE version ----
    if args.gencode and args.gencode not in AVAILABLE_GENCODE_VERSIONS:
        print(f"[ERROR] GENCODE version {args.gencode} not in available range "
              f"({min(AVAILABLE_GENCODE_VERSIONS)}-{max(AVAILABLE_GENCODE_VERSIONS)})", file=sys.stderr)
        sys.exit(1)

    # ---- Import pyBigWig ----
    try:
        import pyBigWig
    except ImportError:
        print("[ERROR] pyBigWig required: pip install pyBigWig", file=sys.stderr)
        sys.exit(1)

    # ---- Get GTF file ----
    if args.gencode:
        annotation_type = "comprehensive" if args.gencode_comprehensive else "basic"
        gtf_path = str(download_gencode(args.gencode, cache_dir, annotation_type))
    else:
        gtf_path = args.gtf

    # ---- Get mappability file ----
    if args.kmer:
        map_path = download_mappability(args.kmer, cache_dir)
        use_bigbed = True
        kmer_label = f"k{args.kmer}"
    elif args.bigbed:
        map_path = Path(args.bigbed)
        use_bigbed = True
        kmer_label = map_path.stem
    else:
        map_path = Path(args.bigwig)
        use_bigbed = False
        kmer_label = map_path.stem

    # ---- Open mappability file ----
    print(f"[INFO] Opening: {map_path}", file=sys.stderr)
    fh_map = pyBigWig.open(str(map_path))
    if fh_map is None:
        print(f"[ERROR] Could not open: {map_path}", file=sys.stderr)
        sys.exit(1)

    is_bigbed = fh_map.isBigBed()
    is_bigwig = fh_map.isBigWig()

    if use_bigbed and not is_bigbed:
        print(f"[WARN] Expected bigBed but file reports isBigBed()=False. Proceeding anyway.",
              file=sys.stderr)
    if not use_bigbed and not is_bigwig:
        print(f"[WARN] Expected bigWig but file reports isBigWig()=False. Proceeding anyway.",
              file=sys.stderr)

    mode_str = "bigBed (binary unique mappability)" if use_bigbed else "bigWig (continuous scores)"
    print(f"[INFO] Mode: {mode_str}", file=sys.stderr)

    # ---- Parse GTF ----
    gene_exons, tx_exons, tx_to_gene = parse_gtf_exons(gtf_path)

    # ---- Define scoring function ----
    def score_fn(merged):
        if use_bigbed:
            return score_intervals_bigbed(fh_map, merged)
        else:
            return score_intervals_bigwig(fh_map, merged, args.mappability_threshold)

    # ---- Score transcripts ----
    print("[INFO] Scoring transcripts...", file=sys.stderr)
    tx_scores = {}
    done = 0
    for tx_id, exons in tx_exons.items():
        merged = merge_intervals(exons)
        total, unmappable = score_fn(merged)
        frac = unmappable / total if total > 0 else 0.0
        tx_scores[tx_id] = (total, unmappable, frac)
        done += 1
        if done % 10_000 == 0:
            print(f"  ... {done:,} / {len(tx_exons):,} transcripts", file=sys.stderr)

    # ---- Score genes ----
    print("[INFO] Scoring genes...", file=sys.stderr)
    gene_scores = {}
    for gene_id, exons in gene_exons.items():
        merged = merge_intervals(exons)
        total, unmappable = score_fn(merged)
        frac = unmappable / total if total > 0 else 0.0
        gene_scores[gene_id] = (total, unmappable, frac)

    fh_map.close()

    # ---- Write full table ----
    print(f"[INFO] Writing: {out_full}", file=sys.stderr)
    with open(out_full, "w") as fh:
        fh.write("level\tid\tgene_id\ttotal_exonic_bases\tunmappable_bases\tfrac_unmappable\n")
        for gene_id in sorted(gene_scores):
            total, unmappable, frac = gene_scores[gene_id]
            fh.write(f"gene\t{gene_id}\t{gene_id}\t{total}\t{unmappable}\t{frac:.6f}\n")
        for tx_id in sorted(tx_scores):
            total, unmappable, frac = tx_scores[tx_id]
            gene_id = tx_to_gene.get(tx_id, "NA")
            fh.write(f"transcript\t{tx_id}\t{gene_id}\t{total}\t{unmappable}\t{frac:.6f}\n")

    # ---- Write removal lists ----
    cutoff = args.frac_unmappable_cutoff

    genes_remove = sorted(
        gid for gid, (t, u, f) in gene_scores.items() if f >= cutoff
    )
    tx_remove = sorted(
        tid for tid, (t, u, f) in tx_scores.items() if f >= cutoff
    )

    with open(out_genes, "w") as fh:
        fh.write(f"# Genes with >= {cutoff:.0%} unmappable exonic bases ({kmer_label} unique mappability)\n")
        fh.write(f"# Input: {map_path}\n")
        fh.write(f"# GTF: {gtf_path}\n")
        fh.write(f"# Total genes scored: {len(gene_scores)}\n")
        fh.write(f"# Genes flagged for removal: {len(genes_remove)}\n")
        fh.write("# gene_id\tfrac_unmappable\n")
        for gid in genes_remove:
            total, unmappable, frac = gene_scores[gid]
            fh.write(f"{gid}\t{frac:.4f}\n")

    with open(out_tx, "w") as fh:
        fh.write(f"# Transcripts with >= {cutoff:.0%} unmappable exonic bases ({kmer_label} unique mappability)\n")
        fh.write(f"# Input: {map_path}\n")
        fh.write(f"# GTF: {gtf_path}\n")
        fh.write(f"# Total transcripts scored: {len(tx_scores)}\n")
        fh.write(f"# Transcripts flagged for removal: {len(tx_remove)}\n")
        fh.write("# transcript_id\tgene_id\tfrac_unmappable\n")
        for tid in tx_remove:
            total, unmappable, frac = tx_scores[tid]
            gene_id = tx_to_gene.get(tid, "NA")
            fh.write(f"{tid}\t{gene_id}\t{frac:.4f}\n")

    # ---- Summary ----
    print(f"\n{'=' * 60}", file=sys.stderr)
    print("SUMMARY", file=sys.stderr)
    print(f"{'=' * 60}", file=sys.stderr)
    print(f"Mappability file:  {map_path}", file=sys.stderr)
    print(f"GTF file:          {gtf_path}", file=sys.stderr)
    print(f"Mode:              {'bigBed' if use_bigbed else 'bigWig'}", file=sys.stderr)
    if not use_bigbed:
        print(f"Mappability thr:   {args.mappability_threshold}", file=sys.stderr)
    print(f"Removal cutoff:    {cutoff:.0%} unmappable exonic bases", file=sys.stderr)
    print(f"Genes scored:      {len(gene_scores):>8,}", file=sys.stderr)
    print(f"Genes to remove:   {len(genes_remove):>8,}  "
          f"({100 * len(genes_remove) / max(len(gene_scores), 1):.1f}%)", file=sys.stderr)
    print(f"Tx scored:         {len(tx_scores):>8,}", file=sys.stderr)
    print(f"Tx to remove:      {len(tx_remove):>8,}  "
          f"({100 * len(tx_remove) / max(len(tx_scores), 1):.1f}%)", file=sys.stderr)
    print(f"\nOutputs:", file=sys.stderr)
    print(f"  {out_genes}", file=sys.stderr)
    print(f"  {out_tx}", file=sys.stderr)
    print(f"  {out_full}", file=sys.stderr)
    print(f"{'=' * 60}", file=sys.stderr)


if __name__ == "__main__":
    main()
