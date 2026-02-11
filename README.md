# Mappability Filter

Filter genes and transcripts based on exonic mappability for RNA-seq and Ribo-seq analysis.

This tool computes the fraction of exonic bases that fall in low-mappability regions and flags genes/transcripts that may produce unreliable quantification due to multi-mapping reads.

## Features

- **Auto-download** UCSC Hoffman lab mappability tracks (k24, k36, k50, k100)
- Supports **bigBed** (unique mappability) and **bigWig** (continuous scores)
- Works with **Ensembl** and **GENCODE** GTF annotations
- Configurable thresholds for flexible filtering

## Installation

```bash
pip install pyBigWig
```

## Quick Start

```bash
# Download data files (mappability + GTF)
./download_data.sh 36 48    # k36 + GENCODE v48
./download_data.sh          # all k-mers + GENCODE v48

# Run the filter
python mappability_filter.py --kmer 36 --gtf gencode.v48.basic.annotation.gtf.gz

# Or auto-download on the fly
python mappability_filter.py --kmer 36 --gencode 47

# Stricter filtering (flag if >=30% unmappable)
python mappability_filter.py --kmer 36 --gencode 47 --frac-unmappable-cutoff 0.3
```

## Usage

```
usage: mappability_filter.py [-h] (--kmer {24,36,50,100} | --bigbed BIGBED | --bigwig BIGWIG)
                             (--gtf GTF | --gencode VERSION) [--gencode-comprehensive]
                             [--out-dir OUT_DIR] [--out-prefix OUT_PREFIX]
                             [--mappability-threshold MAPPABILITY_THRESHOLD]
                             [--frac-unmappable-cutoff FRAC_UNMAPPABLE_CUTOFF]
                             [--cache-dir CACHE_DIR]

options:
  -h, --help            show this help message and exit

Mappability input (choose one):
  --kmer, -k {24,36,50,100}
                        K-mer length for auto-download from UCSC
  --bigbed, -bb BIGBED  Local bigBed (.bb) file
  --bigwig, -bw BIGWIG  Local bigWig (.bw) file

GTF input (choose one):
  --gtf, -g GTF         Local GTF file (Ensembl/GENCODE, plain or .gz)
  --gencode VERSION     GENCODE version to auto-download (v38-v47)
  --gencode-comprehensive
                        Download comprehensive annotation instead of basic

Output options:
  --out-dir, -o         Output directory (default: current directory)
  --out-prefix, -p      Prefix for output files
  --cache-dir           Directory to cache downloaded files

Thresholds:
  --mappability-threshold
                        [bigWig only] Per-base score threshold (default: 0.5)
  --frac-unmappable-cutoff, -c
                        Flag if >= this fraction unmappable (default: 0.5)
```

## Choosing K-mer Length

Match the k-mer to your read length:

| Read Length | Recommended K-mer |
|-------------|-------------------|
| 25-30 bp    | k24               |
| 30-40 bp    | k36               |
| 50-75 bp    | k50               |
| 100+ bp     | k100              |

For **Ribo-seq** (typically 28-32 nt footprints), use **k36**.

## Output Files

| File | Description |
|------|-------------|
| `genes_to_remove.txt` | Gene IDs to filter from analysis |
| `transcripts_to_remove.txt` | Transcript IDs to filter |
| `mappability_scores.tsv` | Full scores for custom filtering |

## Pre-generated Results

The `output/` directory contains pre-generated removal lists for **GENCODE v48 protein-coding genes** (hg38) at the default 50% unmappable threshold:

| K-mer | Genes to Remove | Transcripts to Remove |
|-------|-----------------|----------------------|
| k24   | 821 (4.1%)      | 2,923 (2.6%)         |
| k36   | 729 (3.6%)      | 2,247 (2.0%)         |
| k50   | 646 (3.2%)      | 1,922 (1.7%)         |
| k100  | 494 (2.5%)      | 1,423 (1.3%)         |

Shorter k-mers flag more genes because shorter reads have more multi-mapping potential.

## How It Works

1. Parses exon coordinates from the GTF file
2. For each gene/transcript, merges overlapping exons
3. Queries the mappability track to count unmappable bases
4. Flags entries where unmappable fraction ≥ cutoff (default 50%)

**bigBed mode**: Intervals in the file represent uniquely mappable regions. Bases *not* covered are unmappable.

**bigWig mode**: Per-base scores; bases below threshold (default 0.5) are unmappable.

## Data Sources

### Mappability Tracks (UCSC)

[Browse all mappability files](https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/)

Download with curl:
```bash
# k24 (130 MB) - for ~25-30 bp reads
curl -O https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k24.Unique.Mappability.bb

# k36 (67 MB) - for Ribo-seq / ~30-40 bp reads
curl -O https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k36.Unique.Mappability.bb

# k50 (38 MB) - for ~50-75 bp reads
curl -O https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k50.Unique.Mappability.bb

# k100 (6 MB) - for 100+ bp reads
curl -O https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k100.Unique.Mappability.bb
```

### GENCODE GTF Annotations

| Version | Release Date | Links |
|---------|--------------|-------|
| v47 | 2024-10 | [Basic](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.basic.annotation.gtf.gz) · [Comprehensive](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz) |
| v46 | 2024-05 | [Basic](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.basic.annotation.gtf.gz) · [Comprehensive](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz) |
| v45 | 2024-01 | [Basic](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.basic.annotation.gtf.gz) · [Comprehensive](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz) |
| v44 | 2023-08 | [Basic](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.basic.annotation.gtf.gz) · [Comprehensive](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz) |
| v43 | 2023-02 | [Basic](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.basic.annotation.gtf.gz) · [Comprehensive](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz) |

[Browse all GENCODE releases](https://www.gencodegenes.org/human/releases.html)

Download with curl:
```bash
# GENCODE v47 basic (recommended)
curl -O https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.basic.annotation.gtf.gz

# GENCODE v47 comprehensive (includes all transcripts)
curl -O https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz
```

**Basic vs Comprehensive**: Basic annotation excludes readthrough transcripts and other complex loci. Use basic unless you need everything.

## Example Output

```
============================================================
SUMMARY
============================================================
Mappability file:  k36.Unique.Mappability.bb
GTF file:          gencode.v48.annotation.gtf.gz
Mode:              bigBed
Removal cutoff:    50% unmappable exonic bases
Genes scored:        62,703
Genes to remove:        847  (1.4%)
Tx scored:          252,891
Tx to remove:        4,231  (1.7%)
============================================================
```

## License

MIT
