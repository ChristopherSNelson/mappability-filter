# CLAUDE.md

This file provides context for Claude Code when working on this project.

## Project Overview

A Python tool to filter genes/transcripts with low exonic mappability for RNA-seq and Ribo-seq analysis. Uses UCSC Hoffman lab mappability tracks to identify genes where a high fraction of exonic bases fall in multi-mapping regions.

## Key Files

- `mappability_filter.py` - Main script. Handles GTF parsing, mappability scoring, and output generation.
- `download_data.sh` - Shell script for manual data download (mappability + GTF files).
- `output/` - Pre-generated removal lists for GENCODE v48 protein-coding genes at k24/36/50/100.

## Architecture

The script supports two input modes:
- **bigBed** (.bb) - Binary intervals representing uniquely mappable regions. Bases NOT covered = unmappable.
- **bigWig** (.bw) - Per-base float scores. Bases below threshold = unmappable.

Flow: Parse GTF exons → Merge overlapping intervals → Query mappability track → Calculate unmappable fraction → Output removal lists.

## Dependencies

- `pyBigWig` - Only external dependency. Used to read bigBed/bigWig files.

## Data Sources

- Mappability: https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/
- GTF: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/

## Common Tasks

```bash
# Run with auto-download
python mappability_filter.py --kmer 36 --gencode 47

# Run with local files
python mappability_filter.py --bigbed k36.Unique.Mappability.bb --gtf annotations.gtf

# Test both modes
python mappability_filter.py --bigbed k36.Unique.Mappability.bb --gtf file.gtf --out-prefix test_bb
python mappability_filter.py --bigwig k36.Umap.MultiTrackMappability.bw --gtf file.gtf --out-prefix test_bw
```

## Conventions

- Gene/transcript IDs are stripped of version numbers (ENSG00000123456.1 → ENSG00000123456)
- Chromosome names handle both UCSC (chr1) and Ensembl (1) formats
- All coordinates are 0-based half-open internally (GTF 1-based converted on parse)
