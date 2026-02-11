#!/bin/bash
# Download mappability and GTF data files
#
# Usage:
#   ./download_data.sh              # Download all k-mers + GENCODE v48
#   ./download_data.sh 36           # Download only k36
#   ./download_data.sh 36 47        # Download k36 + GENCODE v47

set -e

UCSC_BASE="https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability"
GENCODE_BASE="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human"

KMER=${1:-all}
GENCODE_VERSION=${2:-48}

download() {
    local url=$1
    local dest=$2
    if [ -f "$dest" ]; then
        echo "✓ Already exists: $dest"
    else
        echo "↓ Downloading: $dest"
        curl -L -o "$dest" "$url"
        echo "✓ Downloaded: $dest"
    fi
}

# Download mappability files
if [ "$KMER" = "all" ]; then
    for k in 24 36 50 100; do
        download "$UCSC_BASE/k${k}.Unique.Mappability.bb" "k${k}.Unique.Mappability.bb"
    done
else
    download "$UCSC_BASE/k${KMER}.Unique.Mappability.bb" "k${KMER}.Unique.Mappability.bb"
fi

# Download GENCODE GTF
GTF_FILE="gencode.v${GENCODE_VERSION}.basic.annotation.gtf.gz"
download "$GENCODE_BASE/release_${GENCODE_VERSION}/$GTF_FILE" "$GTF_FILE"

echo ""
echo "Done! Data files downloaded."
echo ""
echo "Example usage:"
if [ "$KMER" = "all" ]; then
    echo "  python mappability_filter.py --kmer 36 --gtf $GTF_FILE"
else
    echo "  python mappability_filter.py --kmer $KMER --gtf $GTF_FILE"
fi
