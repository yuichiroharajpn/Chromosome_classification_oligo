#!/usr/bin/env bash
# Download the assembly and matching NCBI report beside this script.
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
command -v wget >/dev/null || { echo "Error: wget is required" >&2; exit 1; }
assembly="GCF_003957565.2_bTaeGut1.4.pri"
base_url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/957/565/$assembly"
for suffix in genomic.fna.gz assembly_report.txt; do
    wget --continue --directory-prefix="$script_dir" "$base_url/${assembly}_${suffix}"
done
