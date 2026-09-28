#!/usr/bin/env bash
# All execution control lives here. Python is used only to parse the report.
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    echo "Usage: $(basename "$0") ASSEMBLY.fna.gz [ASSEMBLY_REPORT.txt] [--output-dir DIR] [--create-mid | --create-mid-only] [--cores 1-8]"
    echo "  Default: full chromosomes only; --create-mid: full and mid; --create-mid-only: mid only."
}
fail() { echo "Error: $*" >&2; exit 1; }
absolute_file() {
    local directory filename
    directory=$(dirname -- "$1")
    filename=$(basename -- "$1")
    printf '%s/%s\n' "$(cd -- "$directory" && pwd -P)" "$filename"
}
# Preserve the original Perl floating-point multiplication and truncation.
mid_region() {
    awk -v accession="$1" -v sequence_length="$2" 'BEGIN {
        printf "%s:%.0f-%.0f\n", accession, int(sequence_length * 0.05), int(sequence_length * 0.95)
    }'
}
# Each worker writes only files belonging to its k value.
kmer_worker() {
    local stage=$1 k=$2 path=$3 input_prefix=$4 output_suffix=$5 database
    case "$stage" in
        count)
            database="${path%.fa.gz}_${k}mer.jf"
            gzip -dc "$path" | jellyfish count -m "$k" -s 10000000 /dev/fd/0 -o "$database"
            jellyfish dump "$database" -c -t > "${database%.jf}.txt" ;;
        summary)
            Rscript "$script_dir/summary_kmerfreq.R" "$k" "${input_prefix}${k}mer" "${k}mer${output_suffix}" ;;
        analyze)
            Rscript "$script_dir/analyze_freq.R" "${k}mer${output_suffix}_relfreq.txt" ;;
        *) fail "unknown worker stage: $stage" ;;
    esac
}
run_kmers() {
    local stage=$1 path=$2 input_prefix=$3 output_suffix=$4 cores=$5 k
    if [[ "$cores" == 1 ]]; then
        for k in 2 3 4 5 6 7 8 9; do
            kmer_worker "$stage" "$k" "$path" "$input_prefix" "$output_suffix"
        done
    else
        # --halt 1 is the backward-compatible spelling of soon,fail=1.
        # Explicit Bash workers retain pipefail. Quote each argument, including paths with spaces.
        for k in 2 3 4 5 6 7 8 9; do
            printf '%q ' bash "$script_dir/run_01_kmer_classification.sh" --internal-worker "$stage" "$k" "$path" "$input_prefix" "$output_suffix"
            printf '\n'
        done | SHELL="$(command -v bash)" parallel --jobs "$cores" --keep-order --halt 1
    fi
}
main() {
    local assembly="" report="" output="$PWD" argument executable
    local fasta name accession length region path database k list
    local mode="full" input_prefix="" output_suffix="" cores=1
    local -a positional=()
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -h|--help) usage; return 0 ;;
            --create-mid|--create-mid-only)
                [[ "$mode" == "full" || "$mode" == "$1" ]] || fail "--create-mid and --create-mid-only are mutually exclusive"
                mode=$1; shift ;;
            --cores)
                [[ $# -ge 2 && "$2" =~ ^[1-8]$ ]] || fail "--cores must be an integer from 1 to 8"
                cores=$2; shift 2 ;;
            --output-dir)
                [[ $# -ge 2 && -n "$2" ]] || fail "--output-dir requires a directory"
                output=$2; shift 2 ;;
            --) shift; positional+=("$@"); break ;;
            -*) fail "unknown option: $1" ;;
            *) positional+=("$1"); shift ;;
        esac
    done
    [[ ${#positional[@]} -ge 1 && ${#positional[@]} -le 2 ]] || { usage >&2; return 1; }
    assembly=${positional[0]}
    [[ -f "$assembly" && "$assembly" == *.fna.gz ]] || fail "assembly must be an existing .fna.gz file"
    assembly=$(absolute_file "$assembly")
    report=${positional[1]:-"${assembly%_genomic.fna.gz}_assembly_report.txt"}
    [[ -f "$report" ]] || fail "assembly report not found: $report; supply it as the second argument"
    report=$(absolute_file "$report")
    for executable in python3 awk gzip bgzip samtools jellyfish Rscript; do
        command -v "$executable" >/dev/null || fail "required executable not found: $executable"
    done
    if [[ "$cores" != 1 ]]; then
        command -v parallel >/dev/null || fail "GNU parallel is required for --cores > 1"
        [[ "$(parallel --version)" == *"GNU parallel"* ]] || fail "parallel must be GNU parallel"
    fi
    Rscript -e 'for (p in c("ape","phangorn","dplyr","tibble","ca","readr","stringr")) if (!requireNamespace(p,quietly=TRUE)) stop(paste("Missing R package:",p))'
    mkdir -p -- "$output"
    output=$(cd -- "$output" && pwd -P)
    [[ ! -e "$output/chromosomes" && ! -e "$output/chromosome_list.txt" ]] || fail "output directory already contains pipeline results; choose a fresh --output-dir"
    mkdir -p -- "$output/assembly"
    fasta="$output/assembly/$(basename -- "$assembly")"
    [[ ! "$assembly" -ef "$fasta" ]] || fail "input assembly must be outside output-dir/assembly"
    list="$output/chromosome_list.txt"
    python3 "$script_dir/get_chrlist.py" "$report" > "$list"
    [[ -s "$list" ]] || fail "assembly report contains no selected chromosomes"
    awk -F '\t' '
        { if ($1 == "." || $1 == ".." || index($1, "/") || seen[$1]++ || seen[$1 "_mid"]++) exit 1 }
    ' "$list" || fail "chromosome names cause unsafe or overlapping output paths"
    mkdir -- "$output/chromosomes"
    gzip -dc "$assembly" | bgzip -c > "$fasta"
    samtools faidx "$fasta"
    if [[ "$mode" != "--create-mid-only" ]]; then
        while IFS=$'\t' read -r name accession length; do
            samtools faidx "$fasta" "$accession" | gzip -c > "$output/chromosomes/$name.fa.gz"
        done < "$list"
    fi
    if [[ "$mode" != "full" ]]; then
        while IFS=$'\t' read -r name accession length; do
            region=$(mid_region "$accession" "$length")
            samtools faidx "$fasta" "$region" | gzip -c > "$output/chromosomes/${name}_mid.fa.gz"
        done < "$list"
        input_prefix="mid_"
        output_suffix="_mid"
    fi
    for path in "$output"/chromosomes/*fa.gz; do
        run_kmers count "$path" "$input_prefix" "$output_suffix" "$cores"
    done
    cd -- "$output"
    # Each stage waits for all k values before starting its dependent stage.
    run_kmers summary "" "$input_prefix" "$output_suffix" "$cores"
    run_kmers analyze "" "$input_prefix" "$output_suffix" "$cores"
    Rscript "$script_dir/oligo_coordinates.R" "6mer${output_suffix}_relfreq.txt"
    Rscript "$script_dir/top_oligos.R" "6mer${output_suffix}_relfreq.coord.txt"
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    if [[ "${1:-}" == --internal-worker ]]; then
        shift
        [[ $# == 5 ]] || fail "invalid internal worker arguments"
        kmer_worker "$@"
    else
        main "$@"
    fi
fi
