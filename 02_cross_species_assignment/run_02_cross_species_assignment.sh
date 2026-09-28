#!/usr/bin/env bash
# Pair orientation: supplied One2OnePairs = target (zebra finch), reference (chicken).
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
fail() { echo "Error: $*" >&2; exit 1; }
usage() {
    cat <<'HELP'
Usage: run_02_cross_species_assignment.sh [--config FILE] [options]
With no arguments, load the bundled config.cfg. Command-line paths override config values.
  --config FILE             Plain key = path config (paths relative to config)
Required (via config or command line):
  --kmer-dir DIR             Step 01 *mer_mid_relfreq.dendro.grp.h.txt files
  --pairs FILE               One2OnePairs: target ID, reference ID; no header
  --target-bed FILE          Target canonical protein BED (four columns)
  --reference-bed FILE       Reference canonical protein BED (four columns)
  --ancestral-table FILE     Huang_ancestral.txt
  --reference-classes FILE   Reference chromosome, Larger/Smaller; no header
  --ancestral-classes FILE   Ancestral linkage group, Larger/Smaller; no header
  --output-dir DIR           New/empty output directory
Ancestral-synteny source (default: bundled macrosynteny_ancestral.R):
  --ancestral-script FILE    Override bundled ancestral R script (reference = sp1)
  --ancestral-synteny FILE   Existing macrosynteny.anc_jv.rev.txt (ancestor = sp1)
Optional:
  --reference-synteny FILE   Existing macrosynteny.txt (target = sp1)
  -h, --help
HELP
}
absolute_file() { [[ -f "$1" ]] || fail "file not found: $1"; printf '%s/%s\n' "$(cd -- "$(dirname -- "$1")" && pwd -P)" "$(basename -- "$1")"; }
main() {
    local kmer_dir="" pairs="" target_bed="" reference_bed="" ancestral_table=""
    local reference_classes="" ancestral_classes="" output="" ancestral_script=""
    local ancestral_synteny="" reference_synteny="" value k height
    local config_file="" config_dir="" index line key seen="|" line_number=0
    if [[ $# -eq 0 ]]; then set -- --config "$script_dir/config.cfg"; fi
    local -a arguments=("$@")
    # Load config first so explicit command-line options always take precedence.
    for ((index=0; index<${#arguments[@]}; index++)); do
        if [[ "${arguments[index]}" == "--config" ]]; then
            ((index + 1 < ${#arguments[@]})) || fail "--config requires a file"
            [[ -z "$config_file" ]] || fail "specify --config only once"
            config_file=${arguments[index+1]}
        fi
    done
    if [[ -n "$config_file" ]]; then
        config_file=$(absolute_file "$config_file")
        config_dir=$(dirname -- "$config_file")
        # Parse data only; never source/evaluate config content as shell commands.
        while IFS= read -r line || [[ -n "$line" ]]; do
            line_number=$((line_number + 1))
            line=${line%$'\r'}
            line="${line#"${line%%[![:space:]]*}"}"
            line="${line%"${line##*[![:space:]]}"}"
            [[ -n "$line" && "$line" != \#* ]] || continue
            [[ "$line" == *=* ]] || fail "$config_file:$line_number: expected key = path"
            key=${line%%=*}
            key="${key%"${key##*[![:space:]]}"}"
            value=${line#*=}
            value="${value#"${value%%[![:space:]]*}"}"
            case "$key" in
                kmer_dir|pairs|target_bed|reference_bed|ancestral_table|reference_classes|ancestral_classes|output|ancestral_script|ancestral_synteny|reference_synteny) ;;
                *) fail "$config_file:$line_number: unknown key: $key" ;;
            esac
            [[ "$seen" != *"|$key|"* ]] || fail "$config_file:$line_number: duplicate key: $key"
            seen="$seen$key|"
            if [[ -n "$value" && "$value" != /* ]]; then value="$config_dir/$value"; fi
            printf -v "$key" '%s' "$value"
        done < "$config_file"
    fi
    while [[ $# -gt 0 ]]; do
        [[ "$1" != -h && "$1" != --help ]] || { usage; return; }
        [[ $# -ge 2 ]] || fail "missing value for $1"
        value=$2
        case "$1" in
            --config) ;;
            --kmer-dir) kmer_dir=$value ;;
            --pairs) pairs=$value ;;
            --target-bed) target_bed=$value ;;
            --reference-bed) reference_bed=$value ;;
            --ancestral-table) ancestral_table=$value ;;
            --reference-classes) reference_classes=$value ;;
            --ancestral-classes) ancestral_classes=$value ;;
            --output-dir) output=$value ;;
            --ancestral-script) ancestral_script=$value ;;
            --ancestral-synteny) ancestral_synteny=$value ;;
            --reference-synteny) reference_synteny=$value ;;
            *) fail "unknown option: $1" ;;
        esac
        shift 2
    done
    [[ -n "$kmer_dir" && -n "$pairs" && -n "$target_bed" && -n "$reference_bed" && -n "$ancestral_table" && -n "$reference_classes" && -n "$ancestral_classes" && -n "$output" ]] || { usage >&2; return 1; }
    if [[ -z "$ancestral_script" && -z "$ancestral_synteny" ]]; then
        ancestral_script="$script_dir/macrosynteny_ancestral.R"
    fi
    [[ -z "$ancestral_script" || -z "$ancestral_synteny" ]] || fail "choose only one ancestral-synteny source"
    kmer_dir=$(cd -- "$kmer_dir" && pwd -P)
    compgen -G "$kmer_dir/*mer_mid_relfreq.dendro.grp.h.txt" >/dev/null || fail "no mid dendrogram files found"
    pairs=$(absolute_file "$pairs")
    target_bed=$(absolute_file "$target_bed")
    reference_bed=$(absolute_file "$reference_bed")
    ancestral_table=$(absolute_file "$ancestral_table")
    reference_classes=$(absolute_file "$reference_classes")
    ancestral_classes=$(absolute_file "$ancestral_classes")
    if [[ -n "$ancestral_script" ]]; then ancestral_script=$(absolute_file "$ancestral_script"); fi
    if [[ -n "$ancestral_synteny" ]]; then ancestral_synteny=$(absolute_file "$ancestral_synteny"); fi
    if [[ -n "$reference_synteny" ]]; then reference_synteny=$(absolute_file "$reference_synteny"); fi
    python3 -c 'import pandas, numpy'
    if [[ -z "$reference_synteny" || -n "$ancestral_script" ]]; then
        Rscript -e 'for (p in c("dplyr", "macrosyntR")) if (!requireNamespace(p, quietly=TRUE)) stop(paste("Missing R package:", p)); if (!all(c("sp1_bed", "sp2_bed") %in% names(formals(macrosyntR::load_orthologs)))) stop("macrosyntR API differs from legacy code: restore the original compatible version")'
    fi
    python3 "$script_dir/validate_inputs.py" "$pairs" "$target_bed" "$reference_bed" "$reference_classes" "$ancestral_classes"
    mkdir -p -- "$output"
    output=$(cd -- "$output" && pwd -P)
    [[ -z "$(ls -A -- "$output")" ]] || fail "output directory must be empty"
    mkdir -- "$output/inparanoid_out"
    cp -- "$pairs" "$output/inparanoid_out/One2OnePairs"
    awk '{print $2"\t"$1}' "$pairs" > "$output/inparanoid_out/One2OnePairs.rev"
    if [[ -n "$reference_synteny" ]]; then
        cp -- "$reference_synteny" "$output/inparanoid_out/macrosynteny.txt"
    else
        Rscript "$script_dir/macrosynteny.R" "$output/inparanoid_out/One2OnePairs" "$target_bed" "$reference_bed" > "$output/inparanoid_out/macrosynteny.txt"
    fi
    if [[ -n "$ancestral_synteny" ]]; then
        cp -- "$ancestral_synteny" "$output/inparanoid_out/macrosynteny.anc_jv.rev.txt"
    else
        Rscript "$ancestral_script" "$output/inparanoid_out/One2OnePairs.rev" "$reference_bed" "$target_bed" "$ancestral_table" > "$output/inparanoid_out/macrosynteny.anc_jv.rev.txt"
    fi
    python3 "$script_dir/rank_reference.py" "$kmer_dir" "$output/inparanoid_out/macrosynteny.txt" "$reference_classes" --selection-file "$output/reference.best.tsv" > "$output/stats_mid_relfreq.dendro.grp.refgalGalT2T.h.txt"
    python3 "$script_dir/rank_ancestral.py" "$kmer_dir" "$output/inparanoid_out/macrosynteny.anc_jv.rev.txt" "$ancestral_classes" --selection-file "$output/ancestral.best.tsv" > "$output/stats_mid_relfreq.dendro.grp.refanc_jv.h.rev.txt"
    IFS=$'\t' read -r k height < "$output/reference.best.tsv"
    python3 "$script_dir/match_dendro.py" "$kmer_dir/${k}mer_mid_relfreq.dendro.grp.h.txt" "$height" > "$output/match_mid_relfreq.dendro.grp.h.best.txt"
    IFS=$'\t' read -r k height < "$output/ancestral.best.tsv"
    python3 "$script_dir/match_anc_chrom.py" "$kmer_dir/${k}mer_mid_relfreq.dendro.grp.h.txt" "$output/inparanoid_out/macrosynteny.anc_jv.rev.txt" "$ancestral_classes" "$height" "$ancestral_table" > "$output/match_mid_relfreq.dendro.grp.refanc_jv.h.best.refgalGalT2T.txt"
}
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then main "$@"; fi
