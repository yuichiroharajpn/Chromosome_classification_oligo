# Chromosome oligonucleotide classification

Classify chromosomes by k-mer composition and assign Larger/Smaller classes using cross-species macrosynteny. Run the commands below from `chrom_oligo_classification/`.

## Environment

### Requirements

- Python ≥3.9; pandas and NumPy
- R 4.3.0; ape, phangorn, dplyr, tibble, ca, readr and stringr
- **`jellyfish=2.2.10`**
- samtools and HTSlib (provides `bgzip`)
- GNU parallel for parallel k-mer analysis
- **macrosyntR 0.2.19**; remotes for installation

### Installation

Create the conda environment using the YAML for your OS:

```bash
# Linux
conda env create -f environment.linux.yml

# macOS: use this instead
conda env create -f environment.macos.yml
```

Activate it and install macrosyntR:

```bash
conda activate chrom-oligo-classification
Rscript install_macrosyntR.R
```

The installer uses `remotes::install_version("macrosyntR", version = "0.2.19")`. Other dependencies are specified in the YAML files. To update an existing environment, replace `conda env create` with `conda env update`.

## 01: K-mer classification

Count strand-specific 2–9mers and perform chromosome clustering and correspondence analysis. Inputs are a genomic FASTA (`.fna.gz`) and its NCBI assembly report.

### Example: zebra finch

Download the assembly and report into `01_kmer_classification/ref_assembly/`:

```bash
bash 01_kmer_classification/ref_assembly/download_assembly.sh
```

Analyze the central regions of the chromosomes:

```bash
bash 01_kmer_classification/run_01_kmer_classification.sh \
  01_kmer_classification/ref_assembly/GCF_003957565.2_bTaeGut1.4.pri_genomic.fna.gz \
  01_kmer_classification/ref_assembly/GCF_003957565.2_bTaeGut1.4.pri_assembly_report.txt \
  --create-mid-only \
  --cores 8 \
  --output-dir 01_kmer_classification/kmer_results
```

| Option | Extracted sequences | Analysis target |
| --- | --- | --- |
| No region option | Full chromosomes | Full chromosomes |
| `--create-mid` | Full chromosomes and central 5–95% regions | Central regions |
| `--create-mid-only` | Central 5–95% regions only | Central regions |

`--cores` accepts 1–8 (default: 1) and controls concurrent k-mer jobs through GNU parallel. Use a fresh output directory for each run.

### Outputs

Results are written to `01_kmer_classification/kmer_results/` in this example:

- `chromosome_list.txt`, `assembly/` and `chromosomes/`: chromosome list, indexed assembly, extracted sequences and k-mer counts.
- `Kmer_mid_freq.txt`, `Kmer_mid_relfreq.txt`: count and relative-frequency tables (K = 2–9).
- `Kmer_mid_relfreq.*`: dendrograms, correspondence-analysis results and chromosome groups; `*.dendro.grp.h.txt` files are used in step 02.
- `6mer_mid_relfreq.coord*` and `6mer_mid_relfreq.contr.txt`: oligo coordinates, top oligos and eigenvalue summaries.

Full-chromosome analysis uses output names without `_mid`.

## 02: Cross-species assignment

Assign Larger/Smaller classes using reference chromosomes and ancestral linkage groups. Each comparison selects the k ≥3 and `height_class` combination with the highest unrounded AutoXZ Informedness. Ties favor smaller k, then the existing height-column order.

### Example: zebra finch versus chicken

The example files are included under `02_cross_species_assignment/`:

| File | Description |
| --- | --- |
| `inparanoid_out/One2OnePairs` | InParanoid one-to-one pairs from the canonical proteins below; column 1: zebra finch, column 2: chicken. |
| `ref_chick/Huang_ancestral.txt` | Chicken gene–ancestral linkage group correspondence from Huang 2023. |
| `ref_chick/chicken.v23.pep.canonical.bed` | GGswu chicken annotation from Huang 2023, restricted to proteins from canonical transcripts. |
| `ref_chick/chicken.v23.pep.canonical.fa.gz` | Corresponding chicken protein sequences. |
| `ref_zebrafinch/GCF_003957565.2_bTaeGut1.4.pri_protein.canonical.bed` | Zebra finch annotation, restricted to proteins from canonical transcripts. |
| `ref_zebrafinch/GCF_003957565.2_bTaeGut1.4.pri_protein.canonical.fa.gz` | Corresponding zebra finch protein sequences. |
| `ref_chick/chromosome_list.class.txt` | Chicken chromosome Larger/Smaller classes. |
| `ref_chick/Huang_ancestral.class.txt` | Ancestral linkage group Larger/Smaller classes. |

BED files have four columns: chromosome, start, end and protein ID. Pair and class tables have two columns without headers. Protein FASTAs document the InParanoid inputs; this workflow starts from `One2OnePairs` and does not rerun InParanoid.

Edit `02_cross_species_assignment/config.cfg` for your dataset. The example uses the bundled inputs and connects the two workflows:

```ini
kmer_dir = ../01_kmer_classification/kmer_results
output = cross_species_results
```

Paths are relative to the config file, or may be absolute. Write `key = path` without quotes; leave optional entries empty.

```bash
bash 02_cross_species_assignment/run_02_cross_species_assignment.sh \
  --config 02_cross_species_assignment/config.cfg
```

Use an empty output directory. Command-line path options override config values; `--help` lists the available options.

### Outputs

Results are written to `02_cross_species_assignment/cross_species_results/`:

- `inparanoid_out/`: ortholog pairs and macrosynteny tables.
- `stats_mid_relfreq.dendro.grp.refgalGalT2T.h.txt`: reference-chromosome scores.
- `stats_mid_relfreq.dendro.grp.refanc_jv.h.rev.txt`: ancestral-linkage-group scores.
- `match_mid_relfreq.dendro.grp.h.best.txt`: chromosome classes at the reference optimum.
- `match_mid_relfreq.dendro.grp.refanc_jv.h.best.refgalGalT2T.txt`: ancestral matches at the ancestral optimum.
- `reference.best.tsv`, `ancestral.best.tsv`: selected k and height class.
