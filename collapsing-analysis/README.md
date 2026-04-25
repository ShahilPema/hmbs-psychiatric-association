# Rare-variant burden analysis in SCHEMA and BipEx for HMBS

**Author:** Blake A. Weido (Blake.Weido@bcm.edu · GitHub [@bweido](https://github.com/bweido)). Direct code-level questions about this pipeline here.

## Purpose

Rare-variant case-control "collapsing" burden analysis for HMBS, supporting
the manuscript *A Genetic and Biochemical Reassessment of Acute Intermittent
Porphyria in Psychiatric Disorders* (Pema, Issa, Weido, Zoghbi). Tests
enrichment of qualifying variant classes (LoF, AlphaMissense tiers, MPC
tiers, ClinVar P/LP, synonymous) in the Schizophrenia Exome Meta-analysis
(SCHEMA) and Bipolar Exomes (BipEx) cohorts against gnomAD controls using
the exact Cochran-Mantel-Haenszel test, ancestry-stratified. The output of
`run_collapsing.py` is what populates Table 1 of the manuscript (odds
ratios, 95% confidence intervals, and p-values across every qualifier
mask listed below).

This directory depends only on the Python scientific stack, R + rpy2, and the
files vendored in `data/`. It is derived from
`src/Tools/HMBS_cmh_r_backend.py` and
`src/Validation/SCHEMA/HMBS/run_hmbs_collapsing.py` of the upstream
`Missense_Predictor` repo, with PEXT, multi-consequence, predictions-join, and
AoU-whitelist logic either baked into the annotated inputs or removed entirely.

## Directory layout

```
HMBS_Collapsing/
├── README.md                          — this file
├── run_collapsing.py                  — end-user driver (runs the CMH test)
├── hmbs_backend.py                    — CMH backend (polars + rpy2)
├── write_excel.py                     — vendored Excel writer
├── annotate_inputs.py                 — maintainer-only; rebuilds inputs/ from raw parquets
├── data/
│   ├── clinvar_20251013.vcf.gz        — ClinVar snapshot baked into inputs/
│   ├── hgnc_for_mapping_pc.txt        — HGNC ENSG → symbol mapping
│   └── HMBS_Combined_AoUFilt.parquet  — AoU v8 variant whitelist
└── inputs/
    ├── HMBS_SCHEMA_ONLY.annotated.parquet
    ├── HMBS_Bipex_ONLY.annotated.parquet
    ├── HMBS_Combined.annotated.parquet
    └── annotation_manifest.json       — provenance of the three parquets above
```

## Environment

Two environment specs ship with this package:

- `environment.yml` — minimal top-level spec with version pins from the reference
  rerun (good for `mamba env create`). Human-readable.
- `environment.lock.yml` — fully pinned, transitively-complete export of the
  environment that was verified to reproduce
  `HMBS_Prob_with_AouFilt_HQ_vcf_20251013_rerun_20260421` bit-for-bit. Use this
  if you need byte-identical reproduction of the reference run.

Create and activate:

```bash
# minimal
mamba env create -n hmbs_collapsing -f environment.yml

# or exact lock
mamba env create -n hmbs_collapsing -f environment.lock.yml

mamba activate hmbs_collapsing
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
export R_HOME=$CONDA_PREFIX/lib/R
```

`conda` works too — the YAMLs are plain conda env files and `conda env create -f
environment.yml` behaves identically. On conda ≥ 23.10 the default solver is
already `libmamba`, so solve time is comparable to `mamba`; on older conda
installs set `conda config --set solver libmamba` first, or just use `mamba`.

The LD_LIBRARY_PATH and R_HOME exports are required because rpy2 needs to find the
conda-managed libR. The package was verified with:

- Python 3.12.13
- polars 1.40.0, pandas 3.0.2, numpy 2.4.3, pyarrow 23.0.1
- joblib 1.5.3, rpy2 3.6, r-base 4.5.3
- xlsxwriter 3.2.9, openpyxl 3.1.5

## Typical workflow (end user)

```bash
cd HMBS_Collapsing/
$PY run_collapsing.py \
    --run-name Prob_with_AouFilt \
    --output-root Results/
```

Defaults reproduce the reference HMBS workbook (`--clinvar-mode hq_plp`,
`--row-ac-max 5`). Pass `--clinvar-mode plp` if you want the permissive
ClinVar mask.

This produces, under `Results/HMBS_Prob_with_AouFilt/`:

- `HMBS_SCHEMA_ONLY_Prob_with_AouFilt_merged.{xlsx,parquet}` — SCHEMA-only results
- `HMBS_Bipex_ONLY_Prob_with_AouFilt_merged.{xlsx,parquet}` — BipEx-only results
- `HMBS_Combined_Prob_with_AouFilt_merged.{xlsx,parquet}` — pooled results
- `HMBS_All_Runs_Prob_with_AouFilt_merged.xlsx` — merged summary workbook (one
  row per cohort × qualifier, matches the format used in the upstream pipeline)

CLI flags:

- `--run-name STR` — tag appended to output filenames and folder.
- `--clinvar-mode {plp, hq_plp}` — `plp` flags any CLNSIG containing
  `Pathogenic|Likely_pathogenic`; `hq_plp` additionally requires a
  `criteria_provided,_multiple_submitters,_no_conflicts` review status.
  Default: `hq_plp` (matches the reference workbook).
- `--row-ac-max N` — drop variants whose total case + control AC across strata
  exceeds N. Default 5; pass `-1` to disable.
- `--cohort {schema,bipex,combined}` — restrict to one cohort (repeat the flag
  to select multiple). All three run by default.

## ClinVar snapshot used in the manuscript

The ClinVar release frozen for the manuscript is **`clinvar_20251013.vcf.gz`
(release date 2025-10-13)**, downloaded from the NCBI ClinVar GRCh38 VCF
archive. The exact provenance is recorded in
`inputs/annotation_manifest.json`:

| Field         | Value                                              |
| ------------- | -------------------------------------------------- |
| File          | `data/clinvar_20251013.vcf.gz`                     |
| Release date  | 2025-10-13                                         |
| MD5           | `5524de2068e04e1441e68530baf3efb7`                 |
| Size          | 174,939,681 bytes                                  |
| HMBS records  | 688                                                |
| Annotated on  | 2026-04-21                                         |

`data/clinvar_*.vcf.gz` is **not** tracked in this Git repository — the file
exceeds GitHub's per-file size limit and is excluded by `.gitignore`. You do
**not** need it to run `run_collapsing.py`: the relevant ClinVar fields
(`CLNSIG`, `CLNREVSTAT`) are already baked into `inputs/*.annotated.parquet`
at annotation time, so the manuscript's ClinVar snapshot is preserved
inside those parquets even though the raw VCF is not in the repo.

You only need to fetch a ClinVar VCF if you want to **re-annotate** the
inputs against a different snapshot. See the next section.

## Bumping the ClinVar snapshot (maintainer)

The annotated parquets in `inputs/` have the ClinVar VCF and AoU whitelist
baked in. `AlphaMissense.am_pathogenicity` is read straight from the raw cohort
parquet (verified identical to the predictions checkpoint on all HMBS
variants where both are non-null). When a newer ClinVar release comes out:

1. Drop the new `clinvar_YYYYMMDD.vcf.gz` into `data/`. The NCBI archive for a
   specific date lives at
   `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/YYYY/clinvar_YYYYMMDD.vcf.gz`.
2. Re-run the annotation step:

   ```bash
   cd HMBS_Collapsing/
   $PY annotate_inputs.py \
       --schema-parquet <path to raw SCHEMA-only HMBS parquet> \
       --bipex-parquet <path to raw BipEx-only HMBS parquet> \
       --combined-parquet <path to raw Combined HMBS parquet> \
       --clinvar-vcf data/clinvar_YYYYMMDD.vcf.gz \
       --aou-filter data/HMBS_Combined_AoUFilt.parquet \
       --output-dir inputs/
   ```

3. The script overwrites `inputs/*.annotated.parquet` and updates
   `inputs/annotation_manifest.json` with the new ClinVar MD5 and timestamp.
4. Re-run `run_collapsing.py` as usual. Output will reflect the new ClinVar
   snapshot without any code change.

## Rebuilding the raw input parquets from scratch (reference)

The three `HMBS_SCHEMA_ONLY.parquet`, `HMBS_Bipex_ONLY.parquet`, and
`HMBS_Combined.parquet` inputs come from upstream case-control summary VCFs.
For the manuscript, both the SCHEMA and BipEx VCFs were downloaded from the
SCHEMA downloads page:

- Downloads page (canonical source for both cohorts): <https://schema.broadinstitute.org/downloads>

The same files are also mirrored on the ATGU exome-browser S3 bucket and
can be fetched directly:

- SCHEMA: <https://atgu-exome-browser-data.s3.amazonaws.com/SCHEMA/SCHEMA_variant_results.vcf.bgz>
- BipEx: <https://atgu-exome-browser-data.s3.amazonaws.com/BipEx/BipEx_variant_results.vcf.bgz>

To regenerate from scratch you would need to:

1. Download the VCFs and tabix them.
2. Lift over to GRCh38 (`CrossMap` or `Picard LiftoverVcf` with the GRCh37 →
   GRCh38 chain file).
3. Annotate with VEP (MANE Select transcripts only) and with dbNSFP +
   AlphaMissense for missense tiers.
4. Join each variant row with gnomAD ancestry-stratified counts (`info.ac_case`,
   `info.ac_ctrl`, `info.an_case`, `info.an_ctrl`, `info.group` where `group` is
   one of `AFR (exomes/genomes)`, `AMR (exomes)`, `EAS (exomes)`,
   `EUR (exomes)`, `EUR-N (exomes)`, `FIN (exomes/genomes)`, `SAS (exomes)`,
   `EST (genomes)`, `Bipolar Disorder`).
5. Restrict to HMBS (`info.gene_id == ENSG00000256269`) and save as a parquet.

There is no turn-key script in this repo for steps 1-5 — the parquets were
produced by a separate upstream pipeline in the Zoghbi lab. Ask the maintainer
for the current canonical paths if you need to rebuild.

## Column reference

Columns consumed by the backend on each annotated parquet:

| Column                              | Source           | Used for                                           |
| ----------------------------------- | ---------------- | -------------------------------------------------- |
| `ID_38`                             | input            | join keys, dedup, recurrent-variant tracking       |
| `info.gene_id`                      | input            | HMBS group-by                                       |
| `info.group`                        | input            | ancestry / study stratum for CMH                    |
| `info.ac_case`, `info.ac_ctrl`      | input            | stratum-level case/control allele counts            |
| `info.an_case`, `info.an_ctrl`      | input            | stratum-level case/control allele numbers           |
| `info.consequence`                  | input            | qualifier mask construction (missense/syn/splice)   |
| `mane_select.consequence_terms`     | input            | kept for backward compat, not used in logic         |
| `mane_select.transcript_id`         | input            | predictions join key (bake-in only)                 |
| `alleles`                           | input            | source for `SNV` at annotation time                 |
| `clinvar_vcf.CLNSIG`                | annotation step  | ClinVar P/LP flag                                   |
| `clinvar_vcf.CLNREVSTAT`            | annotation step  | HQ review-status filter in `hq_plp` mode            |
| `AlphaMissense.am_pathogenicity`    | input            | AM missense tiers (0.906 / 0.990 / …)               |
| `dbNSFP.MPC_score`                  | input            | MPC tiers                                           |
| `SpliceAI_max`                      | input            | SpliceAI-gated syn/mis masks                        |
| `SNV`                               | annotation step  | SNV vs indel gating for splice/stop masks           |

## Qualifier masks

Defined in `hmbs_backend.create_missense_columns`. The renames applied in
`run_collapsing.py` (via `VARIANT_MAPPING`) expose these names in the merged
output:

| Mask (backend)                        | Report name            |
| ------------------------------------- | ---------------------- |
| `syn_only`                            | `syn`                  |
| `mis_906_AM`                          | `mis_am_gt906`         |
| `mis_990_AM`                          | `mis_am_gt990`         |
| `lof_or_mis_greater_than_3_MPC`       | `lof_or_mis_MPC_gt3`   |
| `lof_or_mis_990_AM`                   | `lof_or_mis_AM_gt990`  |
| `lof_or_mis_906_AM`                   | `lof_or_mis_AM_gt906`  |
| `lof_or_Clinvar_Path`                 | `lof_or_ClinVar_P_LP`  |
| `lof_all`                             | `lof`                  |
| `frameshift`                          | `lof_frameshift_only`  |
| `Stop_gained`                         | `lof_stop_gained_only` |
| `Splice_SNV`                          | `lof_splice_SNV_only`  |
| `Clinvar_Path`                        | `ClinVar_P_LP`         |

Mask details:

- `mis_###_AM` fires on missense variants with `AlphaMissense.am_pathogenicity`
  strictly greater than the listed threshold. 0.906 and 0.990 correspond to
  ACMG Moderate and Strong evidence thresholds respectively.
- LoF atoms (`Splice_SNV`, `Stop_gained`, `frameshift`) require SNV/indel gating
  as appropriate. Splice/stop indel masks (`Splice_nonSNV`, `Stop_nonSNV`) were
  dropped because no such variants exist in the HMBS input parquets; if this
  changes, revive those masks in the backend.
- `ClinVar_P_LP` has no consequence filter — it fires on any variant whose
  CLNSIG contains `Pathogenic|Likely_pathogenic` (and, in `hq_plp` mode, whose
  review status is `criteria_provided,_multiple_submitters,_no_conflicts`).

## Verification

To reproduce the reference rerun
(`Results/Validation/Collapsing/HMBS_Prob_with_AouFilt_HQ_vcf_20251013_rerun_20260421/`
in the upstream repo):

```bash
cd HMBS_Collapsing/

$PY annotate_inputs.py \
    --schema-parquet /storage/zoghbi/home/u241617/Missense_Model/Missense_Predictor/Data/HMBS/HMBS_SCHEMA_ONLY.parquet \
    --bipex-parquet /storage/zoghbi/home/u241617/Missense_Model/Missense_Predictor/Data/HMBS/HMBS_Bipex_ONLY.parquet \
    --combined-parquet /storage/zoghbi/home/u241617/Missense_Model/Missense_Predictor/Data/HMBS/HMBS_Combined.parquet \
    --clinvar-vcf data/clinvar_20251013.vcf.gz \
    --aou-filter data/HMBS_Combined_AoUFilt.parquet \
    --output-dir inputs/

$PY run_collapsing.py \
    --run-name standalone_verification \
    --clinvar-mode hq_plp \
    --row-ac-max 5 \
    --output-root Results/
```

Note the `--clinvar-mode hq_plp`: the reference rerun directory name
`HMBS_Prob_with_AouFilt_HQ_vcf_20251013_rerun_20260421` encodes HQ ClinVar mode
(CLNSIG P/LP **and** CLNREVSTAT `criteria_provided,_multiple_submitters,_no_conflicts`).

The resulting
`Results/HMBS_standalone_verification/HMBS_All_Runs_standalone_verification_merged.xlsx`
matches
`Results/Validation/Collapsing/HMBS_Prob_with_AouFilt_HQ_vcf_20251013_rerun_20260421/HMBS_All_Runs_Prob_with_AouFilt_HQ_vcf_20251013_rerun_20260421_merged.xlsx`
row-for-row on every numeric column. The three list-valued columns
(`fisher_p_values`, `recurrent_variants`, `strata_specific_ors`) report the same
multisets in different orders — verified by parsing and sorting the contents
before comparison.
