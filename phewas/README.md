# PheWAS in All of Us for pathogenic HMBS variant carriers

**Author:** Shahil Pema (shahil.pema@bcm.edu · GitHub [@ShahilPema](https://github.com/ShahilPema)). Direct code-level questions about this notebook here.

## Purpose

Phenome-wide association study (PheWAS) testing whether carriers of
pathogenic HMBS variants — the gene whose loss-of-function causes acute
intermittent porphyria — show enrichment for any of hundreds of EHR-derived
phenotypes in the All of Us Research Program (v8 release, n = 317,969).

The notebook handles variant identification and annotation, mapping ICD codes
to phecodes (PhecodeX), Firth penalized logistic regression across phenotypes
(via an R helper script), FDR correction, and figure generation (forest plot
and summary statistics table).

## Where this can run

**This notebook only runs inside the All of Us Researcher Workbench.** All of
Us individual-level data cannot leave the secure environment, so the notebook
cannot be executed locally. It is published here for transparency and to allow
researchers with their own Workbench access to reproduce the analysis.

To request Workbench access, see <https://www.researchallofus.org/>.

## Source workspace

The reference run lives in the All of Us Researcher Workbench workspace
titled **"Identifying Phenotypes Associated with Pathogenic HMBS
Variants"**, owned by the Zoghbi lab. Other Workbench users with Controlled
Tier access can request to be added as collaborators on that workspace, or
can copy/duplicate the workspace into their own billing project to obtain
a self-owned replica. Either path produces a runnable copy of the analysis
without leaving the secure environment.

## How to reproduce (with Workbench access)

1. Log in to the All of Us Researcher Workbench at
   <https://workbench.researchallofus.org/>.
2. Either duplicate the source workspace **"Identifying Phenotypes
   Associated with Pathogenic HMBS Variants"** (preferred — preserves the
   exact dataset configuration the manuscript was run against), or create
   a new workspace against the **Controlled Tier, v8** dataset (the
   analysis depends on short-read WGS variant calls and the curated EHR
   domain, both of which are Controlled Tier resources).
3. Start a Jupyter environment using the **Hail Genomic Analysis** image.
   This image ships pre-installed with every Python and R dependency the
   notebook needs (see `requirements.txt` for the list).
4. Upload `AoU_HMBS_PheWAS.ipynb` into the workspace via the Jupyter
   file browser, or `gsutil cp` it into the workspace bucket and pull
   it into the running notebook environment. (If you duplicated the
   source workspace, the notebook is already there.)
5. Open the notebook and execute the cells **in order, top to bottom**.
   Several cells download supporting files from public sources
   (e.g. PhecodeX label tables from <https://github.com/PheWAS/PhecodeX>)
   and stage them in the workspace bucket, so the workspace must have
   outbound network access enabled in the security policy.
6. The notebook shells out to an R script (`Rscript RunFirths.R`) for the
   Firth penalized logistic regression step. That script is **not bundled
   in this directory** — it is staged into the working directory by an
   earlier notebook cell (or supplied separately by the analyst). If you
   are reproducing the analysis from this repository, place a copy of the
   `RunFirths.R` helper alongside the notebook before running the
   regression cell, or add the cell that materializes it.

## ClinVar snapshot

The notebook fetches **`clinvar_20251013.vcf.gz`** (NCBI ClinVar release
dated 2025-10-13) from the public NCBI archive at
`https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/` and uses it to
identify pathogenic / likely pathogenic HMBS variants. This is the same
ClinVar release used by the SCHEMA/BipEx collapsing analysis (see
`collapsing-analysis/inputs/annotation_manifest.json` for the matching MD5).
If you re-run the notebook today, the wget cell will pull whichever release
is currently on the NCBI weekly mirror — pin the URL to
`vcf_GRCh38/archive_2.0/2025/clinvar_20251013.vcf.gz` if you need to
reproduce the manuscript exactly.

## Expected outputs

- A summary statistics table of phecode-level associations, including
  odds ratios, 95% confidence intervals, raw and FDR-adjusted p-values,
  and case/control counts for each phecode.
- A forest plot summarizing effect sizes across phecodes for the
  manuscript figure.

Intermediate artifacts (covariate tables, phecode case tables, phecode
metadata) are written to the workspace bucket under `tmp/`.

## Dependencies

The notebook imports:

- `google.cloud.bigquery` — pulls cohort and EHR tables from the Curated Data
  Repository.
- `pyspark` (with `pyspark.sql`) — Spark session used by Hail and for joins
  against AoU's denormalized EHR tables.
- `hail` — variant filtering and annotation against the AoU short-read WGS
  call set.
- `pyarrow.parquet` — local I/O for staged tables (e.g. the AoU variant
  whitelist parquet).
- `pandas`, `numpy` — tabular wrangling and downstream summarization.
- `matplotlib` (`pyplot`, `patches`) — forest plot and figure rendering.
- `scipy.stats.fisher_exact` — sanity-check 2×2 tests alongside the Firth
  regression.
- An R step invoked as `Rscript RunFirths.R`, which requires the `logistf` R
  package (Firth penalized logistic regression).

All Python packages and R + `logistf` are pre-installed on the Hail Genomic
Analysis Workbench image; nothing needs to be installed manually under the
default configuration. See `requirements.txt` for the full list.
