# A genetic and biochemical reassessment of acute intermittent porphyria in psychiatric disorders

**Authors:** Shahil Pema\*, Dany A. Issa\*, Blake A. Weido, Anthony W. Zoghbi
(\* co-first authors)
Baylor College of Medicine.

## Summary

This repository contains the analysis code for a manuscript reassessing
whether pathogenic HMBS variants — the cause of acute intermittent porphyria
(AIP) — are enriched in psychiatric populations, as prior literature has
claimed. We test this with a phenome-wide association study in the All of Us
Research Program (v8, n = 317,969) and a rare-variant burden ("collapsing")
analysis in the SCHEMA (24,248 cases / 97,322 controls) and BipEx (13,933
cases / 14,422 controls) cohorts. We do not replicate the dramatic enrichment
reported in earlier work; effect estimates are consistently risk-increasing
but do not reach statistical significance.

## Repository structure

```
hmbs-psychiatric-association/
├── README.md                — this file
├── LICENSE                  — MIT license
├── CITATION.cff             — machine-readable citation metadata
├── .zenodo.json             — Zenodo archival metadata
├── .gitignore
├── phewas/                  — All of Us PheWAS notebook + docs
└── collapsing-analysis/     — SCHEMA + BipEx rare-variant burden code
```

The two analysis directories are independent and self-documenting:

- [`phewas/`](phewas/) — Jupyter notebook implementing the All of Us PheWAS.
  Runs only inside the All of Us Researcher Workbench. See
  [`phewas/README.md`](phewas/README.md).
- [`collapsing-analysis/`](collapsing-analysis/) — End-to-end Python + R
  pipeline implementing the rare-variant burden test in SCHEMA and BipEx
  using the exact Cochran-Mantel-Haenszel test, ancestry-stratified. See
  [`collapsing-analysis/README.md`](collapsing-analysis/README.md).

The third arm of the manuscript — a literature review of historical
biochemical screening studies — has no associated code and is not
represented in this repository.

## Code attribution

The two analysis directories were written by different members of the
author team. The manuscript itself is co-authored by all four authors
listed at the top of this README; the breakdown below applies only to
the source code in this repository.

| Component | Author | Contact |
| --- | --- | --- |
| [`phewas/`](phewas/) | Shahil Pema | shahil.pema@bcm.edu · GitHub [@ShahilPema](https://github.com/ShahilPema) |
| [`collapsing-analysis/`](collapsing-analysis/) | Blake A. Weido | Blake.Weido@bcm.edu · GitHub [@bweido](https://github.com/bweido) |

Direct code-level questions to the relevant author. Questions about the
manuscript as a whole should go to the corresponding author.

## ClinVar snapshot

Both analyses use the same ClinVar snapshot — **`clinvar_20251013.vcf.gz`
(release date 2025-10-13)**, downloaded from the NCBI ClinVar GRCh38 VCF
archive. This is the version frozen for the manuscript and is the version
recorded in `collapsing-analysis/inputs/annotation_manifest.json` (with
md5 `5524de2068e04e1441e68530baf3efb7`). The same file path is referenced
by the PheWAS notebook. Anyone re-running the pipeline against a newer
ClinVar release should expect minor differences in the ClinVar P/LP mask;
see `collapsing-analysis/README.md` for the procedure.

## Data availability

### All of Us
All of Us individual-level data cannot be redistributed and cannot leave the
secure Researcher Workbench. The reference analysis was run in a workspace
titled **"Identifying Phenotypes Associated with Pathogenic HMBS Variants"**.
To reproduce, an investigator must:

1. Apply for All of Us Researcher Workbench access at
   <https://www.researchallofus.org/>, including completion of the registered
   and controlled tier credentialing.
2. Either request collaborator access on the source workspace
   "Identifying Phenotypes Associated with Pathogenic HMBS Variants",
   or duplicate it into their own billing project, or create a fresh
   Controlled Tier (v8) workspace and copy the notebook in.

See [`phewas/README.md`](phewas/README.md) for execution details.

### SCHEMA and BipEx
The SCHEMA and BipEx case-control summary statistics used by the collapsing
analysis are released by the consortia maintained by the Analytic and
Translational Genetics Unit at MGH. Both cohorts' per-variant case/control
count VCFs were downloaded from the SCHEMA downloads page:

- Downloads (used for both SCHEMA and BipEx): <https://schema.broadinstitute.org/downloads>
- Companion variant browsers (interactive lookup, not used by this pipeline):
  - SCHEMA: <https://schema.broadinstitute.org/>
  - BipEx: <https://bipex.broadinstitute.org/>

See [`collapsing-analysis/README.md`](collapsing-analysis/README.md) for
the exact filenames and the rebuild procedure.

## How to reproduce

### PheWAS

Requires All of Us Researcher Workbench access. Copy
`phewas/AoU_HMBS_PheWAS.ipynb` into a v8 Controlled Tier workspace running the
**Hail Genomic Analysis** image and execute the cells in order. Full
instructions in [`phewas/README.md`](phewas/README.md).

### Collapsing analysis

Runs end-to-end on a workstation given the SCHEMA and BipEx variant data.
Create the conda environment, activate it, and invoke the driver:

```bash
cd collapsing-analysis/
mamba env create -n hmbs_collapsing -f environment.yml
mamba activate hmbs_collapsing
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
export R_HOME=$CONDA_PREFIX/lib/R

python run_collapsing.py --run-name my_run --output-root Results/
```

This produces SCHEMA-only, BipEx-only, and pooled odds ratio /
confidence-interval / p-value tables across every variant qualifier mask
(the contents of Table 1 in the manuscript). Full details, including how
to refresh the bundled ClinVar snapshot or rebuild the annotated input
parquets from the upstream SCHEMA/BipEx VCFs, are in
[`collapsing-analysis/README.md`](collapsing-analysis/README.md).

## Citation

If you use this code, please cite both the software (via the Zenodo DOI in
the section below) and the article. A BibTeX entry for the article will be
finalized at acceptance; a placeholder is provided here:

```bibtex
@article{pema2026hmbs,
  title   = {A Genetic and Biochemical Reassessment of Acute Intermittent
             Porphyria in Psychiatric Disorders},
  author  = {Pema, Shahil and Issa, Dany A. and Weido, Blake A. and
             Zoghbi, Anthony W.},
  journal = {British Journal of Psychiatry},
  year    = {2026},
  note    = {Pema and Issa contributed equally},
  doi     = {TBD}
}
```

The machine-readable [`CITATION.cff`](CITATION.cff) carries the same
information and will be updated with the final DOI on acceptance.

## Code archival

[![DOI](https://zenodo.org/badge/DOI/XXXXXX.svg)](https://doi.org/XXXXXX)

This repository is archived on Zenodo with a persistent DOI. The badge above
and the DOI fields in [`CITATION.cff`](CITATION.cff) and
[`.zenodo.json`](.zenodo.json) will be populated after the first GitHub
release is created and the GitHub-Zenodo integration mints a DOI. See the
"For maintainers" section below for the one-time setup.

## License

Released under the MIT License. See [`LICENSE`](LICENSE) for the full text.

## For maintainers

### One-time Zenodo setup

1. Log in to Zenodo with the lab's GitHub credentials at
   <https://zenodo.org/account/settings/github/>.
2. Flip the toggle for the `hmbs-psychiatric-association` repository to **on**.
3. After paper acceptance, create a GitHub release (e.g. `v1.0.0`).
   Zenodo will automatically archive the release and mint a DOI.
4. Update the README badge, the DOI field in `CITATION.cff`, and any
   placeholder DOI in `.zenodo.json` with the assigned DOI.

### Before the first push

- Fill in real ORCID iDs for all four authors in
  [`.zenodo.json`](.zenodo.json) (placeholders are `0000-0000-0000-0000`).
- Update the `repository-code` field in [`CITATION.cff`](CITATION.cff) to the
  actual GitHub URL once the repo is created.

### At paper acceptance

- Replace the BibTeX `doi = {TBD}` in this README with the article DOI.
- Uncomment and fill the `doi:` field in [`CITATION.cff`](CITATION.cff)
  (both at the top level for the software and inside `preferred-citation`
  for the article).
- Replace the `10.0000/XXXXX` placeholder in `.zenodo.json` under
  `related_identifiers` with the article DOI.
- Tag a GitHub release to trigger Zenodo archival; populate the badge and
  DOI fields once the Zenodo DOI is assigned.
