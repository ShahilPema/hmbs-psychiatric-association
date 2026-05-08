#!/usr/bin/env python3
"""Standalone HMBS collapsing driver.

Consumes the pre-annotated parquets produced by `annotate_inputs.py` and runs
the ancestry-stratified CMH burden test for HMBS.

Variants in the input parquets are already restricted to the AoU rare-variant
whitelist (AC <= 5 in All of Us). At run time, an additional row-level filter
drops any variant whose total case + control AC across SCHEMA / BipEx strata
exceeds `--row-ac-max` (default 5), so only rare variants in both the
population reference and the case-control cohorts contribute to the burden
test.
"""

from __future__ import annotations

import argparse
import glob
import json
import logging
import sys
from pathlib import Path

import pandas as pd
import polars as pl

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))

from hmbs_backend import run_hmbs_collapsing  # noqa: E402
from write_excel import write_excel  # noqa: E402


DEFAULT_RUN_NAME = "Prob_with_AouFilt"
DEFAULT_OUTPUT_ROOT = _HERE / "Results"
DEFAULT_INPUTS_DIR = _HERE / "inputs"
DEFAULT_COHORTS = {
    "schema": ("HMBS_SCHEMA_ONLY", DEFAULT_INPUTS_DIR / "HMBS_SCHEMA_ONLY.annotated.parquet"),
    "bipex": ("HMBS_Bipex_ONLY", DEFAULT_INPUTS_DIR / "HMBS_Bipex_ONLY.annotated.parquet"),
    "combined": ("HMBS_Combined", DEFAULT_INPUTS_DIR / "HMBS_Combined.annotated.parquet"),
}

# Mapping from backend mask column names → report names used downstream. Matches
# VARIANT_MAPPING in the original src/Validation/SCHEMA/HMBS/run_hmbs_collapsing.py.
VARIANT_MAPPING = {
    "syn_only": "syn",
    "mis_990_AM": "mis_am_gt990",
    "mis_906_AM": "mis_am_gt906",
    "lof_or_mis_greater_than_3_MPC": "lof_or_mis_MPC_gt3",
    "lof_or_mis_990_AM": "lof_or_mis_AM_gt990",
    "lof_or_mis_906_AM": "lof_or_mis_AM_gt906",
    "lof_or_Clinvar_Path": "lof_or_ClinVar_P_LP",
    "lof_all": "lof",
    "frameshift": "lof_frameshift_only",
    "Stop_gained": "lof_stop_gained_only",
    "Splice_SNV": "lof_splice_SNV_only",
    "Clinvar_Path": "ClinVar_P_LP",
}
SIGNIFICANCE_THRESHOLD = 1.4e-6


logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run-name", default=DEFAULT_RUN_NAME, help=f"Run label. Default: {DEFAULT_RUN_NAME}")
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT,
                        help=f"Root directory for outputs. Default: {DEFAULT_OUTPUT_ROOT}")
    parser.add_argument("--schema-parquet", type=str, default=str(DEFAULT_COHORTS["schema"][1]),
                        help="Annotated SCHEMA-only parquet or glob.")
    parser.add_argument("--bipex-parquet", type=str, default=str(DEFAULT_COHORTS["bipex"][1]),
                        help="Annotated BipEx-only parquet or glob.")
    parser.add_argument("--combined-parquet", type=str, default=str(DEFAULT_COHORTS["combined"][1]),
                        help="Annotated Combined parquet or glob.")
    parser.add_argument("--cohort", dest="cohorts", action="append", choices=list(DEFAULT_COHORTS),
                        help="Optional cohort subset to run. Repeat for multiple.")
    parser.add_argument("--n-parallel-mis-cols", type=int, default=1,
                        help="Number of masks to process in parallel inside the backend.")
    parser.add_argument("--clinvar-mode", choices=["plp", "hq_plp"], default="hq_plp",
                        help="ClinVar mask definition (see hmbs_backend._with_clinvar_pathogenic_flag). "
                             "Default hq_plp matches the reference HMBS workbook (SCHEMA Clinvar_Path 3/3).")
    parser.add_argument("--row-ac-max", type=int, default=5,
                        help="Max case+ctrl row-AC (summed across strata) per variant. Default 5; pass -1 to disable.")
    return parser.parse_args()


def validate_input_pattern(file_pattern: str) -> None:
    if not glob.glob(file_pattern):
        raise FileNotFoundError(f"No parquet files matched: {file_pattern}")


def add_model_column(results: pl.DataFrame) -> pl.DataFrame:
    return results.with_columns(
        pl.when(pl.col("Name").str.contains("NDD"))
          .then(pl.lit("NDD"))
          .otherwise(pl.lit("AM"))
          .alias("Model")
    )


def summarize_significant_hits(results: pl.DataFrame, model_name: str | None = None) -> pd.DataFrame:
    filtered = results.filter(pl.col("exact_cmh_pvalue") < SIGNIFICANCE_THRESHOLD)
    if model_name is not None:
        filtered = filtered.filter(pl.col("Name").str.contains(model_name))

    if model_name is None:
        summary = filtered.group_by("Name").agg(
            pl.col("Gene_Name").n_unique().alias("count"),
            pl.col("Gene_Name").unique().sort().alias("Genes"),
        )
    else:
        summary = filtered.group_by("Gene_Name").agg(
            pl.col("Name").n_unique().alias("count"),
            pl.col("Name").unique().sort().alias("Name"),
        )
    return summary.sort("count", descending=True).to_pandas()


def build_excel_payload(results: pl.DataFrame) -> dict[str, pd.DataFrame]:
    results_with_model = add_model_column(results)
    return {
        "OR_pass_No_Tier": results_with_model.to_pandas(),
        "Gene_sum": summarize_significant_hits(results_with_model),
        "Model_sum_NDD": summarize_significant_hits(results_with_model, model_name="NDD"),
        "Model_sum_AM": summarize_significant_hits(results_with_model, model_name="AM"),
    }


def remap_variant_names(results: pl.DataFrame, run_label: str) -> pl.DataFrame:
    run_name = run_label.replace("HMBS_", "").replace("_ONLY", "")
    remapped = (
        results.with_columns(
            pl.col("Name").map_elements(lambda v: VARIANT_MAPPING.get(v), return_dtype=pl.String).alias("Name")
        )
        .drop_nulls("Name")
        .with_columns(pl.lit(run_name).alias("Run"))
    )
    return remapped.select(["Run"] + [c for c in remapped.columns if c != "Run"])


def run_single_cohort(
    run_label: str,
    file_pattern: str,
    output_dir: Path,
    run_name: str,
    n_parallel_mis_cols: int,
    clinvar_mode: str,
    row_ac_max: int | None,
) -> pl.DataFrame | None:
    logger.info(f"Running {run_label} from {file_pattern} (clinvar_mode={clinvar_mode}, row_ac_max={row_ac_max})")
    validate_input_pattern(file_pattern)

    out_cmh, _ = run_hmbs_collapsing(
        file_pattern,
        n_parallel_mis_cols=n_parallel_mis_cols,
        clinvar_mode=clinvar_mode,
        row_ac_max=row_ac_max,
    )

    if "OR_pass_No_Tier" not in out_cmh:
        logger.warning(f"No OR_pass_No_Tier output for {run_label}")
        return None

    raw_results = out_cmh["OR_pass_No_Tier"]
    excel_payload = build_excel_payload(raw_results)

    excel_path = output_dir / f"{run_label}_{run_name}_merged.xlsx"
    parquet_path = output_dir / f"{run_label}_{run_name}_merged.parquet"
    logger.info(f"Writing {excel_path}")
    write_excel(excel_payload, str(excel_path))
    raw_results.write_parquet(parquet_path)

    return remap_variant_names(raw_results, run_label)


def log_manifest_if_present(annotated_parquets: list[str]) -> None:
    """If the annotation manifest sits next to the annotated parquets, log key fields."""
    for path_str in annotated_parquets:
        parent = Path(path_str).expanduser().resolve().parent
        manifest_file = parent / "annotation_manifest.json"
        if manifest_file.exists():
            try:
                info = json.loads(manifest_file.read_text())
                logger.info(
                    "Annotation manifest: generated_at=%s, clinvar_vcf=%s (md5=%s, %d HMBS records)",
                    info.get("generated_at"),
                    info.get("clinvar_vcf", {}).get("path"),
                    info.get("clinvar_vcf", {}).get("md5"),
                    info.get("clinvar_vcf", {}).get("hmbs_records"),
                )
            except Exception as exc:
                logger.warning(f"Could not read manifest {manifest_file}: {exc}")
            return


def main() -> None:
    args = parse_args()
    row_ac_max: int | None = None if args.row_ac_max is not None and args.row_ac_max < 0 else args.row_ac_max

    output_dir = args.output_root / f"HMBS_{args.run_name}"
    output_dir.mkdir(parents=True, exist_ok=True)

    cohort_paths = {
        "schema": (DEFAULT_COHORTS["schema"][0], args.schema_parquet),
        "bipex": (DEFAULT_COHORTS["bipex"][0], args.bipex_parquet),
        "combined": (DEFAULT_COHORTS["combined"][0], args.combined_parquet),
    }
    selected = args.cohorts or list(DEFAULT_COHORTS)

    log_manifest_if_present([cohort_paths[c][1] for c in selected])

    merged_runs: list[pl.DataFrame] = []
    for cohort in selected:
        run_label, file_pattern = cohort_paths[cohort]
        remapped = run_single_cohort(
            run_label=run_label,
            file_pattern=file_pattern,
            output_dir=output_dir,
            run_name=args.run_name,
            n_parallel_mis_cols=args.n_parallel_mis_cols,
            clinvar_mode=args.clinvar_mode,
            row_ac_max=row_ac_max,
        )
        if remapped is not None and remapped.height > 0:
            merged_runs.append(remapped)

    if not merged_runs:
        logger.warning("No remapped HMBS runs were generated; skipping merged workbook")
        return

    merged = pl.concat(merged_runs, how="diagonal_relaxed").to_pandas()
    merged_excel = output_dir / f"HMBS_All_Runs_{args.run_name}_merged.xlsx"
    logger.info(f"Writing merged workbook {merged_excel}")
    write_excel({"All Runs": merged}, str(merged_excel))


if __name__ == "__main__":
    main()
