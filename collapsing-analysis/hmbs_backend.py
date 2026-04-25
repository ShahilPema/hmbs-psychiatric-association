"""HMBS CMH burden-test backend (standalone).

Simplified from src/Tools/HMBS_cmh_r_backend.py of the Missense_Predictor repo.

What was dropped and why (all validated to be equivalent on HMBS inputs):

- PEXT filtering. Previously every mask was gated on PEXT.exp_prop_mean > 0.1, but the
  input df set the column to pl.lit(1) immediately, so the predicate was always True.
- multi_consequence logic. Previously every splice/stop/frameshift mask was gated on
  multi_consequence == False, but the input df set the column to pl.lit(False)
  immediately, so the predicate was always True.
- Splice_nonSNV and Stop_nonSNV masks. Both fire only when SNV=False; verified on
  SCHEMA, BipEx, and Combined parquets that no splice or stop_gained variant has
  SNV=False, so they never contribute to the LoF-indel set for HMBS. This also
  means frameshift is the only LoF-indel consequence, so a separate `lof_indel`
  mask would be a duplicate of `frameshift` — it is omitted.
- Predictions parquet join. The annotation step bakes AlphaMissense.am_pathogenicity
  straight into the annotated input parquet.
- AoU variant whitelist. Applied at annotation time, not runtime.
- SNV flag computation. Pre-computed at annotation time from the alleles list.
- Dead dicts AC_case_filter / AC_ctrl_filter. Never referenced anywhere.
- Case_OR_pass_* / Control_OR_pass_* columns. Never referenced anywhere.
"""

from __future__ import annotations

import gc
import glob
import multiprocessing as mp
import os
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
import polars.selectors as cs
from joblib import Parallel, delayed


_HERE = Path(__file__).resolve().parent
_DEFAULT_HGNC = _HERE / "data" / "hgnc_for_mapping_pc.txt"


def _safe_activate_pandas2ri(pandas2ri_module):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        try:
            pandas2ri_module.activate()
        except DeprecationWarning:
            pass


def _init_and_compute_r(payload):
    if not hasattr(_init_and_compute_r, "r_initialized"):
        os.environ.setdefault("OMP_NUM_THREADS", "1")
        os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
        os.environ.setdefault("MKL_NUM_THREADS", "1")

        global ro, pandas2ri, stats
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
        from rpy2.robjects.packages import importr

        _safe_activate_pandas2ri(pandas2ri)
        stats = importr("stats")
        _init_and_compute_r.r_initialized = True

    return _compute_exact_cmh_for_gene_r(payload)


def _run_pool_r_joblib(payloads, n_workers=8):
    print(f"Processing {len(payloads)} genes with joblib ({n_workers} workers)")
    results = Parallel(
        n_jobs=n_workers,
        backend="loky",
        batch_size="auto",
        verbose=10,
    )(delayed(_init_and_compute_r)(p) for p in payloads)
    errors = [r for r in results if "error" in r and r["error"]]
    return results, errors


def _with_clinvar_pathogenic_flag(df: pl.DataFrame, clinvar_mode: str = "hq_plp") -> pl.DataFrame:
    """Derive the HMBS ClinVar flag from raw ClinVar fields.

    Prefers the clinvar_vcf.* columns baked in at annotation time; falls back to any
    clinvar_38.* columns that may still be present in the input parquet.
    """
    if clinvar_mode not in {"plp", "hq_plp"}:
        raise ValueError(f"Unsupported clinvar_mode: {clinvar_mode}")

    if "clinvar_vcf.CLNSIG" in df.columns:
        clnsig_col = "clinvar_vcf.CLNSIG"
        revstat_col = "clinvar_vcf.CLNREVSTAT"
    else:
        clnsig_col = "clinvar_38.CLNSIG"
        revstat_col = "clinvar_38.CLNREVSTAT"

    has_clnsig = clnsig_col in df.columns
    has_revstat = revstat_col in df.columns

    if clinvar_mode == "hq_plp" and not (has_clnsig and has_revstat):
        raise ValueError("HQ ClinVar mode requires raw ClinVar CLNSIG and CLNREVSTAT fields.")

    if has_clnsig:
        plp_expr = pl.col(clnsig_col).fill_null("").str.contains("Pathogenic|Likely_pathogenic")
        if clinvar_mode == "hq_plp":
            revstat_expr = (
                pl.col(revstat_col).fill_null("").str.contains("criteria_provided")
                & pl.col(revstat_col).fill_null("").str.contains("_multiple_submitters")
                & pl.col(revstat_col).fill_null("").str.contains("_no_conflicts")
            )
            return df.with_columns((plp_expr & revstat_expr).alias("Clinvar_Pathogenic"))
        return df.with_columns(plp_expr.alias("Clinvar_Pathogenic"))

    if "Clinvar_Pathogenic" in df.columns:
        return df.with_columns(pl.col("Clinvar_Pathogenic").fill_null(False).alias("Clinvar_Pathogenic"))
    return df.with_columns(pl.lit(False).alias("Clinvar_Pathogenic"))


def create_missense_columns(df: pl.DataFrame, clinvar_mode: str = "hq_plp") -> tuple[pl.DataFrame, list[str]]:
    """Build the qualifier mask columns used by the CMH test.

    Assumes the following columns exist on the annotated input parquet:
    - info.consequence (str)
    - AlphaMissense.am_pathogenicity (Float64 or castable)
    - dbNSFP.MPC_score (Float64 or castable)
    - SpliceAI_max (Float64 or castable) — pre-computed in the upstream parquet
    - SNV (Bool) — pre-computed at annotation time
    - CLNSIG/CLNREVSTAT (from the clinvar_vcf.* columns baked in at annotation time)
    """
    df = df.with_columns([
        pl.col("AlphaMissense.am_pathogenicity").cast(pl.Float64),
        pl.col("dbNSFP.MPC_score").cast(pl.Float64),
        pl.col("SpliceAI_max").cast(pl.Float64, strict=False),
    ])
    df = _with_clinvar_pathogenic_flag(df, clinvar_mode=clinvar_mode)

    df = df.with_columns([
        # AlphaMissense tiers
        (pl.col("info.consequence").str.contains("missense") & (pl.col("AlphaMissense.am_pathogenicity") > 0.906)).alias("mis_906_AM"),
        (pl.col("info.consequence").str.contains("missense") & (pl.col("AlphaMissense.am_pathogenicity") > 0.990)).alias("mis_990_AM"),
        (pl.col("info.consequence").str.contains("missense") & (pl.col("AlphaMissense.am_pathogenicity") > 0.995)).alias("mis_995_AM"),
        (pl.col("info.consequence").str.contains("missense") & (pl.col("AlphaMissense.am_pathogenicity") > 0.997)).alias("mis_997_AM"),
        (pl.col("info.consequence").str.contains("missense") & (pl.col("AlphaMissense.am_pathogenicity") > 0.999)).alias("mis_999_AM"),

        # MPC tiers
        (pl.col("info.consequence").str.contains("missense") & (pl.col("dbNSFP.MPC_score") > 1.92)).alias("mis_youden_MPC"),
        (pl.col("info.consequence").str.contains("missense") & (pl.col("dbNSFP.MPC_score") > 2)).alias("mis_greater_than_2_MPC"),
        (pl.col("info.consequence").str.contains("missense") & (pl.col("dbNSFP.MPC_score") > 3)).alias("mis_greater_than_3_MPC"),

        # Synonymous
        pl.col("info.consequence").str.contains("synonymous_variant").alias("syn_only"),

        # SpliceAI
        ((pl.col("info.consequence").str.contains("missense") | pl.col("info.consequence").str.contains("synonymous"))
         & (pl.col("SpliceAI_max") > 0.5)).alias("syn_mis_mid_SpliceAI"),
        ((pl.col("info.consequence").str.contains("missense") | pl.col("info.consequence").str.contains("synonymous"))
         & (pl.col("SpliceAI_max") > 0.9)).alias("syn_mis_high_SpliceAI"),
        (pl.col("info.consequence").str.contains("synonymous") & (pl.col("SpliceAI_max") < 0.1)).alias("syn_benign_SpliceAI"),

        # LoF atoms
        (pl.col("info.consequence").is_in(["splice_acceptor_variant", "splice_donor_variant"]) & pl.col("SNV")).alias("Splice_SNV"),
        (pl.col("info.consequence").is_in(["stop_gained"]) & pl.col("SNV")).alias("Stop_gained"),
        pl.col("info.consequence").is_in(["frameshift_variant"]).alias("frameshift"),

        # ClinVar Pathogenic as a qualifier
        (pl.col("Clinvar_Pathogenic") == True).alias("Clinvar_Path"),
    ])

    df = df.with_columns([
        (pl.col("Splice_SNV") | pl.col("Stop_gained")).alias("lof_snv"),
    ])

    df = df.with_columns((pl.col("lof_snv") | pl.col("frameshift")).alias("lof_all"))

    df = df.with_columns([
        (pl.col("lof_all") | pl.col("mis_906_AM")).alias("lof_or_mis_906_AM"),
        (pl.col("lof_all") | pl.col("mis_990_AM")).alias("lof_or_mis_990_AM"),
        (pl.col("lof_all") | pl.col("mis_995_AM")).alias("lof_or_mis_995_AM"),
        (pl.col("lof_all") | pl.col("mis_997_AM")).alias("lof_or_mis_997_AM"),
        (pl.col("lof_all") | pl.col("mis_999_AM")).alias("lof_or_mis_999_AM"),
        (pl.col("lof_all") | pl.col("mis_youden_MPC")).alias("lof_or_youden_MPC"),
        (pl.col("lof_all") | pl.col("mis_greater_than_3_MPC")).alias("lof_or_mis_greater_than_3_MPC"),
        (pl.col("lof_all") | (pl.col("Clinvar_Pathogenic") == True)).alias("lof_or_Clinvar_Path"),
    ])

    mask_cols = [
        "mis_906_AM", "mis_990_AM", "mis_995_AM", "mis_997_AM", "mis_999_AM",
        "mis_youden_MPC", "mis_greater_than_2_MPC", "mis_greater_than_3_MPC",
        "syn_only", "syn_mis_mid_SpliceAI", "syn_mis_high_SpliceAI", "syn_benign_SpliceAI",
        "Splice_SNV", "Stop_gained", "frameshift", "Clinvar_Path",
        "lof_snv", "lof_all",
        "lof_or_mis_906_AM", "lof_or_mis_990_AM", "lof_or_mis_995_AM",
        "lof_or_mis_997_AM", "lof_or_mis_999_AM",
        "lof_or_youden_MPC", "lof_or_mis_greater_than_3_MPC", "lof_or_Clinvar_Path",
    ]
    return df, mask_cols


def _weighted_or(a, b, c, d):
    num = 0.0
    den = 0.0
    for ai, bi, ci, di in zip(a, b, c, d):
        if any(pd.isna(x) for x in (ai, bi, ci, di)):
            continue
        n = ai + bi + ci + di
        if n > 0:
            num += (ai * di) / n
            den += (bi * ci) / n
    return (num / den) if den != 0 else None


def _filter_valid_indices(a, b, c, d):
    valid_idx = []
    for i, (ai, bi, ci, di) in enumerate(zip(a, b, c, d)):
        if any(pd.isna(x) for x in (ai, bi, ci, di)):
            continue
        if any(x < 0 for x in (ai, bi, ci, di)):
            continue
        if (ai > 0 or ci > 0) and (bi > 0 or di > 0) and (ai > 0 or bi > 0) and (ci > 0 or di > 0):
            valid_idx.append(i)
    return valid_idx


def _compute_exact_cmh_for_gene_r(payload):
    gene_id = payload["gene_id"]
    a = payload["a"]; b = payload["b"]; c = payload["c"]; d = payload["d"]

    try:
        valid_idx = _filter_valid_indices(a, b, c, d)
        if not valid_idx:
            return {
                "info.gene_id": gene_id,
                "exact_cmh_pvalue": None,
                "fisher_p_values": [],
                "weighted_or": None,
                "or_ci_lower": None,
                "or_ci_upper": None,
                "py_weighted_or": None,
                "or_source": None,
                "total_strata": 0,
            }

        fisher_p = []
        fisher_results = []
        tables = []
        for i in valid_idx:
            tab = np.array([[a[i], b[i]], [c[i], d[i]]], dtype=float)
            tables.append(tab)
            try:
                res = stats.fisher_test(
                    ro.r.matrix(ro.FloatVector(tab.ravel()), nrow=2),
                    alternative="two.sided",
                )
                fisher_p.append(float(res.rx2("p.value")[0]))
                fisher_results.append(res)
            except Exception:
                fisher_p.append(1.0)
                fisher_results.append(None)

        valid_a = [a[i] for i in valid_idx]
        valid_b = [b[i] for i in valid_idx]
        valid_c = [c[i] for i in valid_idx]
        valid_d = [d[i] for i in valid_idx]
        py_weighted_or = _weighted_or(valid_a, valid_b, valid_c, valid_d)

        # Single stratum → Fisher OR + CI
        if len(tables) == 1:
            exact_p = fisher_p[0] if len(fisher_p) == 1 else None
            single_tab = tables[0]
            num = single_tab[0, 0] * single_tab[1, 1]
            den = single_tab[0, 1] * single_tab[1, 0]
            if den == 0 and num > 0:
                weighted_or = float("inf")
            elif num == 0 and den > 0:
                weighted_or = 0.0
            elif den > 0:
                weighted_or = num / den
            else:
                weighted_or = None

            or_ci_lower = None
            or_ci_upper = None
            if fisher_results[0] is not None:
                try:
                    ci = fisher_results[0].rx2("conf.int")
                    or_ci_lower = float(ci[0])
                    or_ci_upper = float(ci[1])
                except Exception:
                    pass

            return {
                "info.gene_id": gene_id,
                "exact_cmh_pvalue": exact_p,
                "fisher_p_values": fisher_p,
                "weighted_or": weighted_or,
                "or_ci_lower": or_ci_lower,
                "or_ci_upper": or_ci_upper,
                "py_weighted_or": py_weighted_or,
                "or_source": "fisher",
                "total_strata": len(fisher_p),
            }

        # ≥2 strata → exact CMH with CI
        flat = np.array(tables, dtype=float).ravel()
        r_arr = ro.r.array(ro.FloatVector(flat), dim=ro.IntVector([2, 2, len(tables)]))

        try:
            cmh = stats.mantelhaen_test(r_arr, exact=True, correct=False, alternative="two.sided")
            exact_p = float(cmh.rx2("p.value")[0])
            mh_or = float(cmh.rx2("estimate")[0])
            ci = cmh.rx2("conf.int")
            or_ci_lower = float(ci[0])
            or_ci_upper = float(ci[1])
            weighted_or = mh_or
            or_source = "cmh"
        except Exception as e:
            return {
                "info.gene_id": gene_id,
                "error": f"mantelhaen_test failed: {e}",
                "exact_cmh_pvalue": None,
                "fisher_p_values": fisher_p,
                "weighted_or": py_weighted_or,
                "or_ci_lower": None,
                "or_ci_upper": None,
                "py_weighted_or": py_weighted_or,
                "or_source": "mh_python",
                "total_strata": len(fisher_p),
            }

        return {
            "info.gene_id": gene_id,
            "exact_cmh_pvalue": exact_p,
            "fisher_p_values": fisher_p,
            "weighted_or": weighted_or,
            "or_ci_lower": or_ci_lower,
            "or_ci_upper": or_ci_upper,
            "py_weighted_or": py_weighted_or,
            "or_source": or_source,
            "total_strata": len(fisher_p),
        }

    except Exception as e:
        return {
            "info.gene_id": gene_id,
            "error": f"R worker error: {e}",
            "exact_cmh_pvalue": None,
            "fisher_p_values": [],
            "weighted_or": None,
            "or_ci_lower": None,
            "or_ci_upper": None,
            "py_weighted_or": None,
            "or_source": None,
            "total_strata": 0,
        }


def exact_cmh_test_polars_parallel(ht: pl.DataFrame, n_workers: int = 30) -> pl.DataFrame:
    grouped_ht = ht.group_by(["info.gene_id", "info.group"]).agg([
        pl.col("info.ac_case").sum().alias("ac_case"),
        pl.col("info.ac_ctrl").sum().alias("ac_ctrl"),
        pl.col("info.an_case").max().alias("an_case"),
        pl.col("info.an_ctrl").max().alias("an_ctrl"),
    ])

    totals_ht = grouped_ht.group_by("info.gene_id").agg([
        pl.col("ac_case").sum().alias("total_ac_case"),
        pl.col("ac_ctrl").sum().alias("total_ac_ctrl"),
        pl.col("an_case").sum().alias("total_an_case"),
        pl.col("an_ctrl").sum().alias("total_an_ctrl"),
    ])

    wide_ht = grouped_ht.pivot(
        index="info.gene_id",
        on="info.group",
        values=["ac_case", "ac_ctrl", "an_case", "an_ctrl"],
        separator="_",
    )

    per_gene_payloads = (
        grouped_ht.with_columns([
            pl.col("ac_case").alias("a"),
            pl.col("ac_ctrl").alias("b"),
            (pl.col("an_case") - pl.col("ac_case")).alias("c"),
            (pl.col("an_ctrl") - pl.col("ac_ctrl")).alias("d"),
        ])
        .group_by("info.gene_id")
        .agg([pl.col("a"), pl.col("b"), pl.col("c"), pl.col("d")])
        .select([pl.col("info.gene_id").alias("gene_id"), "a", "b", "c", "d"])
        .to_dicts()
    )

    results, errors = _run_pool_r_joblib(per_gene_payloads, n_workers=n_workers)

    failed = [e["info.gene_id"] for e in errors if e.get("error")]
    if failed:
        print(f"[ERROR] Exact CMH failed for {len(failed)} gene(s). Examples: {failed[:10]}")
        raise RuntimeError(f"Exact CMH (R) failed for {len(failed)} genes. First: {failed[0]}")

    cmh_df = pl.from_dicts(results)

    return (
        cmh_df
        .join(totals_ht, on="info.gene_id", how="left")
        .join(wide_ht, on="info.gene_id", how="left")
    )


def load_gene_map(hgnc_path: str | os.PathLike = _DEFAULT_HGNC, return_dict: bool = True):
    hgnc = pd.read_csv(str(hgnc_path), sep="\t")
    if return_dict:
        return dict(zip(hgnc["Ensembl gene ID"], hgnc["Approved symbol"]))
    return hgnc


def _process_single_mis_col(args):
    mis, df, map_ensg = args
    print(f"[{mis}] Starting processing")

    filtered_df = df.filter(pl.col(mis) == True)
    if filtered_df.height == 0:
        print(f"[{mis}] No variants found")
        return None, None

    variant_counts = (
        filtered_df.group_by(["info.gene_id", "ID_38"])
        .agg([
            pl.len().alias("variant_occurrence"),
            pl.col("info.ac_case").sum().alias("total_ac_case"),
            pl.col("info.ac_ctrl").sum().alias("total_ac_ctrl"),
        ])
    )

    recurrent_variants = (
        variant_counts.filter(pl.col("variant_occurrence") > 1)
        .group_by("info.gene_id")
        .agg([pl.struct(["ID_38", "variant_occurrence"]).alias("recurrent_variants")])
    )

    gene_metrics = (
        filtered_df.group_by("info.gene_id").agg([
            pl.n_unique("ID_38").alias("unique_variants"),
            (pl.col("info.ac_case").sum() + pl.col("info.ac_ctrl").sum()).alias("total_qualifying_counts"),
        ])
    )

    strata_counts = (
        filtered_df.group_by("info.gene_id").agg([
            pl.n_unique("info.group").alias("total_strata"),
            pl.col("info.group").str.contains("EUR").sum().alias("euro_strata"),
        ])
    )

    strata_ors = (
        filtered_df.group_by(["info.gene_id", "info.group"])
        .agg([
            ((pl.col("info.ac_case") / pl.col("info.an_case"))
             / (pl.col("info.ac_ctrl") / pl.col("info.an_ctrl"))).alias("strata_or")
        ])
        .group_by("info.gene_id")
        .agg([pl.struct(["info.group", "strata_or"]).alias("strata_specific_ors")])
    )

    test = exact_cmh_test_polars_parallel(filtered_df, n_workers=os.cpu_count())
    if test.height == 0:
        print(f"[{mis}] No CMH results")
        return None, filtered_df

    test = (
        test
        .join(gene_metrics, on="info.gene_id", how="left")
        .join(recurrent_variants, on="info.gene_id", how="left")
        .join(strata_counts, on="info.gene_id", how="left")
        .join(strata_ors, on="info.gene_id", how="left")
    )
    test = test.with_columns(pl.col("info.gene_id").replace(map_ensg).alias("Gene_Name"))

    # Rename per-stratum AC/AN columns from "ac_case_<stratum>" → "<stratum>_ac_case"
    other_cols = test.select(~cs.matches("^ac|^an")).columns
    renames = {}
    for c in test.columns:
        if c in other_cols:
            continue
        m = re.search(r"^a._c..._(.+)", c)
        if m:
            renames[c] = f"{m.group(1)}_{c.split('_')[0]}_{c.split('_')[1]}"
    if renames:
        test = test.rename(renames)

    test = test.select(other_cols + sorted([c for c in test.columns if c not in other_cols]))
    test = test.with_columns(pl.lit(mis).alias("Name"))
    test = test.select(
        pl.col(["info.gene_id", "Name"]),
        pl.col([c for c in test.columns if c not in ("info.gene_id", "Name")]),
    )

    print(f"[{mis}] Completed: {test.shape}")
    return test, filtered_df


# Column-keep regex — only retains what the CMH backend actually reads.
_KEEP_RE = re.compile(
    r"^ID_38$|^info\.|^mane_select\.(transcript_id|consequence_terms)$"
    r"|^AlphaMissense\.am_pathogenicity$|^dbNSFP\.MPC_score$|^SpliceAI_max$|^SNV$"
    r"|^clinvar_vcf\.(CLNSIG|CLNREVSTAT)$|^Clinvar_Pathogenic$"
)


def run_hmbs_collapsing(
    file_pattern: str,
    n_parallel_mis_cols: int = 1,
    clinvar_mode: str = "hq_plp",
    row_ac_max: int | None = 5,
):
    """Run the HMBS CMH burden test against one or more annotated parquets.

    Parameters
    ----------
    file_pattern
        Glob pattern pointing at one or more annotated input parquets. All matched
        parquets are concatenated (diagonal_relaxed) before the analysis.
    n_parallel_mis_cols
        Number of masks to process in parallel.
    clinvar_mode
        "plp" or "hq_plp". See _with_clinvar_pathogenic_flag. Default "hq_plp"
        matches the reference HMBS workbook (SCHEMA Clinvar_Path 3/3).
    row_ac_max
        If not None, drop variants whose total case + control AC across strata exceeds
        this value. Default 5 matches the upstream pipeline.

    Returns
    -------
    (out_cmh, all_dfs)
        out_cmh: dict with a single key "OR_pass_No_Tier" → pl.DataFrame of CMH results.
        all_dfs: list of per-mask filtered source DataFrames (useful for debugging).
    """
    matches = glob.glob(file_pattern)
    if not matches:
        raise FileNotFoundError(f"No parquets matched: {file_pattern}")

    df = pl.concat([pl.read_parquet(p) for p in matches], how="diagonal_relaxed").lazy()

    cols = [c for c in df.collect_schema().names() if _KEEP_RE.search(c)]
    df = df.select(pl.col(cols))

    df, mis_cols = create_missense_columns(df, clinvar_mode=clinvar_mode)

    # Row-AC filter: drop variants with total case+ctrl AC (summed across strata) > row_ac_max
    if row_ac_max is not None:
        df = df.with_columns((pl.col("info.ac_case") + pl.col("info.ac_ctrl")).alias("info.row_ac"))
        df = df.with_columns(pl.col("info.row_ac").sum().over("ID_38").alias("info.row_ac"))
        df = df.filter(pl.col("info.row_ac") <= row_ac_max)

    df = df.collect()
    map_ensg = load_gene_map()

    print(f"Processing {len(mis_cols)} missense columns")
    print(f"Running {min(n_parallel_mis_cols, len(mis_cols))} mis_cols in parallel")

    args_list = [(mis, df, map_ensg) for mis in mis_cols]

    results = Parallel(
        n_jobs=min(n_parallel_mis_cols, 4),
        backend="threading",
        prefer="threads",
        verbose=0,
    )(delayed(_process_single_mis_col)(args) for args in args_list)

    out_cmh: dict[str, pl.DataFrame] = {}
    all_dfs = []
    key = "OR_pass_No_Tier"
    for test, filtered_df in results:
        if test is not None:
            if key not in out_cmh:
                out_cmh[key] = test
            else:
                out_df = pl.concat([out_cmh[key], test], how="diagonal_relaxed")
                out_cmh[key] = out_df.select(pl.col(["Gene_Name"] + out_df.drop("Gene_Name").columns))
        if filtered_df is not None:
            all_dfs.append(filtered_df)

    for k in out_cmh:
        out_cmh[k] = out_cmh[k].sort("Name", "exact_cmh_pvalue", descending=[True, False])

    gc.collect()
    return out_cmh, all_dfs
