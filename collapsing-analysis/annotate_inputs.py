"""Pre-annotate the raw SCHEMA / BipEx / Combined HMBS parquets.

Runs ONCE per ClinVar snapshot bump. Produces self-contained parquets that the
end-user `run_collapsing.py` driver consumes without needing access to the ClinVar
VCF or the AoU whitelist again.

Each output parquet has the following bits baked in:
- Variants restricted to the AoU whitelist (HMBS_Combined_AoUFilt.parquet).
- ClinVar annotations from the supplied VCF, joined onto every variant as
  `clinvar_vcf.CLNSIG`, `clinvar_vcf.CLNREVSTAT`, `clinvar_vcf.ALLELEID`,
  `clinvar_vcf.RS`.
- `SNV` boolean column precomputed from the `alleles` list length.

`AlphaMissense.am_pathogenicity` in the input parquet is used as-is — verified to
match the upstream predictions checkpoint exactly for all HMBS variants where
both are non-null (max |diff| = 0, identical null patterns).

A sidecar `annotation_manifest.json` records the provenance (paths, MD5s,
timestamp, git SHA) so `run_collapsing.py` can log what snapshot the analysis
was run against.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import logging
import subprocess
import time
from pathlib import Path

import polars as pl

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


def _file_md5(path: Path, block_size: int = 1 << 20) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(block_size), b""):
            h.update(chunk)
    return h.hexdigest()


def _git_sha(here: Path) -> str | None:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=str(here),
            stderr=subprocess.DEVNULL,
        ).decode().strip()
    except Exception:
        return None


def load_hmbs_clinvar_from_vcf(vcf_path: Path) -> pl.DataFrame:
    """Parse a ClinVar VCF and extract HMBS records.

    Keeps only variants whose INFO contains `GENEINFO=HMBS:3145`. Returns one row per
    unique ID_38 (`{chrom}-{pos}-{ref}-{alt}`) with CLNSIG, CLNREVSTAT, ALLELEID, RS.
    """
    records = []
    with gzip.open(str(vcf_path), "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t", 8)
            if len(parts) < 8:
                continue
            chrom, pos, _vid, ref, alt, _qual, _flt, info = parts[:8]
            if "GENEINFO=HMBS:3145" not in info:
                continue
            info_map = {}
            for field in info.split(";"):
                if "=" not in field:
                    continue
                key, value = field.split("=", 1)
                info_map[key] = value
            records.append({
                "ID_38": f"{chrom.replace('chr', '')}-{pos}-{ref}-{alt}",
                "clinvar_vcf.CLNSIG": info_map.get("CLNSIG"),
                "clinvar_vcf.CLNREVSTAT": info_map.get("CLNREVSTAT"),
                "clinvar_vcf.ALLELEID": info_map.get("ALLELEID"),
                "clinvar_vcf.RS": info_map.get("RS"),
            })

    if not records:
        return pl.DataFrame(schema={
            "ID_38": pl.Utf8,
            "clinvar_vcf.CLNSIG": pl.Utf8,
            "clinvar_vcf.CLNREVSTAT": pl.Utf8,
            "clinvar_vcf.ALLELEID": pl.Utf8,
            "clinvar_vcf.RS": pl.Utf8,
        })
    return pl.DataFrame(records).unique(subset=["ID_38"], keep="first")


def annotate_one(
    input_parquet: Path,
    aou_ids: list[str],
    clinvar_df: pl.DataFrame,
    output_parquet: Path,
) -> dict:
    """Annotate one HMBS cohort parquet and write the annotated output."""
    logger.info(f"[{input_parquet.name}] reading")
    df = pl.read_parquet(input_parquet)
    n_in = df.height

    df = df.filter(pl.col("ID_38").is_in(aou_ids))
    n_aou = df.height
    logger.info(f"[{input_parquet.name}] after AoU filter: {n_in} -> {n_aou}")

    df = df.join(clinvar_df, on="ID_38", how="left", coalesce=True)
    logger.info(f"[{input_parquet.name}] ClinVar VCF annotations joined")

    df = df.with_columns(
        pl.when(pl.col("alleles").list[0].len() != pl.col("alleles").list[1].len())
          .then(False).otherwise(True).alias("SNV")
    )

    output_parquet.parent.mkdir(parents=True, exist_ok=True)
    df.write_parquet(output_parquet)
    logger.info(f"[{input_parquet.name}] wrote {output_parquet} ({df.height} rows, {len(df.columns)} cols)")

    return {
        "input_parquet": str(input_parquet),
        "output_parquet": str(output_parquet),
        "rows_in": n_in,
        "rows_after_aou": n_aou,
        "rows_out": df.height,
        "cols_out": len(df.columns),
    }


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--schema-parquet", required=True, type=Path, help="Raw SCHEMA-only HMBS parquet.")
    p.add_argument("--bipex-parquet", required=True, type=Path, help="Raw BipEx-only HMBS parquet.")
    p.add_argument("--combined-parquet", required=True, type=Path, help="Raw Combined HMBS parquet.")
    p.add_argument("--clinvar-vcf", required=True, type=Path, help="ClinVar VCF snapshot (.vcf.gz).")
    p.add_argument("--aou-filter", required=True, type=Path, help="Parquet with ID_38 column listing AoU-passing variants.")
    p.add_argument("--output-dir", required=True, type=Path, help="Where to write the annotated parquets.")
    return p.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Loading AoU whitelist from {args.aou_filter}")
    aou_ids = pl.read_parquet(args.aou_filter)["ID_38"].unique().to_list()
    logger.info(f"AoU whitelist: {len(aou_ids)} unique variant IDs")

    logger.info(f"Parsing ClinVar VCF {args.clinvar_vcf}")
    clinvar_df = load_hmbs_clinvar_from_vcf(args.clinvar_vcf)
    logger.info(f"ClinVar VCF: {clinvar_df.height} HMBS records")

    jobs = [
        (args.schema_parquet, args.output_dir / "HMBS_SCHEMA_ONLY.annotated.parquet"),
        (args.bipex_parquet, args.output_dir / "HMBS_Bipex_ONLY.annotated.parquet"),
        (args.combined_parquet, args.output_dir / "HMBS_Combined.annotated.parquet"),
    ]

    manifest = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "git_sha": _git_sha(Path(__file__).resolve().parent),
        "clinvar_vcf": {
            "path": str(args.clinvar_vcf),
            "md5": _file_md5(args.clinvar_vcf),
            "size_bytes": args.clinvar_vcf.stat().st_size,
            "hmbs_records": clinvar_df.height,
        },
        "aou_filter": {
            "path": str(args.aou_filter),
            "md5": _file_md5(args.aou_filter),
            "unique_ids": len(aou_ids),
        },
        "cohorts": [],
    }

    for input_parquet, output_parquet in jobs:
        manifest["cohorts"].append(
            annotate_one(input_parquet, aou_ids, clinvar_df, output_parquet)
        )

    manifest_path = args.output_dir / "annotation_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))
    logger.info(f"Wrote manifest to {manifest_path}")


if __name__ == "__main__":
    main()
