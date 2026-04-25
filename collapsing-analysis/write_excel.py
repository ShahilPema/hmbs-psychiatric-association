"""Standalone Excel writer vendored from src/utils.py of Missense_Predictor."""

import pandas as pd


def _col_idx_to_excel(col_idx):
    """Convert 0-based column index to Excel column letter(s) (A, B, ..., Z, AA, AB, ...)."""
    result = ""
    idx = col_idx
    while True:
        result = chr(65 + idx % 26) + result
        idx = idx // 26 - 1
        if idx < 0:
            break
    return result


def _write_column(sheet, writer, df, max_allowed_length):
    """Optimized column width calculation using DataFrame data directly."""
    worksheet = writer.sheets[sheet]

    for idx, column in enumerate(df.columns):
        header_length = len(str(column))

        sample_size = min(1000, len(df))
        if len(df) > sample_size:
            sample_data = df[column].sample(n=sample_size, random_state=42)
        else:
            sample_data = df[column]

        if pd.api.types.is_numeric_dtype(sample_data):
            valid_values = [val for val in sample_data.dropna().iloc[:100] if pd.notna(val)]
            max_data_length = max((len(str(val)) for val in valid_values), default=5)
        else:
            max_data_length = sample_data.astype(str).str.len().max()
            if pd.isna(max_data_length):
                max_data_length = 5

        column_width = min(max(header_length, max_data_length) + 2, max_allowed_length)
        worksheet.set_column(idx, idx, column_width)


def write_excel(write_dict, name="", max_allowed_length=40, color_columns=None):
    """Write dict of DataFrames to Excel with auto-fitted columns and optional color scales."""
    with pd.ExcelWriter(name, engine="xlsxwriter") as writer:
        writer.book.set_size(1920, 1080)

        for k, v in write_dict.items():
            v.to_excel(writer, sheet_name=k, index=False)
            _write_column(k, writer, v, max_allowed_length)

            if color_columns and k in color_columns:
                worksheet = writer.sheets[k]
                sheet_spec = color_columns[k]

                if isinstance(sheet_spec, list):
                    col_directions = {col: "higher" for col in sheet_spec}
                else:
                    col_directions = sheet_spec

                valid_columns = {col: d for col, d in col_directions.items() if col in v.columns}
                for col_name, direction in valid_columns.items():
                    col_idx = v.columns.get_loc(col_name)
                    col_letter = _col_idx_to_excel(col_idx)
                    col_range = f"{col_letter}2:{col_letter}{len(v) + 1}"

                    if direction == "lower":
                        worksheet.conditional_format(col_range, {
                            "type": "3_color_scale",
                            "min_color": "#63BE7B",
                            "mid_color": "#FFEB84",
                            "max_color": "#F8696B",
                        })
                    else:
                        worksheet.conditional_format(col_range, {
                            "type": "3_color_scale",
                            "min_color": "#F8696B",
                            "mid_color": "#FFEB84",
                            "max_color": "#63BE7B",
                        })

    print(f"Data written to {name} and columns autofitted!")
