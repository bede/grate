#!/usr/bin/env python3

import argparse
import os
import re
import sys

import altair as alt
import pandas as pd


def natural_sort_key(text):
    """Generate a sort key for natural sorting of strings with numbers."""
    return tuple(int(c) if c.isdigit() else c.lower() for c in re.split(r"(\d+)", str(text)))


def main():
    parser = argparse.ArgumentParser(description="Plot Grate CSV file (single CSV; one or many samples)")
    parser.add_argument("input_csv", help="Input CSV file containing one or many samples (must include a sample column)")
    parser.add_argument("--output-plot", default="containment.png", help="Output plot filename")
    parser.add_argument("--output-csv", default="containment.csv", help="(Re)written CSV filename")
    parser.add_argument("--debug", action="store_true", help="Enable debug output")
    parser.add_argument("--title", default="Containment analysis (Grate)", help="Plot title")
    parser.add_argument("--short-names", action="store_true",
                        help="Remove accession prefix (before first space) from target names")
    parser.add_argument("--sample-column", default="sample",
                        help="Name of the column that identifies samples (default: 'sample')")

    args = parser.parse_args()

    try:
        # --- Load & validate ---
        if not os.path.exists(args.input_csv):
            print(f"ERROR: File {args.input_csv} does not exist")
            sys.exit(1)
        if os.path.getsize(args.input_csv) == 0:
            print(f"ERROR: File {args.input_csv} is empty")
            sys.exit(1)

        df = pd.read_csv(args.input_csv)

        if args.debug:
            print(f"\n{args.input_csv}:")
            print(f"  Columns: {list(df.columns)}")
            print(f"  Shape: {df.shape}")
            print("  First few rows:")
            print(df.head())

        if df.empty:
            print("ERROR: Input CSV has no data")
            sys.exit(1)

        # Ensure required columns exist
        required_cols = {"target", "containment", "median_abundance"}
        missing = [c for c in required_cols if c not in df.columns]
        if missing:
            print(f"ERROR: Input CSV missing required columns: {', '.join(missing)}")
            sys.exit(1)

        sample_col = args.sample_column
        if sample_col not in df.columns:
            print(f"ERROR: Input CSV missing sample column '{sample_col}'. "
                  f"Use --sample-column to specify the correct column.")
            sys.exit(1)

        # Coerce numeric columns in case they're strings
        for c in ("containment", "median_abundance", "length_bp", "contained_minimizers"):
            if c in df.columns:
                df[c] = pd.to_numeric(df[c], errors="coerce")

        # Keep a (re)written copy for reproducibility
        df.to_csv(args.output_csv, index=False)
        print(f"CSV saved to: {args.output_csv}")

        # --- Prep for plotting ---
        plot_df = df.copy()

        # Optional: drop 'ALL' row if present
        if "target" in plot_df.columns:
            plot_df = plot_df[plot_df["target"] != "ALL"]

        if plot_df.empty:
            print("ERROR: No data to plot after filtering")
            sys.exit(1)

        # Display name tweaks
        if args.short_names:
            plot_df["display_name"] = plot_df["target"].apply(
                lambda x: str(x).split(" ", 1)[1] if isinstance(x, str) and " " in x else x
            )
            plot_df["sort_key"] = plot_df["display_name"].apply(natural_sort_key)
        else:
            plot_df["display_name"] = plot_df["target"]
            plot_df["sort_key"] = plot_df["target"].apply(natural_sort_key)

        # Ordering
        plot_df = plot_df.sort_values(["sort_key", sample_col])
        plot_df = plot_df.drop("sort_key", axis=1)

        sample_order = sorted(plot_df[sample_col].dropna().astype(str).unique(), key=natural_sort_key)
        target_order = sorted(plot_df["display_name"].dropna().astype(str).unique(), key=natural_sort_key)

        # Text labels for median depth
        plot_df["depth_label"] = plot_df["median_abundance"].apply(
            lambda x: f"med(depth): {x:.0f}" if pd.notna(x) and float(x) > 0 else ""
        )

        # --- Plot ---
        alt.data_transformers.enable("json")

        bars = (
            alt.Chart(plot_df)
            .mark_bar(size=6)
            .encode(
                y=alt.Y("display_name:N", title="", sort=target_order),
                x=alt.X("containment:Q", title="Containment", scale=alt.Scale(domain=[0, 1])),
                color=alt.Color(
                    f"{sample_col}:N",
                    sort=sample_order,
                    title="",
                    scale=alt.Scale(scheme="category20"),
                ),
                yOffset=alt.YOffset(f"{sample_col}:N", sort=sample_order),
                tooltip=[
                    "target:N",
                    alt.Tooltip(f"{sample_col}:N", title="sample"),
                    "containment:Q",
                    alt.Tooltip("median_abundance:Q", title="median_abundance"),
                    alt.Tooltip("length_bp:Q", title="length_bp", format=",.0f") if "length_bp" in plot_df.columns else alt.value(None),
                    alt.Tooltip("contained_minimizers:Q", title="contained_minimizers", format=",.0f")
                    if "contained_minimizers" in plot_df.columns else alt.value(None),
                ],
            )
        )

        text_labels = (
            alt.Chart(plot_df)
            .mark_text(align="left", baseline="middle", dx=5, fontSize=7, color="black")
            .encode(
                y=alt.Y("display_name:N", sort=target_order),
                yOffset=alt.YOffset(f"{sample_col}:N", sort=sample_order),
                x=alt.value(0),
                text="depth_label:N",
            )
        )

        chart = (
            (bars + text_labels)
            .properties(title=args.title, width=450, height=alt.Step(7))
            .configure_legend(titleFontSize=14, labelFontSize=12, symbolSize=150)
            .configure_axis(labelFontSize=12, titleFontSize=14)
            .configure_title(fontSize=14)
            .resolve_scale(y="shared")
        )

        chart.save(args.output_plot, scale_factor=2.0)
        print(f"Plot saved to: {args.output_plot}")

    except Exception as e:
        print(f"Error creating plot: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
