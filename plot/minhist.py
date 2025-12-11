#!/usr/bin/env -S uv run
# /// script
# dependencies = [
#   "pandas",
#   "altair",
#   "vl-convert-python",
# ]
# ///

import argparse
import json
import os
import re
import sys

import altair as alt
import pandas as pd


def natural_sort_key(text):
    """Generate a sort key for natural sorting of strings with numbers."""
    return tuple(int(c) if c.isdigit() else c.lower() for c in re.split(r"(\d+)", str(text)))


def load_data(input_file):
    """Load data from JSON or CSV file."""
    with open(input_file, 'r') as f:
        # Try to detect if it's JSON by checking first character
        first_char = f.read(1)
        f.seek(0)

        if first_char == '{':
            # Parse JSON
            data = json.load(f)
            rows = []

            # Extract data from JSON structure
            for sample in data.get("samples", []):
                sample_name = sample.get("sample_name", "unknown")
                for target_data in sample.get("targets", []):
                    # Convert abundance_histogram from [[depth, count], ...] to dict
                    hist_list = target_data.get("abundance_histogram", [])
                    hist_dict = {str(depth): count for depth, count in hist_list}

                    row = {
                        "sample": sample_name,
                        "target": target_data.get("target", ""),
                        "length": target_data.get("length", 0),
                        "total_minimizers": target_data.get("total_minimizers", 0),
                        "contained_minimizers": target_data.get("contained_minimizers", 0),
                        "containment1": target_data.get("containment1", 0.0),
                        "median_nz_abundance": target_data.get("median_nz_abundance", 0.0),
                        "abundance_histogram": hist_dict
                    }
                    rows.append(row)

            return pd.DataFrame(rows)
        else:
            # Parse CSV
            return pd.read_csv(input_file)


def parse_histogram(hist_data):
    """Parse histogram data to dictionary."""
    if pd.isna(hist_data):
        return {}
    if isinstance(hist_data, dict):
        return hist_data
    if isinstance(hist_data, list):
        # Handle [[depth, count], ...] format
        return {str(depth): count for depth, count in hist_data}
    try:
        # Try parsing as JSON string
        parsed = json.loads(hist_data)
        if isinstance(parsed, list):
            return {str(depth): count for depth, count in parsed}
        return parsed
    except (json.JSONDecodeError, TypeError):
        # Try parsing as string representation
        try:
            return eval(hist_data)
        except:
            return {}


def main():
    parser = argparse.ArgumentParser(description="Plot Grate minimizer abundance histograms")
    parser.add_argument("input_file", help="Input JSON or CSV file containing histogram data")
    parser.add_argument("--output-plot", help="Output plot filename (default: <input_prefix>_hist.png)")
    parser.add_argument("--force", action="store_true", help="Overwrite existing output files")
    parser.add_argument("--debug", action="store_true", help="Enable debug output")
    parser.add_argument("--title", default="Minimizer abundance histograms",
                        help="Plot title (default: %(default)s)")
    parser.add_argument("--short-names", action="store_true",
                        help="Remove accession prefix (before first space) from target names")
    parser.add_argument("--max-depth", type=int, default=None,
                        help="Maximum depth to display (default: auto-detect)")
    parser.add_argument("--min-depth", type=int, default=1,
                        help="Minimum depth to display (default: 1, excludes zero-depth)")
    parser.add_argument("--linear", action="store_true",
                        help="Use linear scale for y-axis instead of log scale")

    args = parser.parse_args()

    # Auto-generate output filename from input prefix if not specified
    input_filename = os.path.basename(args.input_file)
    input_prefix = os.path.splitext(input_filename)[0]
    if args.output_plot is None:
        args.output_plot = f"{input_prefix}-hist.png"

    # Check if output file exists and require --force to overwrite
    if os.path.exists(args.output_plot) and not args.force:
        print(f"ERROR: Output file already exists: {args.output_plot}")
        print("Use --force to overwrite existing files")
        sys.exit(1)

    try:
        # --- Load & validate ---
        if not os.path.exists(args.input_file):
            print(f"ERROR: File {args.input_file} does not exist")
            sys.exit(1)
        if os.path.getsize(args.input_file) == 0:
            print(f"ERROR: File {args.input_file} is empty")
            sys.exit(1)

        df = load_data(args.input_file)

        if args.debug:
            print(f"\n{args.input_file}:")
            print(f"  Columns: {list(df.columns)}")
            print(f"  Shape: {df.shape}")
            print("  First few rows:")
            print(df.head())

        if df.empty:
            print("ERROR: Input file has no data")
            sys.exit(1)

        # Ensure required columns exist
        required_cols = {"target", "abundance_histogram", "sample"}
        missing = [c for c in required_cols if c not in df.columns]
        if missing:
            print(f"ERROR: Input data missing required columns: {', '.join(missing)}")
            sys.exit(1)

        # Display name tweaks
        if args.short_names:
            df["display_name"] = df["target"].apply(
                lambda x: str(x).split(" ", 1)[1] if isinstance(x, str) and " " in x else x
            )
            df["sort_key"] = df["display_name"].apply(natural_sort_key)
        else:
            df["display_name"] = df["target"]
            df["sort_key"] = df["target"].apply(natural_sort_key)

        # Parse histogram data and expand to long format
        rows = []
        for _, row in df.iterrows():
            hist = parse_histogram(row["abundance_histogram"])
            if not hist:
                continue

            for depth_str, count in hist.items():
                try:
                    depth = int(depth_str)
                    # Apply depth filters
                    if depth < args.min_depth:
                        continue
                    if args.max_depth is not None and depth > args.max_depth:
                        continue

                    rows.append({
                        "target": row["target"],
                        "display_name": row["display_name"],
                        "sort_key": row["sort_key"],
                        "sample": row["sample"],
                        "depth": depth,
                        "count": count
                    })
                except (ValueError, TypeError):
                    continue

        if not rows:
            print("ERROR: No valid histogram data found")
            sys.exit(1)

        hist_df = pd.DataFrame(rows)

        # Calculate frequency (normalize counts within each target-sample group)
        hist_df["frequency"] = hist_df.groupby(["target", "sample"])["count"].transform(
            lambda x: x / x.sum()
        )

        if args.debug:
            print(f"\nExpanded histogram data:")
            print(f"  Shape: {hist_df.shape}")
            print(f"  Depth range: {hist_df['depth'].min()} to {hist_df['depth'].max()}")
            print(f"  First few rows:")
            print(hist_df.head())

        # Ordering
        hist_df = hist_df.sort_values(["sort_key", "sample", "depth"])
        sample_order = sorted(hist_df["sample"].dropna().astype(str).unique(), key=natural_sort_key)
        target_order = sorted(hist_df["display_name"].dropna().astype(str).unique(), key=natural_sort_key)

        # --- Plot ---
        alt.data_transformers.enable("json")

        # Determine y-axis scale
        y_scale = alt.Scale(type="linear") if args.linear else alt.Scale(type="log")

        base = (
            alt.Chart(hist_df)
            .mark_line(point=True, strokeWidth=2)
            .encode(
                x=alt.X("depth:Q",
                       title="Abundance",
                       scale=alt.Scale(zero=False)),
                y=alt.Y("frequency:Q",
                       title="Frequency",
                       scale=y_scale),
                color=alt.Color("sample:N",
                              title="Sample",
                              sort=sample_order,
                              scale=alt.Scale(scheme="category10")),
                tooltip=[
                    alt.Tooltip("target:N", title="target"),
                    alt.Tooltip("sample:N", title="sample"),
                    alt.Tooltip("depth:Q", title="depth"),
                    alt.Tooltip("count:Q", title="count", format=",.0f"),
                    alt.Tooltip("frequency:Q", title="frequency", format=".4f"),
                ],
            )
            .properties(
                width=600,
                height=80
            )
        )

        chart = (
            base.facet(
                row=alt.Row("display_name:N",
                           title=None,
                           sort=target_order,
                           header=alt.Header(labelAlign="left", labelAnchor="start", labelFontSize=11))
            )
            .resolve_scale(y="independent")
            .properties(title=args.title)
            .configure_legend(titleFontSize=12, labelFontSize=11)
            .configure_axis(labelFontSize=10, titleFontSize=12)
            .configure_title(fontSize=14)
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
