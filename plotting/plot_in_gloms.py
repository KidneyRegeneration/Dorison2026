"""
Generate boxplots and/or ridge plots from a conditions CSV file.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde


COLUMN_SUFFIX_TO_LABEL = {
    "_total_marker_per_glom_divided_by_glom_size": "_Normalised_by_Glom_Size",
    "total_marker_norm_by_minmax_of_nuclei_channel_divided_by_glom_size": "_Normalised_by_Glom_Size+DAPI",
    "_total_marker_norm_by_minmax_of_nuclei_channel": "_Normalised_by_DAPI_only",
}

TARGET_SUFFIX = "total_marker_per_glom_divided_by_glom_size"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot boxplots and/or ridge plots from a conditions CSV."
    )
    parser.add_argument(
        "-i", "--input-csv",
        required=True,
        help="Path to the input CSV file.",
    )
    parser.add_argument(
        "-b", "--boxplot-dir",
        default=None,
        help="Output directory for boxplots. Omit to skip boxplots.",
    )
    parser.add_argument(
        "-r", "--ridge-dir",
        default=None,
        help="Output directory for ridge plots. Omit to skip ridge plots.",
    )
    parser.add_argument(
        "--add-outlier-labels",
        action="store_true",
        default=False,
        help="Annotate outlier points with their replicate label.",
    )
    return parser.parse_args()


def resolve_column_label(column: str) -> str:
    """Map a raw column name to a human-readable axis label."""
    for suffix, label_suffix in COLUMN_SUFFIX_TO_LABEL.items():
        if column.endswith(suffix):
            return column.split("_")[0] + label_suffix
    return column


def annotate_outliers(
    ax: plt.Axes, df: pd.DataFrame, column: str, conditions: list[str]
) -> None:
    """Label outlier data points (IQR method) with their replicate name."""
    for x_pos, condition in enumerate(conditions):
        subset = df[df["condition"] == condition]
        q1, q3 = subset[column].quantile(0.25), subset[column].quantile(0.75)
        iqr = q3 - q1
        lower, upper = q1 - 1.5 * iqr, q3 + 1.5 * iqr
        outliers = subset[(subset[column] < lower) | (subset[column] > upper)]
        for _, row in outliers.iterrows():
            ax.text(
                x_pos, row[column], row["replicate"],
                fontsize=8, ha="left", va="bottom",
            )


def plot_boxplot(
    df: pd.DataFrame,
    column: str,
    column_label: str,
    conditions: list[str],
    output_dir: Path,
    add_outlier_labels: bool,
) -> None:
    """Save a boxplot with individual data points overlaid."""
    fig, ax = plt.subplots()
 
    has_experiment = "experiment" in df.columns and df["experiment"].notna().any()
    experiments = df["experiment"].unique().tolist() if has_experiment else []
    condition_colours = plt.cm.tab10(np.linspace(0, 1, len(conditions)))
    experiment_colours = plt.cm.tab10(np.linspace(0, 1, len(experiments))) if experiments else []
 
    # --- Boxplots ---
    grouped_data = [
        df[df["condition"] == condition][column].dropna().to_numpy()
        for condition in conditions
    ]
    bp = ax.boxplot(
        grouped_data,
        positions=range(len(conditions)),
        patch_artist=True,
        widths=0.5,
        showfliers=False,
    )
    for patch, colour in zip(bp["boxes"], condition_colours):
        patch.set_facecolor((*colour[:3], 0.3))
        patch.set_edgecolor(colour[:3])
    for element in ("whiskers", "caps", "medians"):
        for line, colour in zip(
            bp[element],
            np.repeat(condition_colours, 2 if element != "medians" else 1, axis=0),
        ):
            line.set_color(colour[:3])
 
    # --- Strip plot ---
    rng = np.random.default_rng(seed=42)
    has_experiment = "experiment" in df.columns and df["experiment"].notna().any()
 
    for x_pos, condition in enumerate(conditions):
        subset = df[df["condition"] == condition]
 
        if has_experiment:
            for exp_idx, experiment in enumerate(experiments):
                exp_subset = subset[subset["experiment"] == experiment][column].dropna()
                if exp_subset.empty:
                    continue
                jitter = rng.uniform(-0.15, 0.15, size=len(exp_subset))
                ax.scatter(
                    x_pos + jitter,
                    exp_subset,
                    alpha=0.5,
                    color=experiment_colours[exp_idx],
                    s=20,
                    zorder=3,
                )
        else:
            values = subset[column].dropna()
            jitter = rng.uniform(-0.15, 0.15, size=len(values))
            ax.scatter(
                x_pos + jitter,
                values,
                alpha=0.5,
                color=condition_colours[x_pos],
                s=20,
                zorder=3,
            )
 
    if add_outlier_labels:
        annotate_outliers(ax, df, column, conditions)
 
    ax.set_xticks(range(len(conditions)))
    ax.set_xticklabels(conditions, rotation=45, ha='right')
    ax.set_ylabel(column_label)
    ax.set_title(column_label)
 
    if experiments:
        legend_handles = [
            mpatches.Patch(color=experiment_colours[i], alpha=0.7, label=exp)
            for i, exp in enumerate(experiments)
        ]
        ax.legend(
            handles=legend_handles,
            title="experiment",
            loc="upper left",
            bbox_to_anchor=(1, 1),
        )
 
    fig.savefig(output_dir / f"{column}_boxplot.png", bbox_inches="tight", dpi=300)
    plt.close(fig)
 

def plot_ridge(
    df: pd.DataFrame,
    column: str,
    column_label: str,
    conditions: list[str],
    output_dir: Path,
) -> None:
    """Save a ridge (KDE overlay) plot."""
    fig, ax = plt.subplots()

    colours = plt.cm.tab10(np.linspace(0, 1, len(conditions)))
    x_min = df[column].min()
    x_max = df[column].max()
    x_grid = np.linspace(x_min, x_max, 300)

    for colour, condition in zip(colours, conditions):
        values = df[df["condition"] == condition][column].dropna().to_numpy()
        if len(values) < 2:
            continue
        kde = gaussian_kde(values)
        density = kde(x_grid)
        ax.plot(x_grid, density, color=colour, label=condition)
        ax.fill_between(x_grid, density, alpha=0.4, color=colour)

    ax.set_xlabel(column_label)
    ax.set_ylabel("Density")
    ax.set_title(column_label)
    ax.legend(title="condition", loc="upper left", bbox_to_anchor=(1, 1))

    fig.savefig(output_dir / f"{column}_ridge.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def main() -> None:
    args = parse_args()

    if not args.boxplot_dir and not args.ridge_dir:
        print(
            "No output directories specified — nothing to do. "
            "Pass --boxplot-dir and/or --ridge-dir."
        )
        return

    df = pd.read_csv(args.input_csv)
    conditions = df["condition"].unique().tolist()
    target_columns = [c for c in df.columns if c.endswith(TARGET_SUFFIX)]

    boxplot_dir = Path(args.boxplot_dir) if args.boxplot_dir else None
    ridge_dir = Path(args.ridge_dir) if args.ridge_dir else None

    if boxplot_dir:
        boxplot_dir.mkdir(parents=True, exist_ok=True)
    if ridge_dir:
        ridge_dir.mkdir(parents=True, exist_ok=True)

    for column in target_columns:
        print(f"Processing column: {column}")
        column_label = resolve_column_label(column)

        if boxplot_dir:
            plot_boxplot(
                df, column, column_label, conditions,
                boxplot_dir, args.add_outlier_labels,
            )

        if ridge_dir:
            plot_ridge(df, column, column_label, conditions, ridge_dir)


if __name__ == "__main__":
    main()