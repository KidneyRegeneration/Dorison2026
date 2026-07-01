#!/usr/bin/env python3
"""
plot_colocalisation.py — Plot MOC_AinB and MOC_BinA by condition.

Produces a box plot with individual data points overlaid, and a ridge plot,
for each directional Manders coefficient. Samples are grouped by condition
on the x-axis. All data is plotted together — no faceting.

Usage
-----
    python plot_colocalisation.py \
        -i combined_coloc.csv \
        --marker-a NPHS1 \
        --marker-b CD31 \
        --boxplot-a-in-b   boxplots_NPHS1_in_CD31 \
        --boxplot-b-in-a   boxplots_CD31_in_NPHS1 \
        --ridgeplot-a-in-b ridgeplots_NPHS1_in_CD31 \
        --ridgeplot-b-in-a ridgeplots_CD31_in_NPHS1 \
        --group-by condition \
        --hue-by treatment
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

# Source - https://stackoverflow.com/a/3496790
# Posted by Matt Williamson, modified by community. See post 'Timeline' for change history
# Retrieved 2026-04-16, License - CC BY-SA 2.5

try:
    import seaborn as sns
    SEABORN = True
except ImportError:
    SEABORN = False

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _condition_order(df: pd.DataFrame, group_col: str) -> list[str]:
    return sorted(df[group_col].dropna().unique().tolist())


def _is_hue_valid(df: pd.DataFrame, hue_col: str | None) -> bool:
    """
    Check if hue column is valid: column exists, has data, and has > 1 unique value.
    """
    if hue_col is None:
        return False
    if hue_col not in df.columns:
        print(f"  Warning: hue_by column '{hue_col}' not found in data. Ignoring hue.")
        return False
    unique_vals = df[hue_col].dropna().unique()
    if len(unique_vals) <= 1:
        print(f"  Warning: hue_by column '{hue_col}' has only {len(unique_vals)} unique value(s). Ignoring hue.")
        return False
    return True


# ---------------------------------------------------------------------------
# Box plot with scatter overlay
# ---------------------------------------------------------------------------

def save_boxplot(
    df: pd.DataFrame,
    moc_col: str,
    group_col: str,
    title: str,
    output_path: Path,
    jitter_on=True,
    hue_col: str | None = None,
) -> None:
    """
    Single box plot with jittered individual points overlaid.
    One box per condition on the x-axis. If hue_col is provided,
    boxes are split by hue within each condition.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    conditions = _condition_order(df, group_col)
    rng = np.random.default_rng(seed=42)

    # Validate hue column
    use_hue = _is_hue_valid(df, hue_col)
    
    if use_hue:
        hues = sorted(df[hue_col].dropna().unique().tolist())
        hue_colors = plt.cm.Set1(np.linspace(0, 1, len(hues)))
        hue_color_map = {hue: color for hue, color in zip(hues, hue_colors)}
    else:
        hues = None
        hue_color_map = None

    # Prepare data structure
    if use_hue:
        # data_per_group: list of dicts, each with (condition, hue) as key
        data_per_group = {}
        x_pos_map = {}
        x_pos = 0
        
        for cond in conditions:
            for hue in hues:
                cond_hue_data = df.loc[
                    (df[group_col] == cond) & (df[hue_col] == hue),
                    moc_col
                ].dropna().values
                data_per_group[(cond, hue)] = cond_hue_data
                x_pos_map[(cond, hue)] = x_pos
                x_pos += 1
        
        x_ticks = []
        x_tick_labels = []
        x_pos = 0
        for cond in conditions:
            # Find middle position for this condition's hue group
            hue_positions = [x_pos_map[(cond, h)] for h in hues]
            x_ticks.append(np.mean(hue_positions))
            x_tick_labels.append(cond)
            x_pos = max(hue_positions) + 2
    else:
        data_per_condition = [
            df.loc[df[group_col] == cond, moc_col].dropna().values
            for cond in conditions
        ]
        condition_colors = plt.cm.Set2(np.linspace(0, 1, len(conditions)))

    fig, ax = plt.subplots(figsize=(max(5, len(conditions) * 1.6 * (len(hues) if use_hue else 1)), 5))

    if use_hue:
        # Draw boxes and points for each (condition, hue) pair
        for (cond, hue), vals in data_per_group.items():
            x_pos = x_pos_map[(cond, hue)]
            color = hue_color_map[hue]
            
            # Draw box
            bp = ax.boxplot(
                [vals] if len(vals) > 0 else [[]],
                positions=[x_pos],
                patch_artist=True,
                widths=0.6,
                medianprops=dict(color="black", linewidth=2),
                flierprops=dict(marker=""),
                zorder=2,
            )
            
            if len(vals) > 0 and len(bp["boxes"]) > 0:
                bp["boxes"][0].set_facecolor(color)
                bp["boxes"][0].set_alpha(0.45)
            
            # Draw scatter points
            if len(vals) > 0:
                if jitter_on:
                    jitter = rng.uniform(-0.12, 0.12, size=len(vals))
                else:
                    jitter = np.zeros(len(vals))
                ax.scatter(
                    x_pos + jitter,
                    vals,
                    color=color,
                    edgecolors="black",
                    linewidths=0.5,
                    s=40,
                    zorder=3,
                    alpha=0.85,
                )
        
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_tick_labels, rotation=30, ha="right", fontsize=9)
        
        # Add legend for hues
        from matplotlib.patches import Patch
        legend_patches = [
            Patch(facecolor=hue_color_map[hue], alpha=0.45, edgecolor="black", label=hue)
            for hue in hues
        ]
        ax.legend(handles=legend_patches, title=hue_col, fontsize=9, title_fontsize=10)
    else:
        # Original single-color boxplot
        bp = ax.boxplot(
            data_per_condition,
            patch_artist=True,
            medianprops=dict(color="black", linewidth=2),
            flierprops=dict(marker=""),
            zorder=2,
        )

        for patch, colour in zip(bp["boxes"], condition_colors):
            patch.set_facecolor(colour)
            patch.set_alpha(0.45)

        for x_pos, (vals, colour) in enumerate(zip(data_per_condition, condition_colors), start=1):
            if len(vals) == 0:
                continue
            if jitter_on:
                jitter = rng.uniform(-0.18, 0.18, size=len(vals))
            else:
                jitter = np.zeros(len(vals))
            ax.scatter(
                x_pos + jitter,
                vals,
                color=colour,
                edgecolors="black",
                linewidths=0.5,
                s=40,
                zorder=3,
                alpha=0.85,
            )

        ax.set_xticks(range(1, len(conditions) + 1))
        ax.set_xticklabels(conditions, rotation=30, ha="right", fontsize=9)

    ax.set_ylabel("Manders coefficient", fontsize=11)
    ax.set_ylim(-0.05, 1.05)
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.2f"))
    ax.set_title(title, fontsize=12)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    print(f"  Saved boxplot → {output_path}")


def save_swarmplot(
    df: pd.DataFrame,
    moc_col: str,
    group_col: str,
    title: str,
    output_path: Path,
) -> None:
    """
    Single box plot with jittered individual points overlaid.
    One box per condition on the x-axis.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    conditions = _condition_order(df, group_col)
    rng = np.random.default_rng(seed=42)

    data_per_condition = [
        df.loc[df[group_col] == cond, moc_col].dropna().values
        for cond in conditions
    ]

    colours = plt.cm.Set2(np.linspace(0, 1, len(conditions)))

    fig, ax = plt.subplots(figsize=(max(5, len(conditions) * 1.6), 5))

    bp = sns.swarmplot(
        data_per_condition,
        zorder=2,
    )

   
    ax.set_xticks(range(1, len(conditions) + 1))
    ax.set_xticklabels(conditions, rotation=30, ha="right", fontsize=9)
    ax.set_ylabel("Manders coefficient", fontsize=11)
    ax.set_ylim(-0.05, 1.05)
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.2f"))
    ax.set_title(title, fontsize=12)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    print(f"  Saved swarmplot → {output_path}")


# ---------------------------------------------------------------------------
# Ridge plot
# ---------------------------------------------------------------------------

def save_ridgeplot(
    df: pd.DataFrame,
    moc_col: str,
    group_col: str,
    title: str,
    output_path: Path,
) -> None:
    """
    Ridge plot: one row per condition, KDE of MOC values across all samples.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    conditions = _condition_order(df, group_col)
    n = len(conditions)
    colours = plt.cm.Set2(np.linspace(0, 1, n))

    fig, axes = plt.subplots(
        n, 1,
        figsize=(7, max(3, n * 1.2)),
        sharex=True,
    )
    if n == 1:
        axes = [axes]

    for ax, cond, colour in zip(axes, conditions, colours):
        vals = df.loc[df[group_col] == cond, moc_col].dropna().values

        if len(vals) >= 2:
            xs = np.linspace(0, 1, 300)
            kde = gaussian_kde(vals, bw_method="scott")
            ys = kde(xs)
            ax.fill_between(xs, ys, alpha=0.6, color=colour)
            ax.plot(xs, ys, color=colour, linewidth=1.5)
        else:
            # Too few points for KDE — fall back to rug plot
            ax.eventplot(vals, colors=[colour], lineoffsets=0.5, linelengths=0.8)

        ax.set_ylabel(cond, fontsize=8, rotation=0, ha="right", va="center")
        ax.set_yticks([])
        ax.spines[["top", "right", "left"]].set_visible(False)

    axes[-1].set_xlabel("Manders coefficient", fontsize=10)
    axes[-1].set_xlim(0, 1)
    fig.suptitle(title, fontsize=12, y=1.01)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved ridgeplot → {output_path}")


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="plot_colocalisation.py",
        description="Plot directional Manders coefficients by condition.",
    )
    parser.add_argument("-i", "--input", required=True,
                        help="Combined colocalisation CSV")
    parser.add_argument("--col-a-in-b", default="MOC_AinB",
                        help="Column name for A-in-B coefficient. Default: MOC_AinB")
    parser.add_argument("--col-b-in-a", default="MOC_BinA",
                        help="Column name for B-in-A coefficient. Default: MOC_BinA")
    parser.add_argument("--marker-a", required=True,
                        help="Human-readable name for marker A (e.g. NPHS1)")
    parser.add_argument("--marker-b", required=True,
                        help="Human-readable name for marker B (e.g. CD31)")
    parser.add_argument("--boxplot-a-in-b", default="boxplots_AinB",
                        help="Output directory for A-in-B box plot.")
    parser.add_argument("--boxplot-b-in-a", default="boxplots_BinA",
                        help="Output directory for B-in-A box plot.")
    parser.add_argument("--ridgeplot-a-in-b", default="ridgeplots_AinB",
                        help="Output directory for A-in-B ridge plot.")
    parser.add_argument("--ridgeplot-b-in-a", default="ridgeplots_BinA",
                        help="Output directory for B-in-A ridge plot.")
    parser.add_argument("--group-by", default="condition",
                        help="Column to group samples by on the x-axis. Default: condition")
    parser.add_argument("--hue-by", default=None,
                        help="Column to split boxes by hue within each condition. Optional.")
    parser.add_argument("--jitter_off", action="store_false",
                        help="Turn jitter off or on.")
    
    return parser


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = build_parser().parse_args()
    name_a = args.marker_a
    name_b = args.marker_b

    print(f"Reading {args.input} …")
    df = pd.read_csv(args.input)

    required = {args.col_a_in_b, args.col_b_in_a, args.group_by}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Combined CSV is missing expected columns: {missing}")

    title_a_in_b = f"{name_a} colocalising in {name_b}"
    title_b_in_a = f"{name_b} colocalising in {name_a}"

    print(f"Generating {name_a}-in-{name_b} plots …")
    save_boxplot(
        df, args.col_a_in_b, args.group_by, title_a_in_b,
        Path(args.boxplot_a_in_b) / f"boxplot_{name_a}_in_{name_b}.png",
        jitter_on=args.jitter_off,
        hue_col=args.hue_by,
    )
    save_ridgeplot(
        df, args.col_a_in_b, args.group_by, title_a_in_b,
        Path(args.ridgeplot_a_in_b) / f"ridgeplot_{name_a}_in_{name_b}.png",
    )
    
    if SEABORN: 
        save_swarmplot(
            df, args.col_a_in_b, args.group_by, title_a_in_b, Path("swarmplot_a_in_b")
        )

    print(f"Generating {name_b}-in-{name_a} plots …")
    save_boxplot(
        df, args.col_b_in_a, args.group_by, title_b_in_a,
        Path(args.boxplot_b_in_a) / f"boxplot_{name_b}_in_{name_a}.png",
        jitter_on=args.jitter_off,
        hue_col=args.hue_by,
    )
    save_ridgeplot(
        df, args.col_b_in_a, args.group_by, title_b_in_a,
        Path(args.ridgeplot_b_in_a) / f"ridgeplot_{name_b}_in_{name_a}.png",
    )
    if SEABORN: 
        save_swarmplot(
            df, args.col_b_in_a, args.group_by, title_b_in_a, Path("swarmplot_b_in_a")
        )
    print("Done.")


if __name__ == "__main__":
    main()