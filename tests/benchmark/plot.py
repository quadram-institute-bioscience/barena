#!/usr/bin/env python3
"""Plot benchmark CSVs from tests/benchmark/{Linux,Darwin}/*.csv."""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCENARIO_LABELS = {
    "children": "With children propagation",
    "no-children": "Without children propagation",
}
OS_COLORS = {"Linux": "#4878d0", "Darwin": "#ee854a"}
OS_DISPLAY = {"Linux": "Linux", "Darwin": "macOS"}


def load_data(base_dir: Path) -> dict:
    data = {}
    for os_name in ("Linux", "Darwin"):
        for scenario in ("children", "no-children"):
            path = base_dir / os_name / f"benchmark-{scenario}.csv"
            if path.exists():
                df = pd.read_csv(path)
                df["mean_ms"] = df["mean"] * 1000
                df["stddev_ms"] = df["stddev"] * 1000
                data[(os_name, scenario)] = df
    return data


def plot(data: dict, out_path: Path) -> None:
    os_names = [o for o in ("Linux", "Darwin") if any(o == k[0] for k in data)]
    scenarios = [s for s in ("children", "no-children") if any(s == k[1] for k in data)]

    # Stable command order: union across all frames
    all_cmds: list[str] = []
    for df in data.values():
        for cmd in df["command"]:
            if cmd not in all_cmds:
                all_cmds.append(cmd)

    fig, axes = plt.subplots(1, len(scenarios), figsize=(7 * len(scenarios), 6), sharey=False)
    if len(scenarios) == 1:
        axes = [axes]

    fig.suptitle("Kraken Filter — Benchmark Comparison", fontsize=14, fontweight="bold")

    bar_w = 0.35
    n_os = len(os_names)
    offsets = np.linspace(-(n_os - 1) * bar_w / 2, (n_os - 1) * bar_w / 2, n_os)

    for ax, scenario in zip(axes, scenarios):
        # Commands present in this scenario (preserve order)
        cmds = [c for c in all_cmds
                if any(c in data[k]["command"].values for k in data if k[1] == scenario)]
        x = np.arange(len(cmds))

        max_top = 0.0
        for i, os_name in enumerate(os_names):
            key = (os_name, scenario)
            if key not in data:
                continue
            df = data[key].set_index("command")

            means = np.array([df.loc[c, "mean_ms"] if c in df.index else np.nan for c in cmds])
            stds  = np.array([df.loc[c, "stddev_ms"] if c in df.index else np.nan for c in cmds])

            bars = ax.bar(
                x + offsets[i], means, bar_w,
                label=OS_DISPLAY[os_name],
                color=OS_COLORS[os_name],
                yerr=stds,
                capsize=4,
                error_kw={"elinewidth": 1.5, "ecolor": "#222", "alpha": 0.75},
                alpha=0.87,
                zorder=3,
            )

            for bar, val, std in zip(bars, means, stds):
                if np.isnan(val):
                    continue
                top = val + (std if not np.isnan(std) else 0)
                max_top = max(max_top, top)
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    top + max_top * 0.01,
                    f"{val:.0f}",
                    ha="center", va="bottom", fontsize=7.5, color="#333",
                )

        ax.set_xticks(x)
        ax.set_xticklabels(cmds, rotation=28, ha="right", fontsize=9)
        ax.set_ylabel("Time (ms)", fontsize=10)
        ax.set_title(SCENARIO_LABELS[scenario], fontsize=11, pad=8)
        ax.yaxis.grid(True, linestyle="--", alpha=0.45, zorder=0)
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend(fontsize=9, framealpha=0.7)
        ax.set_ylim(bottom=0)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {out_path}")


def main() -> None:
    base_dir = Path(__file__).parent
    data = load_data(base_dir)
    if not data:
        print("No CSV files found.", file=sys.stderr)
        sys.exit(1)
    out_path = base_dir / "benchmark_plot.png"
    plot(data, out_path)


if __name__ == "__main__":
    main()
