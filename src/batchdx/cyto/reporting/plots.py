"""Inline SVG figures for the offline report.

Figures are rendered to SVG strings and embedded directly in the HTML, so the report
is a single file that opens without a server and without any external asset.
"""

from __future__ import annotations

import io

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def _to_svg(fig, salt: str) -> str:
    buffer = io.StringIO()
    # Matplotlib otherwise writes the current time and random element IDs into SVG.
    # Use a stable, figure-specific salt so repeated reports have identical bytes and
    # separate inline figures do not share IDs in the same HTML document.
    with matplotlib.rc_context({"svg.hashsalt": salt}):
        fig.savefig(buffer, format="svg", bbox_inches="tight", metadata={"Date": None})
    plt.close(fig)
    svg = buffer.getvalue()
    return svg[svg.index("<svg") :]


def density_by_batch(
    values_by_batch: dict[str, np.ndarray],
    title: str,
    xlabel: str,
    shared_limits: tuple[float, float] | None = None,
) -> str:
    """Per-batch kernel-free density (histogram outline) on comparable axes."""

    fig, ax = plt.subplots(figsize=(4.2, 2.4), dpi=110)
    all_values = np.concatenate([v for v in values_by_batch.values() if len(v)]) \
        if values_by_batch else np.array([0.0])
    lo, hi = shared_limits or (float(all_values.min()), float(all_values.max()))
    if hi <= lo:
        hi = lo + 1.0
    bins = np.linspace(lo, hi, 60)
    for batch, values in sorted(values_by_batch.items()):
        if len(values) == 0:
            continue
        counts, edges = np.histogram(values, bins=bins, density=True)
        centres = 0.5 * (edges[:-1] + edges[1:])
        ax.plot(centres, counts, label=str(batch), linewidth=1.2)
    ax.set_title(title, fontsize=9)
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel("density", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.legend(fontsize=6, frameon=False)
    return _to_svg(fig, f"density:{title}")


def before_after(
    before: dict[str, np.ndarray], after: dict[str, np.ndarray], title: str, xlabel: str
) -> str:
    """Before and after panels on identical axes, so the comparison is honest."""

    pooled = [v for v in list(before.values()) + list(after.values()) if len(v)]
    all_values = np.concatenate(pooled) if pooled else np.array([0.0, 1.0])
    lo, hi = float(all_values.min()), float(all_values.max())
    if hi <= lo:
        hi = lo + 1.0
    bins = np.linspace(lo, hi, 60)

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.4), dpi=110, sharex=True, sharey=True)
    for ax, data, label in zip(axes, (before, after), ("uncorrected", "candidate"), strict=True):
        for batch, values in sorted(data.items()):
            if len(values) == 0:
                continue
            counts, edges = np.histogram(values, bins=bins, density=True)
            ax.plot(0.5 * (edges[:-1] + edges[1:]), counts, linewidth=1.1, label=str(batch))
        ax.set_title(f"{title} — {label}", fontsize=9)
        ax.set_xlabel(xlabel, fontsize=8)
        ax.tick_params(labelsize=7)
    axes[0].set_ylabel("density", fontsize=8)
    axes[0].legend(fontsize=6, frameon=False)
    return _to_svg(fig, f"before_after:{title}")


def frequency_bars(freq: dict[str, dict[str, float]], batch_of: dict[str, str]) -> str:
    """Population frequencies per sample, grouped by batch."""

    samples = sorted(freq)
    populations = sorted({p for row in freq.values() for p in row})
    fig, ax = plt.subplots(figsize=(max(4.5, 0.32 * len(samples)), 2.6), dpi=110)
    bottom = np.zeros(len(samples))
    for pop in populations:
        heights = np.array([freq[s].get(pop, 0.0) for s in samples])
        ax.bar(range(len(samples)), heights, bottom=bottom, label=pop, width=0.8)
        bottom += heights
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels([f"{s}\n{batch_of.get(s, '')}" for s in samples], rotation=90, fontsize=5)
    ax.set_ylabel("fraction of evaluation cells", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.legend(fontsize=6, frameon=False, ncol=3)
    return _to_svg(fig, "frequency_bars")


def drift_lines(series: dict[str, np.ndarray], title: str) -> str:
    """Block medians against acquisition order for a handful of sample-marker pairs."""

    fig, ax = plt.subplots(figsize=(4.4, 2.4), dpi=110)
    for label, values in sorted(series.items()):
        ax.plot(np.linspace(0, 1, len(values)), values, linewidth=1.1, label=label)
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("acquisition order (scaled)", fontsize=8)
    ax.set_ylabel("block median", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.legend(fontsize=6, frameon=False)
    return _to_svg(fig, f"drift:{title}")
