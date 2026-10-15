"""Reproduce the synthetic layout preview; this is not a package renderer."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Rectangle

OUT = Path(__file__).resolve().parent
LABELS = ["EMD H ↓", "EMD V ↓", "Batch R² ↓", "iLISI ↑", "FlowSOM ↑",
          "cLISI ↑", "Peaks ↑", "Wilcoxon ↑", "Δ Cohen d ↓"]
NAMES = ["Candidate A", "Candidate B", "Identity", "Candidate C (partial)", "Candidate D (failed)"]
# Synthetic already-scaled values; no real correction method is represented.
VALUES = np.array([
    [.87, .80, .91, .84, .88, .79, .82, .76, .81],
    [.69, .63, .72, .71, .74, .85, .83, .80, .86],
    [.18, .21, .26, .19, .30, .91, .90, .94, .92],
    [.93, .94, .93, .96, .94, .91, .92, .94, np.nan],
    [np.nan] * 9,
])
assert VALUES.shape == (len(NAMES), len(LABELS))
complete = np.isfinite(VALUES).all(axis=1)
totals = np.array([row.mean() if ok else np.nan for row, ok in zip(VALUES, complete, strict=True)])
assert list(np.flatnonzero(complete)) == [0, 1, 2]
assert totals[0] > totals[1] > totals[2]
assert np.isnan(totals[3:]).all()  # Missing evidence cannot enter the ranking.

with plt.rc_context({"svg.hashsalt": "batchdx-synthetic-benchmark-preview", "svg.fonttype": "none"}):
    fig, ax = plt.subplots(figsize=(17, 7.6))
    ax.set(xlim=(-2.8, 13.4), ylim=(-1.2, 6.1))
    ax.set_aspect("equal")
    ax.axis("off")
    ax.text(-2.7, 5.8, "Benchmark ranking — layout preview", fontsize=17, weight="bold")
    ax.text(-2.7, 5.35, "SYNTHETIC SCORES • invented candidates • no method performance evidence",
            fontsize=11, color="#9b3300", weight="bold")
    for title, left, right in [("Batch / replicate consistency", -.45, 4.45),
                               ("Biological preservation", 4.55, 8.45),
                               ("Score summaries", 8.7, 13.3)]:
        ax.text((left + right) / 2, 4.85, title, ha="center", fontsize=10, weight="bold")
        ax.plot([left, right], [4.6, 4.6], color="#7d8790", lw=.8)
    ax.text(-2.7, 4.1, "Rank / candidate", fontsize=10, weight="bold")
    for x, label in enumerate(LABELS):
        ax.text(x, 4.1, label, fontsize=9, ha="center")
    summary_x = [9.05, 10.5, 11.95]
    for x, label in zip(summary_x, ["Batch", "Biology", "Overall"], strict=True):
        ax.text(x + .5, 4.1, label, fontsize=9, ha="center")
    for i, (name, row) in enumerate(zip(NAMES, VALUES, strict=True)):
        y = 3.5 - i * .75
        rank = str(1 + np.count_nonzero(totals > totals[i])) if complete[i] else "—"
        ax.text(-2.7, y, f"{rank}   {name}", fontsize=10, va="center")
        for x, value in enumerate(row):
            if np.isfinite(value):
                fill = plt.cm.PRGn(value)
                ax.add_patch(Circle((x, y), .27, facecolor=fill, edgecolor="#d0d4d8", lw=.5))
                ax.text(x, y, f"{value:.2f}", ha="center", va="center", fontsize=9,
                        color="white" if value < .25 or value > .75 else "#152026")
            else:
                ax.add_patch(Rectangle((x - .29, y - .23), .58, .46,
                                      facecolor="#f2f2f2", edgecolor="#a2a2a2", hatch="///", lw=.4))
                ax.text(x, y, "ERR" if i == 4 else "NA", ha="center", va="center",
                        fontsize=8, bbox={"facecolor": "white", "edgecolor": "none", "pad": 1})
        groups = [row[:5], row[5:], row]
        for x, group in zip(summary_x, groups, strict=True):
            if np.isfinite(group).all():
                score = group.mean()
                ax.add_patch(Rectangle((x, y - .18), score, .36, facecolor="#31758d"))
                ax.text(x + 1.03, y, f"{score:.2f}", fontsize=9, va="center")
            else:
                ax.text(x + .5, y, "unranked" if x == summary_x[-1] else "incomplete",
                        fontsize=8, ha="center", va="center", color="#666666")
        ax.plot([-2.7, 13.3], [y - .36, y - .36], color="#d2d7db", lw=.6, ls=":")
    ax.text(-2.7, -.55, "Circle values: higher-is-better scaled scores. Arrows: original metric direction. Grey hatching: unavailable / failed.", fontsize=9)
    ax.text(-2.7, -.92, "Overall = equal mean of all 9 scores; complete candidates only. Synthetic values show layout, not measured ranking or approval.", fontsize=9)
    for suffix in ["png", "svg", "pdf"]:
        fig.savefig(OUT / f"benchmark-layout-preview.{suffix}", dpi=160,
                    bbox_inches="tight", metadata={"Date": None} if suffix == "svg" else {})
    plt.close(fig)
