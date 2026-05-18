#!/usr/bin/env python3
"""Plot IBD gut virome composition before and after virome-qc host removal."""

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np

# === INPUT (before QC) ===
categories = ["Viral", "rRNA", "Host DNA", "PhiX", "Reagent\nBacteria"]
input_counts = np.array([482123, 6710, 1693, 287, 108])
input_total = input_counts.sum()
input_pct = input_counts / input_total * 100

# === CLEAN OUTPUT (after QC) ===
clean_counts = np.array([325531, 326, 1, 0, 108])
clean_total = clean_counts.sum()
clean_pct = clean_counts / clean_total * 100

colors_before = "#4C72B0"
colors_after = "#55A868"

fig, axes = plt.subplots(2, 2, figsize=(16, 11))

# ============================================================
# TOP LEFT: Before QC — all categories
# ============================================================
ax = axes[0, 0]
x = np.arange(len(categories))
bars = ax.bar(x, input_pct, 0.6, color=colors_before, edgecolor="white", linewidth=0.5)
ax.set_ylabel("Proportion of reads (%)", fontsize=11)
ax.set_title("Before QC (Input)\nAll Categories", fontsize=13, fontweight="bold")
ax.set_xticks(x)
ax.set_xticklabels(categories, fontsize=10)
ax.yaxis.set_major_formatter(mtick.PercentFormatter(decimals=0))
ax.set_ylim(0, 105)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
for bar, pct, cnt in zip(bars, input_pct, input_counts):
    if pct > 0.5:
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                f"{pct:.2f}%\n({cnt:,})", ha="center", va="bottom", fontsize=8.5, color=colors_before)

# ============================================================
# TOP RIGHT: Before QC — contaminant zoom
# ============================================================
ax = axes[0, 1]
contam_idx = list(range(1, len(categories)))
contam_cats = [categories[i] for i in contam_idx]
contam_pct = input_pct[contam_idx]
contam_counts = input_counts[contam_idx]
x2 = np.arange(len(contam_cats))
bars = ax.bar(x2, contam_pct, 0.6, color=colors_before, edgecolor="white", linewidth=0.5)
ax.set_ylabel("Proportion of reads (%)", fontsize=11)
ax.set_title("Before QC (Input)\nContaminant Detail", fontsize=13, fontweight="bold")
ax.set_xticks(x2)
ax.set_xticklabels(contam_cats, fontsize=10)
ax.yaxis.set_major_formatter(mtick.PercentFormatter(decimals=2))
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
for bar, pct, cnt in zip(bars, contam_pct, contam_counts):
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.015,
            f"{pct:.3f}%\n({cnt:,})", ha="center", va="bottom", fontsize=8.5, color=colors_before)

# ============================================================
# BOTTOM LEFT: After QC — all categories
# ============================================================
ax = axes[1, 0]
bars = ax.bar(x, clean_pct, 0.6, color=colors_after, edgecolor="white", linewidth=0.5)
ax.set_ylabel("Proportion of reads (%)", fontsize=11)
ax.set_title("After QC (Clean Output)\nAll Categories", fontsize=13, fontweight="bold")
ax.set_xticks(x)
ax.set_xticklabels(categories, fontsize=10)
ax.yaxis.set_major_formatter(mtick.PercentFormatter(decimals=0))
ax.set_ylim(0, 105)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
for bar, pct, cnt in zip(bars, clean_pct, clean_counts):
    if pct > 0.01:
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                f"{pct:.2f}%\n({cnt:,})", ha="center", va="bottom", fontsize=8.5, color="#3a7a4a")

# ============================================================
# BOTTOM RIGHT: After QC — contaminant zoom
# ============================================================
ax = axes[1, 1]
clean_contam_pct = clean_pct[contam_idx]
clean_contam_counts = clean_counts[contam_idx]
bars = ax.bar(x2, clean_contam_pct, 0.6, color=colors_after, edgecolor="white", linewidth=0.5)
ax.set_ylabel("Proportion of reads (%)", fontsize=11)
ax.set_title("After QC (Clean Output)\nContaminant Detail", fontsize=13, fontweight="bold")
ax.set_xticks(x2)
ax.set_xticklabels(contam_cats, fontsize=10)
ax.yaxis.set_major_formatter(mtick.PercentFormatter(decimals=2))
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
# Use same y-axis scale as the before-QC zoom for direct comparison
ax.set_ylim(axes[0, 1].get_ylim())
for bar, pct, cnt in zip(bars, clean_contam_pct, clean_contam_counts):
    if cnt > 0:
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.015,
                f"{pct:.4f}%\n({cnt:,})", ha="center", va="bottom", fontsize=8.5, color="#3a7a4a")
    else:
        ax.text(bar.get_x() + bar.get_width() / 2, 0.01,
                f"0\nremoved", ha="center", va="bottom", fontsize=8, color="gray")

fig.suptitle("ViroForge IBD Gut Virome: Before vs After virome-qc\n(per read direction, R1 shown)",
             fontsize=15, fontweight="bold", y=1.01)

summary = (
    f"Input: {input_total:,} reads  |  Clean: {clean_total:,} reads  |  "
    f"Host removed: {1693-1:,} ({(1693-1)/1693*100:.1f}%)  |  "
    f"rRNA removed: {6710-326:,} ({(6710-326)/6710*100:.1f}%)  |  "
    f"PhiX removed: {287:,} (100%)  |  "
    f"Reagent bacteria: 0 removed (below detection)"
)
fig.text(0.5, -0.01, summary, ha="center", fontsize=9, color="gray")

plt.tight_layout(rect=[0, 0.02, 1, 0.98])
outpath = "test_data/ibd_gut_virome/before_after_qc.png"
plt.savefig(outpath, dpi=200, bbox_inches="tight")
print(f"Saved to {outpath}")
plt.close()
