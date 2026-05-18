#!/usr/bin/env python3
"""Plot the composition of ViroForge synthetic data (R1 vs R2) as a grouped barplot."""

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np

# Read counts from ViroForge wastewater virome dataset
categories = ["Viral", "rRNA", "Host DNA", "Reagent\nBacteria", "ERV\n(Endogenous)", "ERV\n(Exogenous)", "PhiX"]
r1_counts = np.array([518300, 73085, 51823, 10365, 5180, 2590, 518])
r2_counts = np.array([518300, 73085, 51823, 10365, 5180, 2590, 518])

total_r1 = r1_counts.sum()
total_r2 = r2_counts.sum()

r1_pct = r1_counts / total_r1 * 100
r2_pct = r2_counts / total_r2 * 100

# Colors matching typical virome paper palettes
colors_r1 = "#4C72B0"
colors_r2 = "#DD8452"

x = np.arange(len(categories))
width = 0.35

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5), gridspec_kw={"width_ratios": [1, 1]})

# --- Left panel: all categories ---
bars1 = ax1.bar(x - width / 2, r1_pct, width, label="R1", color=colors_r1, edgecolor="white", linewidth=0.5)
bars2 = ax1.bar(x + width / 2, r2_pct, width, label="R2", color=colors_r2, edgecolor="white", linewidth=0.5)

ax1.set_ylabel("Proportion of reads (%)", fontsize=12)
ax1.set_title("ViroForge Synthetic Data Composition\n(Wastewater Virome, Heavy Contamination)", fontsize=13, fontweight="bold")
ax1.set_xticks(x)
ax1.set_xticklabels(categories, fontsize=10)
ax1.legend(fontsize=11)
ax1.yaxis.set_major_formatter(mtick.PercentFormatter(decimals=0))
ax1.set_ylim(0, 85)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)

# Add percentage labels on bars
for bar, pct in zip(bars1, r1_pct):
    if pct > 2:
        ax1.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.8,
                 f"{pct:.1f}%", ha="center", va="bottom", fontsize=8, color=colors_r1)
for bar, pct in zip(bars2, r2_pct):
    if pct > 2:
        ax1.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.8,
                 f"{pct:.1f}%", ha="center", va="bottom", fontsize=8, color=colors_r2)

# --- Right panel: zoom into non-viral (contaminant) categories ---
contam_idx = list(range(1, len(categories)))  # skip viral
contam_cats = [categories[i] for i in contam_idx]
contam_r1 = r1_pct[contam_idx]
contam_r2 = r2_pct[contam_idx]

x2 = np.arange(len(contam_cats))
bars3 = ax2.bar(x2 - width / 2, contam_r1, width, label="R1", color=colors_r1, edgecolor="white", linewidth=0.5)
bars4 = ax2.bar(x2 + width / 2, contam_r2, width, label="R2", color=colors_r2, edgecolor="white", linewidth=0.5)

ax2.set_ylabel("Proportion of reads (%)", fontsize=12)
ax2.set_title("Non-Viral Contaminant Detail\n(zoomed in)", fontsize=13, fontweight="bold")
ax2.set_xticks(x2)
ax2.set_xticklabels(contam_cats, fontsize=10)
ax2.legend(fontsize=11)
ax2.yaxis.set_major_formatter(mtick.PercentFormatter(decimals=1))
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)

# Add percentage + count labels
for bar, pct, cnt in zip(bars3, contam_r1, r1_counts[contam_idx]):
    ax2.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.15,
             f"{pct:.2f}%\n({cnt:,})", ha="center", va="bottom", fontsize=7.5, color=colors_r1)
for bar, pct, cnt in zip(bars4, contam_r2, r2_counts[contam_idx]):
    ax2.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.15,
             f"{pct:.2f}%\n({cnt:,})", ha="center", va="bottom", fontsize=7.5, color=colors_r2)

# Add total read count annotation
fig.text(0.5, 0.01,
         f"Total reads: R1 = {total_r1:,}  |  R2 = {total_r2:,}  |  Dataset: Wastewater Virome (Urban Treatment Plant, 10× coverage, NovaSeq, heavy contamination)",
         ha="center", fontsize=9, color="gray")

plt.tight_layout(rect=[0, 0.04, 1, 1])
outpath = "test_data/host_depletion_test/results/synthetic_composition.png"
plt.savefig(outpath, dpi=200, bbox_inches="tight")
print(f"Saved to {outpath}")
plt.close()
