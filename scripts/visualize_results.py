#!/usr/bin/env python
"""
Visualization script for CRISPR gene editing analysis pipeline.
Generates publication-quality figures comparing editing outcomes
across Cas9, base editor, and allele-specific experiments.

Author: Priya D
Project: CRISPR Editing Analysis Pipeline
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

# ============================================================
# GLOBAL STYLE — Professional bioinformatics look
# ============================================================
mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 12,
    'axes.titlesize': 16,
    'axes.titleweight': 'bold',
    'axes.labelsize': 13,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'legend.fontsize': 11,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.linewidth': 1.2,
    'xtick.major.width': 1.0,
    'ytick.major.width': 1.0,
})

# ============================================================
# COLOR PALETTE — Coordinated teal/blue family
# ============================================================
COLORS = {
    'primary':    '#1B4F72',   # Deep navy
    'secondary':  '#2E86C1',   # Medium blue
    'tertiary':   '#85C1E9',   # Light blue
    'accent':     '#D4E6F1',   # Very light blue
    'highlight':  '#E74C3C',   # Red for emphasis
    'insert':     '#2E86C1',   # Blue — insertions
    'delete':     '#1B4F72',   # Dark navy — deletions
    'subst':      '#85C1E9',   # Light blue — substitutions
    'modified':   '#1B4F72',   # Dark — modified
    'unmodified': '#D4E6F1',   # Light — unmodified
    'p23h':       '#1B4F72',   # Dark — mutant allele
    'wt':         '#85C1E9',   # Light — wild-type allele
}

# ============================================================
# CONFIGURATION
# ============================================================
CRISPRESSO_DIR = "results/crispresso"
FIGURES_DIR = "results/figures"
os.makedirs(FIGURES_DIR, exist_ok=True)

EXPERIMENTS = {
    "nhej": {
        "label": "Cas9 (NHEJ)",
        "path": os.path.join(CRISPRESSO_DIR, "CRISPResso_on_nhej")
    },
    "base_editor": {
        "label": "Base Editor (EMX1)",
        "path": os.path.join(CRISPRESSO_DIR, "CRISPResso_on_base_editor")
    },
    "allele_specific": {
        "label": "Allele-Specific (Rho)",
        "path": os.path.join(CRISPRESSO_DIR, "CRISPResso_on_allele_specific")
    }
}

# ============================================================
# HELPER: Parse CRISPResso2 quantification file
# ============================================================
def parse_quantification(filepath):
    df = pd.read_csv(filepath, sep="\t")
    return df

# ============================================================
# FIGURE 1: Editing Efficiency Comparison
# ============================================================
def plot_editing_efficiency():
    labels = []
    modified_pcts = []
    bar_colors = [COLORS['primary'], COLORS['secondary'], COLORS['tertiary'], COLORS['accent']]

    for key, exp in EXPERIMENTS.items():
        qf = os.path.join(exp["path"], "CRISPResso_quantification_of_editing_frequency.txt")
        if os.path.exists(qf):
            df = parse_quantification(qf)
            for _, row in df.iterrows():
                amp = row["Amplicon"] if "Amplicon" in df.columns else key
                labels.append(f"{exp['label']}\n({amp})")
                modified_pcts.append(row["Modified%"])

    fig, ax = plt.subplots(figsize=(10, 6))
    bars = ax.bar(range(len(labels)), modified_pcts,
                  color=bar_colors[:len(labels)], edgecolor='white', width=0.55, linewidth=1.5)

    for bar, val in zip(bars, modified_pcts):
        ax.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 1.5,
                f'{val:.1f}%', ha='center', va='bottom', fontweight='bold', fontsize=12,
                color=COLORS['primary'])

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels)
    ax.set_ylabel("Editing Efficiency (%)")
    ax.set_title("Figure 1 — Editing Efficiency Across Experimental Approaches")
    ax.set_ylim(0, 100)
    ax.yaxis.grid(True, alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "fig1_editing_efficiency.png"))
    plt.savefig(os.path.join(FIGURES_DIR, "fig1_editing_efficiency.svg"))
    print("Figure 1 saved: fig1_editing_efficiency.png")
    plt.close()

# ============================================================
# FIGURE 2: Indel Size Distribution (NHEJ)
# ============================================================
def plot_indel_distribution():
    indel_file = os.path.join(EXPERIMENTS["nhej"]["path"], "Indel_histogram.txt")
    if not os.path.exists(indel_file):
        print("Indel histogram not found, skipping Figure 2")
        return

    df = pd.read_csv(indel_file, sep="\t")
    sizes = df.iloc[:, 0]
    counts = df.iloc[:, 1]

    # Color negative (deletions) darker, positive (insertions) lighter
    colors = [COLORS['delete'] if s < 0 else COLORS['insert'] if s > 0 else '#95a5a6' for s in sizes]

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.bar(sizes, counts, color=colors, edgecolor='none', width=0.8)
    ax.axvline(x=0, color=COLORS['primary'], linestyle='-', alpha=0.5, linewidth=1.5)
    ax.set_xlabel("Indel Size (bp)")
    ax.set_ylabel("Number of Reads")
    ax.set_title("Figure 2 — Indel Size Distribution (Cas9 NHEJ)")
    ax.yaxis.grid(True, alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)

    # Add annotation
    ax.annotate('← Deletions', xy=(-30, ax.get_ylim()[1]*0.85),
                fontsize=11, color=COLORS['delete'], fontweight='bold')
    ax.annotate('Insertions →', xy=(10, ax.get_ylim()[1]*0.85),
                fontsize=11, color=COLORS['insert'], fontweight='bold')

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "fig2_indel_distribution.png"))
    print("Figure 2 saved: fig2_indel_distribution.png")
    plt.close()

# ============================================================
# FIGURE 3: Modification Types (FIXED — uses Only columns)
# ============================================================
def plot_modification_types():
    labels = []
    only_ins = []
    only_del = []
    only_sub = []
    mixed = []

    for key, exp in EXPERIMENTS.items():
        qf = os.path.join(exp["path"], "CRISPResso_quantification_of_editing_frequency.txt")
        if os.path.exists(qf):
            df = parse_quantification(qf)
            for _, row in df.iterrows():
                amp = row["Amplicon"] if "Amplicon" in df.columns else key
                total_mod = row["Modified"]
                if total_mod > 0:
                    labels.append(f"{exp['label']}\n({amp})")
                    oi = row.get("Only Insertions", 0) / total_mod * 100
                    od = row.get("Only Deletions", 0) / total_mod * 100
                    os_ = row.get("Only Substitutions", 0) / total_mod * 100
                    remainder = 100 - oi - od - os_
                    only_ins.append(oi)
                    only_del.append(od)
                    only_sub.append(os_)
                    mixed.append(max(0, remainder))

    fig, ax = plt.subplots(figsize=(10, 6))
    x = range(len(labels))
    w = 0.55

    b1 = ax.bar(x, only_ins, w, label="Insertions Only", color=COLORS['insert'], edgecolor='white', linewidth=1)
    b2 = ax.bar(x, only_del, w, bottom=only_ins, label="Deletions Only", color=COLORS['delete'], edgecolor='white', linewidth=1)
    bot2 = [i+d for i,d in zip(only_ins, only_del)]
    b3 = ax.bar(x, only_sub, w, bottom=bot2, label="Substitutions Only", color=COLORS['subst'], edgecolor='white', linewidth=1)
    bot3 = [b+s for b,s in zip(bot2, only_sub)]
    b4 = ax.bar(x, mixed, w, bottom=bot3, label="Mixed/Complex", color=COLORS['accent'], edgecolor='white', linewidth=1)

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Proportion of Modified Reads (%)")
    ax.set_title("Figure 3 — Types of Editing Modifications")
    ax.set_ylim(0, 105)
    ax.legend(loc='upper right', framealpha=0.9, edgecolor='none')
    ax.yaxis.grid(True, alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "fig3_modification_types.png"))
    print("Figure 3 saved: fig3_modification_types.png")
    plt.close()

# ============================================================
# FIGURE 4: Allele Discrimination
# ============================================================
def plot_allele_discrimination():
    qf = os.path.join(EXPERIMENTS["allele_specific"]["path"],
                       "CRISPResso_quantification_of_editing_frequency.txt")
    if not os.path.exists(qf):
        print("Allele-specific data not found, skipping Figure 4")
        return

    df = parse_quantification(qf)

    fig, ax = plt.subplots(figsize=(7, 6))
    colors = [COLORS['p23h'], COLORS['wt']]
    bars = ax.bar(df["Amplicon"], df["Modified%"], color=colors, edgecolor='white', width=0.45, linewidth=1.5)

    for bar, val in zip(bars, df["Modified%"]):
        ax.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 1.5,
                f'{val:.1f}%', ha='center', va='bottom', fontweight='bold', fontsize=13,
                color=COLORS['primary'])

    ax.set_ylabel("Editing Efficiency (%)")
    ax.set_title("Figure 4 — Allele-Specific Editing Discrimination")
    ax.set_ylim(0, 75)
    ax.yaxis.grid(True, alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)

    # Add fold-change annotation
    p23h_val = df[df["Amplicon"]=="P23H"]["Modified%"].values[0]
    wt_val = df[df["Amplicon"]=="WT"]["Modified%"].values[0]
    fold = p23h_val / wt_val
    ax.annotate(f'{fold:.1f}x higher editing\non target allele',
                xy=(0.5, 60), fontsize=11, ha='center',
                color=COLORS['primary'], fontstyle='italic')

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "fig4_allele_discrimination.png"))
    print("Figure 4 saved: fig4_allele_discrimination.png")
    plt.close()

# ============================================================
# FIGURE 5: Summary Dashboard
# ============================================================
def plot_summary_dashboard():
    """Combined summary of all key metrics."""
    fig, axes = plt.subplots(1, 3, figsize=(24, 7), gridspec_kw={'width_ratios': [1, 1, 1]})
    fig.subplots_adjust(wspace=0.5)

    # --- Panel A: Read Alignment Rates ---
    exp_labels = []
    aligned_pcts = []
    panel_colors = [COLORS['primary'], COLORS['secondary'], COLORS['tertiary']]
    for key, exp in EXPERIMENTS.items():
        qf = os.path.join(exp["path"], "CRISPResso_quantification_of_editing_frequency.txt")
        if os.path.exists(qf):
            df = parse_quantification(qf)
            total_in = df["Reads_in_input"].iloc[0]
            total_al = df["Reads_aligned"].iloc[0]
            exp_labels.append(exp["label"])
            aligned_pcts.append(total_al / total_in * 100)

    bars_a = axes[0].barh(exp_labels, aligned_pcts, color=panel_colors, edgecolor='white', height=0.5, linewidth=1.5)
    for bar, val in zip(bars_a, aligned_pcts):
        axes[0].text(bar.get_width() + 0.3, bar.get_y() + bar.get_height()/2.,
                     f'{val:.1f}%', ha='left', va='center', fontweight='bold', color=COLORS['primary'], fontsize=12)
    axes[0].set_xlabel("Reads Aligned (%)")
    axes[0].set_title("A) Read Alignment Rate", fontweight="bold", pad=15)
    axes[0].set_xlim(0, 110)
    axes[0].xaxis.grid(True, alpha=0.3, linestyle='--')
    axes[0].set_axisbelow(True)

    # --- Panel B: Modified vs Unmodified Pie (NHEJ) ---
    nhej_f = os.path.join(EXPERIMENTS["nhej"]["path"],
                           "CRISPResso_quantification_of_editing_frequency.txt")
    df_nhej = parse_quantification(nhej_f)
    sizes = [df_nhej["Modified%"].iloc[0], df_nhej["Unmodified%"].iloc[0]]
    wedges, texts, autotexts = axes[1].pie(
        sizes, labels=["Modified", "Unmodified"],
        colors=[COLORS['modified'], COLORS['unmodified']],
        autopct='%1.1f%%', startangle=90,
        textprops={'fontsize': 12},
        wedgeprops={'edgecolor': 'white', 'linewidth': 2},
        pctdistance=0.55
    )
    autotexts[0].set_color('white')
    autotexts[0].set_fontweight('bold')
    autotexts[1].set_color(COLORS['primary'])
    autotexts[1].set_fontweight('bold')
    axes[1].set_title("B) Cas9 NHEJ Editing Outcome", fontweight="bold", pad=15)

    # --- Panel C: Allele Discrimination ---
    al_f = os.path.join(EXPERIMENTS["allele_specific"]["path"],
                         "CRISPResso_quantification_of_editing_frequency.txt")
    df_al = parse_quantification(al_f)
    bars_c = axes[2].bar(df_al["Amplicon"], df_al["Modified%"],
                         color=[COLORS['p23h'], COLORS['wt']], edgecolor='white', width=0.45, linewidth=1.5)
    for bar, val in zip(bars_c, df_al["Modified%"]):
        axes[2].text(bar.get_x() + bar.get_width()/2., bar.get_height() + 1.5,
                     f'{val:.1f}%', ha='center', va='bottom', fontweight='bold',
                     fontsize=12, color=COLORS['primary'])
    axes[2].set_ylabel("Editing Efficiency (%)")
    axes[2].set_title("C) Allele Discrimination", fontweight="bold", pad=15)
    axes[2].set_ylim(0, 75)
    axes[2].yaxis.grid(True, alpha=0.3, linestyle='--')
    axes[2].set_axisbelow(True)

    fig.suptitle("CRISPR Editing Analysis — Summary Dashboard",
                 fontsize=20, fontweight="bold", color=COLORS['primary'], y=1.0)
    plt.savefig(os.path.join(FIGURES_DIR, "fig5_summary_dashboard.png"), bbox_inches='tight')
    print("Figure 5 saved: fig5_summary_dashboard.png")
    plt.close()

# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    print("=" * 60)
    print("CRISPR Editing Analysis — Generating Visualizations")
    print("=" * 60)

    plot_editing_efficiency()
    plot_indel_distribution()
    plot_modification_types()
    plot_allele_discrimination()
    plot_summary_dashboard()

    print("=" * 60)
    print("All figures saved to results/figures/")
    print("=" * 60)