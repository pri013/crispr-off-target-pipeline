# ============================================================
# Snakefile — CRISPR Editing Analysis Pipeline
# Automates: QC → CRISPResso2 → Visualization
# ============================================================

import yaml

# Load config
with open("config/config.yaml", "r") as f:
    config = yaml.safe_load(f)

SAMPLES = list(config.keys())

# ============================================================
# RULE ALL — defines the final outputs we want
# ============================================================
rule all:
    input:
        # FastQC reports
        expand("results/qc/{sample}_fastqc.html", sample=["nhej.r1", "nhej.r2", "base_editor", "allele_specific"]),
        # MultiQC report
        "results/qc/multiqc_report.html",
        # CRISPResso2 results
        expand("results/crispresso/CRISPResso_on_{sample}/.done", sample=SAMPLES),
        # Visualization figures
        "results/figures/fig1_editing_efficiency.png"

# ============================================================
# RULE: FastQC — Quality control on raw reads
# ============================================================
rule fastqc:
    input:
        "data/raw/{sample}.fastq.gz"
    output:
        "results/qc/{sample}_fastqc.html"
    shell:
        "fastqc {input} -o results/qc/"

# ============================================================
# RULE: MultiQC — Aggregate all FastQC reports
# ============================================================
rule multiqc:
    input:
        expand("results/qc/{sample}_fastqc.html", sample=["nhej.r1", "nhej.r2", "base_editor", "allele_specific"])
    output:
        "results/qc/multiqc_report.html"
    shell:
        "multiqc results/qc/ -o results/qc/ --force"

# ============================================================
# RULE: CRISPResso2 — NHEJ (paired-end)
# ============================================================
rule crispresso_nhej:
    input:
        r1 = config["nhej"]["r1"],
        r2 = config["nhej"]["r2"]
    output:
        touch("results/crispresso/CRISPResso_on_nhej/.done")
    params:
        amplicon = config["nhej"]["amplicon_seq"],
        name = config["nhej"]["name"]
    shell:
        """
        CRISPResso \
            --fastq_r1 {input.r1} \
            --fastq_r2 {input.r2} \
            --amplicon_seq {params.amplicon} \
            -n {params.name} \
            -o results/crispresso/
        touch {output}
        """

# ============================================================
# RULE: CRISPResso2 — Base Editor (single-end)
# ============================================================
rule crispresso_base_editor:
    input:
        r1 = config["base_editor"]["r1"]
    output:
        touch("results/crispresso/CRISPResso_on_base_editor/.done")
    params:
        amplicon = config["base_editor"]["amplicon_seq"],
        guide = config["base_editor"]["guide_seq"],
        name = config["base_editor"]["name"]
    shell:
        """
        CRISPResso \
            --fastq_r1 {input.r1} \
            --amplicon_seq {params.amplicon} \
            --guide_seq {params.guide} \
            --quantification_window_size 10 \
            --quantification_window_center -10 \
            --base_editor_output \
            -n {params.name} \
            -o results/crispresso/
        touch {output}
        """

# ============================================================
# RULE: CRISPResso2 — Allele-Specific (single-end)
# ============================================================
rule crispresso_allele_specific:
    input:
        r1 = config["allele_specific"]["r1"]
    output:
        touch("results/crispresso/CRISPResso_on_allele_specific/.done")
    params:
        amplicon = config["allele_specific"]["amplicon_seq"],
        guide = config["allele_specific"]["guide_seq"],
        name = config["allele_specific"]["name"]
    shell:
        """
        CRISPResso \
            --fastq_r1 {input.r1} \
            --amplicon_seq {params.amplicon} \
            --amplicon_name P23H,WT \
            --guide_seq {params.guide} \
            -n {params.name} \
            -o results/crispresso/
        touch {output}
        """

# ============================================================
# RULE: Visualization — Generate all figures
# ============================================================
rule visualize:
    input:
        expand("results/crispresso/CRISPResso_on_{sample}/.done", sample=SAMPLES)
    output:
        "results/figures/fig1_editing_efficiency.png"
    shell:
        "python scripts/visualize_results.py"