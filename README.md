CRISPR Gene Editing Analysis Pipeline
An automated Snakemake pipeline for analyzing CRISPR-Cas9 and base editor outcomes from amplicon sequencing data. Quantifies on-target editing efficiency, indel profiles, and allele-specific editing discrimination using CRISPResso2.
Overview
This pipeline processes raw amplicon sequencing data from CRISPR gene editing experiments and produces publication-quality analysis results. It compares three different editing approaches:

Cas9 Nuclease (NHEJ) — Standard double-strand break editing with insertion/deletion outcomes
Adenine Base Editor (EMX1 locus) — Precise A→G single-base conversion without double-strand breaks
Allele-Specific Editing (Rho gene) — Targeted editing of the P23H mutant allele vs wild-type

Key Results
ExperimentTargetEditing EfficiencyKey FindingCas9 NHEJReference amplicon76.6%High efficiency, characteristic indel profileBase EditorEMX1 locus29.0%Precise substitution with minimal indelsAllele-SpecificP23H mutant55.3%2.5x preference over wild-type (22.4%)
Pipeline Architecture
Raw FASTQ → FastQC/MultiQC → CRISPResso2 → Python Visualization
                                    ↓
                            Editing quantification
                            Indel profiling
                            Allele discrimination
All steps automated via Snakemake — full analysis runs with:
bashsnakemake --cores 4
Repository Structure
crispr-off-target-pipeline/
├── README.md                     # This file
├── Snakefile                     # Automated pipeline definition
├── environment.yml               # Reproducible conda environment
├── config/
│   └── config.yaml               # All parameters and file paths
├── data/
│   └── raw/                      # Input FASTQ files
├── results/
│   ├── qc/                       # FastQC and MultiQC reports
│   ├── crispresso/               # CRISPResso2 output per experiment
│   └── figures/                  # Publication-quality visualizations
├── scripts/
│   └── visualize_results.py      # Custom Python visualization code
├── notebooks/
│   └── analysis_report.ipynb     # Detailed Jupyter analysis report
└── docs/
    └── methods.md                # Methods documentation
Quick Start
1. Clone the repository
bashgit clone https://github.com/yourusername/crispr-off-target-pipeline.git
cd crispr-off-target-pipeline
2. Create the conda environment
bashconda env create -f environment.yml
conda activate crispr-offTarget
3. Run the full pipeline
bashsnakemake --cores 4
This will execute all steps: quality control, CRISPResso2 analysis on all three datasets, and visualization generation.
Tools Used
ToolVersionPurposeFastQC0.12.1Raw read quality assessmentMultiQC1.14Aggregated QC reportingCRISPResso22.2.9CRISPR editing quantificationSnakemake7.32.4Workflow automationPython3.10Custom visualization and analysismatplotlib3.7.3Publication-quality figurespandas—Data parsing and manipulation
Data Source
Example datasets from the CRISPResso2 project (Clement et al., 2019, Nature Biotechnology). These contain real amplicon sequencing data from CRISPR editing experiments:

NHEJ dataset: Paired-end Cas9 editing experiment (25,000 reads)
Base editor dataset: EMX1 locus base editing experiment (25,000 reads)
Allele-specific dataset: Rho gene P23H/WT allele discrimination (25,000 reads)

Relevance to Therapeutic Gene Editing
This pipeline demonstrates skills directly applicable to evaluating LNP-delivered gene editing therapies:

Quantification of on-target and off-target editing outcomes
Comparison of editing efficacy across different experimental approaches
Automated, reproducible analysis workflows for high-throughput data
Clear visualization of results for cross-functional communication

Author
Priya D — Graduate Bionformatics student
Northeastern University 

References

Clement, K. et al. (2019). CRISPResso2 provides accurate and rapid genome editing sequence analysis. Nature Biotechnology, 37, 224-226.
Tsai, S.Q. et al. (2015). GUIDE-seq enables genome-wide profiling of off-target cleavage by CRISPR-Cas nucleases. Nature Biotechnology, 33, 187-197.
Rothgangl, T. et al. (2021). In vivo adenine base editing of PCSK9 in macaques reduces LDL cholesterol levels. Nature Biotechnology, 39, 949-957.
