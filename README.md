# MethNet Snakemake workflow

This repository has the code and data required to reproduce the figures for the article:

[MethNet: a robust approach to identify regulatory hubs and their distal targets in cancer](https://doi.org/10.1101/2023.07.07.548142)

## Usage

To run the code you need to:

1. If not available in your system, install
   [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
2. Edit the configuration file `config/config.yaml` to match your computer setup
3. Execute the command `snakemake --profile config/config.yaml` from the directory where the `Snakefile` is

This command should generate 3 html files under `./data/figures/` with all the
programmatically generated figures of the article.

The expected running times and other resources for each `rule` is specified in the `resources` field 
and refer to the whole pipeline (all 24 TCGA studies). For testing/development users can specify fewer studies.

## Reproducibility notes

The pipeline was developed using:

- slurm executor on NYU's UltraViolet HPC cluster
- conda version `23.3.1`
- mamba version `1.4.2`
- snakemake version `8.5.3`

## Data Directory Structure

- `annot/`: annotation for genes and CREs
  + `hugo_hg19_pancan.tsv`: gene coordinates from xenahubs
  + `hugo_hg19_pancan_augmented.tsv`: genes with annotations (biotype)
  + `intergenic_cre.txt`: (temporary file)
  + `meth-cluster_annotated.csv`: potential CRE with chromatin annotation data
  + `meth-cluster_loops.csv`: number of loops per potential CRE aggregated from `Bhattacharyya_loops`
  + `meth-cluster_TFBS.csv`: transcription factors binding each potential CRE
  + `meth-cluster_beta_summary.csv`: mean and standard deviation of methylation status for each CpG cluster.
  + `methyl450k_hg19_ilmn12.csv`: annotation for CpG probe based on Illumina's manifest file
  + `wgEncodeBroadHmm_state_count.csv`: ChrommHMM states downloaded from Genome Browser.
- `experiments/`: (**INPUT**) data generated for this study
  + `Bhattacharyya_loops.csv.gz`: loops matrix compiled from [FitHiChIP data](https://doi.org/10.1038/s41467-019-11950-y)
  + `HiCUP_summary_report.tsv`: QC report for promoter-capture HiC
  + `pcChIP_loops.tsv`: loops recovered from promoter-capture HiC
  + `perturb_seq_results.h5`: CellRanger output from the analysis of perturb-seq experiment.
  + `protospacer_umi_thresholds.csv`: CellRanger output, threshold to detect sgRNA
- `figures/`: figures uses in the article (final output)
- `genomic-networks/`: cre-gene adjacency network and CpG-clustering
  + `absolute_clusters_probes.tsv`
  + `gene_cgi_overlaps.csv`
  + `pcgene-cluster_complete_200bp_network_1e6.tsv`
  + `pcgene-probe450_distances_list.rds`
  + `pcgene-probe450_network_1e6.tsv`
- `matrices/`:
  + `methyl450k`: methylation data (beta values) for each TCGA cancer. Rows = CpG probes, Columns = TCGA samples
  + `meth_clusters`: methylation data (beta values) for each TCGA cancer. Rows = TCGA samples, Columns = CpG clusters
  + `rnaseq`: gene expression data for each tCGA cancer. Rows = TCGA Samples, Columns = Genes
- `models_enet/`: intermediate results directory to store the glmnet models for each TCGA cancer
- `pancan/`: matrices and clinical data for all TCGA cancers combined (PanCancer)
- `results/`: main output of the pipeline
  + `models_enet_hyperparams.csv`: hyper-parameters (alpha and lambda) for each elastic-net model
  + `models_enet_params.csv`: coefficients for each elastic-net model
  + `models_enet_samples.csv`: TCGA patients used to fit each elastic-net models for each TCGA cancer.
  + `models_enet_clusters.csv`: CpG clusters used to fit each elastic-net models for each TCGA cancer.
  + `intercepts.csv`: intercept coefficients from elastic-net regression models.
  + `associations.csv`: all associations recovered from elastic-net regression models.
  + `methnet.csv`: aggregated associations across multiple cancers.
  + `cluster_score.csv`: score for intergenic CpG clusters
  + `cre_cox_glance.csv`: statistics of Cox-PH models fitted for each intergenic CpG cluster
  + `cre_cox_params.csv`: coefficients of Cox-PH models fitted for each intergenic CpG cluster
  + `interaction_examples.csv`: methylation and gene expression data to generate figure S2 and S3
  + `perturb-seq_results.csv`: interactions experimentally recovered from the perturb-seq assay
- `xenahubs/`: xenahub clones for GDC (`gdc/`) and PanCancerAtlas (`pancanatlas/`) hubs.
  `scripts/download_xenahubs.sh` will clone both hubs but the pipeline only requires the following files:
  + `gdc/genomicMatrix/TCGA-{cancer}/Xena_Matrices/TCGA-{cancer}.methylation450.tsv` (where `{cancer}` is a wildcard for all the TCGA projects required)
  + `pancanatlas/genomicMatrix/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.tsv`
  + `pancanatlas/clinicalMatrix/Survival_SupplementalTable_S1_20171025_xena_sp.tsv`
  + `pancanatlas/probeMap/hugo_gencode_good_hg19_V24lift37_probemap.tsv hugo_hg19_pancan.tsv`
