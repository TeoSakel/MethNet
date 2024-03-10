configfile: 'config/params.yaml'

CANCERS=config['tcga_cancers']
CUTOFF=config['cpg_cluster_cutoff']
RADIUS=config['max_association_distance']
LINKAGE=config['linkage']
CLUSTER_SCHEME='{}_{}bp'.format(LINKAGE, CUTOFF)
CRE_CHUNKS=config['cre_chunks']

localrules: gene_annot, probe_annot, combine_results, download_clinical,
            split_intergenic_cre, combine_survival, get_intergenic_cre,
            download_perturb_seq_results, download_protospacer_umi_thresholds

rule all:
    input:
        "data/figures/figures.html",
        "data/figures/figure_6.html",
        "data/figures/figure_7.html"

rule gene_annot:
    input:
        "data/annot/hugo_hg19_pancan.tsv"
    output:
        "data/annot/hugo_hg19_pancan_augmented.tsv"
    conda:
        "envs/preprocess.yaml"
    shell:
        "Rscript scripts/augment_genes.R {input} {output}"

rule probe_annot:
    output:
        "data/annot/methyl450k_hg19_ilmn12.csv"
    conda:
        "envs/preprocess.yaml"
    shell:
        "Rscript scripts/probe450k_annotation.R {output}"

rule cluster_probes:
    input:
        "data/annot/methyl450k_hg19_ilmn12.csv"
    output:
        "data/genomic-networks/absolute_clusters_probes.tsv"
    conda:
        "envs/preprocess.yaml"
    params:
        window=CUTOFF
    shell:
        "Rscript scripts/cluster_probes.R {input} {params.window} {output}"

rule annotate_clusters:
    input:
        manifest="data/annot/methyl450k_hg19_ilmn12.csv",
        clusters="data/genomic-networks/absolute_clusters_probes.tsv",
        loops="data/experiments/Bhattacharyya_loops.csv.gz"
    conda:
        "envs/preprocess.yaml"
    output:
        "data/annot/wgEncodeBroadHmm_state_count.csv",
        "data/annot/meth-cluster_annotated.csv",
        "data/annot/meth-cluster_TFBS.csv",
        "data/annot/meth-cluster_loops.csv"
    params:
        prefix=CLUSTER_SCHEME
    shell:
        "Rscript scripts/annotate_clusters.R {input.manifest} {input.clusters} {params.prefix} data/annot"

rule distances:
    input:
        genes="data/annot/hugo_hg19_pancan_augmented.tsv",
        probes="data/annot/methyl450k_hg19_ilmn12.csv"
    params:
        prefix="data/genomic-networks/pcgene-probe450"
    output:
        "data/genomic-networks/pcgene-probe450_distances_list.rds"
    conda:
        "envs/preprocess.yaml"
    shell:
        """
        Rscript scripts/gene-probe-distances.R \
            --genetype=protein_coding \
            {input.genes} {input.probes} {params.prefix}
        """

rule networks:
    input:
        genes="data/annot/hugo_hg19_pancan_augmented.tsv",
        distances="data/genomic-networks/pcgene-probe450_distances_list.rds"
    output:
        f"data/genomic-networks/pcgene-probe450_network_{RADIUS}.tsv"
    conda:
        "envs/preprocess.yaml"
    params:
        radius=RADIUS
    shell:
        """
        Rscript scripts/gene-probe-network.R \
            --combine {input.distances} {input.genes} \
            {output} {params.radius}
        """

rule cluster_network:
    input:
        net=f"data/genomic-networks/pcgene-probe450_network_{RADIUS}.tsv",
        clusters="data/genomic-networks/absolute_clusters_probes.tsv"
    output:
        f"data/genomic-networks/pcgene-cluster_{CLUSTER_SCHEME}_network_{RADIUS}.tsv"
    conda:
        "envs/preprocess.yaml"
    params:
        scheme=CLUSTER_SCHEME
    shell:
        "Rscript scripts/gene-cluster-network.R {input.net} {input.clusters} {params.scheme} {output}"

rule gene_cluster_overlap:
    input:
        network=f"data/genomic-networks/pcgene-probe450_network_{RADIUS}.tsv",
        clusters="data/genomic-networks/absolute_clusters_probes.tsv"
    output:
        "data/genomic-networks/gene_cluster_overlaps.csv"
    conda:
        "envs/preprocess.yaml"
    params:
        prefix=CLUSTER_SCHEME
    shell:
        """
        Rscript scripts/gene_cluster_overlap.R \
            {input.network} {input.cluster_file} {params.prefix} {output}
        """

rule gene_cgi_overlap:
    input:
        "data/annot/hugo_hg19_pancan_augmented.tsv"
    output:
        "data/genomic-networks/gene_cgi_overlaps.csv"
    conda:
        "envs/preprocess.yaml"
    shell:
        "Rscript scripts/gene_cgi_overlap.R {input} {output}"

rule download_meth:
    output:
        expand("data/matrices/methyl450k/{cancer}.tsv", cancer=CANCERS)
    conda:
        "envs/download_xena.yaml"
    resources:
        runtime=60*24*2
    shell:
        "./scripts/download_xenahubs.sh met data/matrices/methyl450k"

rule download_rna:
    output:
        "data/pancan/rna.tsv"
    conda:
        "envs/download_xena.yaml"
    resources:
        runtime=60*24*2
    shell:
        "./scripts/download_xenahubs.sh rna data/pancan"

rule get_genes:
    input:
        "data/pancan/rna.tsv"
    output:
        "data/annot/hugo_hg19_pancan.tsv"
    shell:
        "(cd data/annot && ln -s ../xenahubs/pancanatlas/probeMap/hugo_gencode_good_hg19_V24lift37_probemap.tsv hugo_hg19_pancan.tsv)"

rule download_clinical:
    input:
        "data/pancan/rna.tsv"
    output:
        "data/pancan/survival_data.tsv"
    conda:
        "envs/download_xena.yaml"
    shell:
        "./scripts/download_xenahubs.sh clinical data/pancan"

rule split_rna_data:
    input:
        rna="data/pancan/rna.tsv",
        clinical="data/pancan/survival_data.tsv"
    output:
        expand("data/matrices/rnaseq/{cancer}.rds", cancer=CANCERS),
        "data/pancan/rna.rds"
    conda:
        "envs/download_xena.yaml"
    resources:
        runtime=60*4,
        mem_mb=8128
    shell:
        "Rscript scripts/split_pancancer.R {input.rna} {input.clinical} data/matrices/rnaseq"

rule aggregate_meth_clusters:
    input:
        cluster_file="data/genomic-networks/absolute_clusters_probes.tsv",
        probes="data/matrices/methyl450k/{cancer}.tsv"
    output:
        "data/matrices/meth_clusters/{cancer}.rds"
    conda:
        "envs/preprocess.yaml"
    resources:
        mem_mb=lambda wc, input: max(2.5 * input.size_mb, 4096),
        cpus_per_task=3,
        runtime=4*60
    params:
        prefix=CLUSTER_SCHEME
    shell:
        """
        Rscript scripts/aggregate_meth_clusters.R mean \
            {input.probes} {input.cluster_file} {params.prefix} {output}
        """

rule combine_meth_clusters:
    input:
        expand("data/matrices/meth_clusters/{cancer}.rds", cancer=CANCERS)
    output:
        "data/pancan/meth.rds"
    conda:
        "envs/preprocess.yaml"
    resources:
        mem_mb=65536
    shell:
        "Rscript scripts/combine_samples.R {output} {input}"

rule summarize_meth_clusters:
    input:
        expand("data/matrices/meth_clusters/{cancer}.rds", cancer=CANCERS)
    output:
        "data/annot/meth-cluster_beta_summary.csv"
    conda:
        "envs/preprocess.yaml"
    resources:
        mem_mb=32768,
        cpus_per_task=8,
        runtime=12*60
    shell:
        "Rscript scripts/meth_cluster_beta_summary.R {output} {input}"

rule interaction_example:
    input:
        expand("data/matrices/meth_clusters/{cancer}.rds", cancer=CANCERS),
        "data/pancan/rna.rds"
    output:
        "data/results/interaction_examples.csv"
    conda:
        "envs/preprocess.yaml"
    shell:
        """
        Rscript scripts/get_example_data.R \
            data/matrices/meth_clusters \
            data/pancan/rna.rds \
            {output} \
            GSTT1 'chr22:24373322-24373323' \
            IFNG 'chr12:68800138-68800139'
        """

rule fit_enet_models:
    input:
        meth="data/matrices/meth_clusters/{cancer}.rds",
        rna="data/matrices/rnaseq/{cancer}.rds",
        net=f"data/genomic-networks/pcgene-cluster_{CLUSTER_SCHEME}_network_{RADIUS}.tsv"
    output:
        models="data/models_enet/{cancer}_enet-models.rds",
        coefs="data/models_enet/{cancer}_enet-params.csv",
        hyper="data/models_enet/{cancer}_enet-hyperparams.csv",
        clusters="data/models_enet/{cancer}_meth-clusters.txt",
        samples="data/models_enet/{cancer}_samples.txt",
        diagnostics="data/models_enet/{cancer}_enet-diagnostics.pdf"
    params:
        var_min=1
    conda:
        "envs/fit_models.yaml"
    resources:
        runtime=60*36,
        mem_mb=65536,
        cpus_per_task=10
    shell:
        """
        Rscript scripts/fit_enet_models.R \
            --norm --gex_var={params.var_min} \
            {input.net} {input.meth} {input.rna} \
            {output.coefs}
        """

rule combine_results:
    input:
        expand(["data/models_enet/{cancer}_enet-params.csv",
                "data/models_enet/{cancer}_enet-hyperparams.csv"], cancer=CANCERS)
    output:
        "data/results/models_enet_params.csv",
        "data/results/models_enet_hyperparams.csv",
        "data/results/models_enet_samples.csv"
    conda:
        "envs/postprocess.yaml"
    shell:
        "Rscript scripts/combine_enet_models.R data/models_enet data/results/models_enet"

rule process_results:
    input:
        net=f"data/genomic-networks/pcgene-cluster_{CLUSTER_SCHEME}_network_{RADIUS}.tsv",
        hyper="data/results/models_enet_hyperparams.csv",
        coef="data/results/models_enet_params.csv",
        samples="data/results/models_enet_samples.csv"
    output:
        associations="data/results/associations.csv",
        methnet="data/results/methnet.csv",
        cre="data/results/cluster_score.csv"
    params:
        min_samples=100,
        min_coef=0.0001,
        prom_up=2000,
        prom_dn=200
    conda:
        "envs/postprocess.yaml"
    shell:
        """
        Rscript scripts/process_results.R \
            --min_samples={params.min_samples} \
            --min_coef={params.min_coef} \
            --promoter_up={params.prom_up} \
            --promoter_down={params.prom_dn} \
            {input.net} {input.hyper} {input.coef} {input.samples} data/results
        """

rule get_intergenic_cre:
    input:
        "data/results/cluster_score.csv"
    output:
        "data/annot/intergenic_cre.txt"
    shell:
        "cut -d, -f1 {input} | tail -n +2 > {output}"

rule split_intergenic_cre:
    input:
        "data/annot/intergenic_cre.txt"
    output:
        ["data/clinical/cre{:02d}".format(n) for n in range(CRE_CHUNKS)]
    params:
        chunks=CRE_CHUNKS
    shell:
        "mkdir -p data/clinical && split -d --number=l/{params.chunks} {input} data/clinical/cre"

rule survival_analysis:
    input:
        survival="data/pancan/survival_data.tsv",
        cre="data/clinical/{cre}",
        meth=expand("data/matrices/meth_clusters/{cancer}.rds", cancer=CANCERS)
    output:
        coefs="data/clinical/{cre}_params.csv",
        stats="data/clinical/{cre}_glance.csv"
    conda:
        "envs/postprocess.yaml"
    resources:
        mem_mb=65536,
        cpus_per_task=4,
        runtime=6*60
    script:
        "scripts/survival_analysis.R"

rule combine_survival:
    input:
        ["data/clinical/cre{:02d}_params.csv".format(n) for n in range(CRE_CHUNKS)],
        ["data/clinical/cre{:02d}_glance.csv".format(n) for n in range(CRE_CHUNKS)]
    output:
        coefs="data/results/cre_cox_params.csv",
        stats="data/results/cre_cox_glance.csv"
    shell:
        """
        awk 'NR == 1 || FNR > 1' data/clinical/*_params.csv > {output.coefs} && \
        awk 'NR == 1 || FNR > 1' data/clinical/*_glance.csv > {output.stats}
        """

rule figures_2to5:
    input:
        net=f"data/genomic-networks/pcgene-cluster_{CLUSTER_SCHEME}_network_{RADIUS}.tsv",
        gene_cgi_net="data/genomic-networks/gene_cgi_overlaps.csv",
        manifest="data/annot/meth-cluster_annotated.csv",
        loops="data/annot/meth-cluster_loops.csv",
        tfbs="data/annot/meth-cluster_TFBS.csv",
        cox="data/results/cre_cox_params.csv",
        meth_summary="data/annot/meth-cluster_beta_summary.csv",
        hyper="data/results/models_enet_hyperparams.csv",
        associations="data/results/associations.csv",
        methnet="data/results/methnet.csv",
        hubs="data/results/cluster_score.csv",
        samples="data/results/models_enet_samples.csv",
        examples="data/results/interaction_examples.csv"
    params:
        min_samples=100
    resources:
        runtime=32
    output:
        "data/figures/figures.html"
    conda:
        "envs/postprocess.yaml"
    script:
        "scripts/figures.Rmd"

rule figure_6:
    input:
        genes="data/annot/hugo_hg19_pancan_augmented.tsv",
        clusters="data/genomic-networks/absolute_clusters_probes.tsv",
        QC="data/experiments/HiCUP_summary_report.tsv",
        loops="data/experiments/pcChIP_loops.tsv.gz",
        associations="data/results/methnet.csv"
    params:
        cluster_scheme=CLUSTER_SCHEME,
        min_score=5.0
    output:
        "data/figures/figure_6.html"
    conda:
        "envs/capture-hic.yaml"
    resources:
        mem_mb=8192,
        cpus_per_task=4
    script:
        "scripts/figure_6.Rmd"

rule download_perturb_seq_results:
    output:
        "data/experiments/perturb_seq_results.h5"
    shell:
        "wget 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE236304&format=file&file=GSE236304%5Ffiltered%5Ffeature%5Fbc%5Fmatrix%2Eh5' -O {output}"

rule download_protospacer_umi_thresholds:
    output:
        "data/experiments/protospacer_umi_thresholds.csv"
    shell:
        "wget 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE236304&format=file&file=GSE236304%5Fprotospacer%5Fumi%5Fthresholds%2Ecsv%2Egz' -O {output}"

rule figure_7:
    input:
        net=f"data/genomic-networks/pcgene-cluster_{CLUSTER_SCHEME}_network_{RADIUS}.tsv",
        methnet="data/results/methnet.csv",
        umi_thres="data/experiments/protospacer_umi_thresholds.csv",
        h5file="data/experiments/perturb_seq_results.h5"
    output:
        report="data/figures/figure_7.html"
    params:
        Nboot=20000
    conda:
        "envs/perturb-seq.yaml"
    resources:
        mem_mb=16384,
        cpus_per_task=8,
        runtime=60*24*1
    script:
        "scripts/figure_7.Rmd"
