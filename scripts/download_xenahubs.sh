#!/bin/bash
set -e
# module add python/cpu/3.8.11 r/4.1.1

CANCERS=( "acc" "blca" "brca" "cesc" "chol" "coad" "dlbc" "esca" "gbm" \
          "hnsc" "kich" "kirc" "kirp" "laml" "lgg" "lihc" "luad" "lusc" \
          "meso" "ov" "paad" "pancan" "pcpg" "prad" "read" "sarc" "skcm" \
          "stad" "tgct" "thca" "thym" "ucec" "ucs" "uvm" )

DATA="$1"
OUTDIR="$2"
mkdir -p "$OUTDIR"

if [[ "$DATA" == meth ]]; then
    ./scripts/clone_xena.py gdc -o data/xenahubs
    gdc="$(realpath data/xenahubs/gdc/genomicMatrix)"
    for cancer in "${CANCERS[@]}" ; do
        CANCER="${cancer^^}"
        mat="$gdc"/TCGA-"$CANCER"/Xena_Matrices/TCGA-"$CANCER".methylation450.tsv
        [[ -f "$mat" ]] && ln -s "$mat" "$OUTDIR"/"$cancer".tsv
    done
elif [[ "$DATA" == rna ]]; then
    pancan="$(realpath data/xenahubs/pancanatlas)"
    rna="$pancan"/genomicMatrix/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.tsv
    if [[ ! -f "$rna" ]]; then
        echo "Downloading Xenahubs PanCancer Atlas Cohort"
        ./scripts/clone_xena.py pancanatlas -o data/xenahubs
    fi
    ln -s "$rna" "$OUTDIR"/pancan.tsv
elif [[ "$DATA" == clinical ]]; then
    pancan="$(realpath data/xenahubs/pancanatlas)"
    clinical="$pancan"/clinicalMatrix/Survival_SupplementalTable_S1_20171025_xena_sp.tsv
    if [[ ! -f "$clinical" ]]; then
        echo "Downloading Xenahubs PanCancer Atlas Cohort"
        ./scripts/clone_xena.py pancanatlas -o data/xenahubs
    fi
    output="$OUTDIR"/survival_data.tsv
    ln -s "$clinical" "$output"
else
    >&2 echo Uknown data type: "$DATA"
    exit 1
fi
