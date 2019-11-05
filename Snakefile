##### Target rules #####

rule all:
    input:
        "output/annotated/all.vcf",
        "output/plots/depths.svg",
        "output/plots/allele-freqs.svg",
        "output/concordance/all.concordance.snp.tsv"

##### Modules #####

include: "rules/common.smk"
include: "rules/processbam.smk"
include: "rules/calling.smk"
include: "rules/filtering.smk"
include: "rules/stats.smk"
include: "rules/qc.smk"
include: "rules/annotation.smk"
