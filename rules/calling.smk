print(config["ref"]["known-variants"])
rule call_variants:
    input:
        bam=get_sample_bams,
        ref=config["ref"]["genome"],
        known=config["ref"]["known-variants"]
    output:
        gvcf=protected("output/called/{sample}.{contig}.g.vcf.gz")
    log:
        "logs/gatk/haplotypecaller/{sample}.{contig}.log"
    params:
        extra=get_call_variants_params
    benchmark:
        "benchmarks/calling/call_variants.{sample}.{contig}.json"
    wrapper:
        "0.27.1/bio/gatk/haplotypecaller"


rule combine_calls:
    input:
        conc=expand("output/concordance/{sample}.{{contig}}.concordance.snp.tsv",sample=samples.index),
        ref=config["ref"]["genome"],
        gvcfs=expand("output/called/{sample}.{{contig}}.g.vcf.gz", sample=samples.index)
    output:
        gvcf="output/called/all.{contig}.g.vcf.gz"
    log:
        "logs/gatk/combinegvcfs.{contig}.log"
    benchmark:
        "benchmarks/calling/combine_calls.{contig}.json"
    wrapper:
        "0.27.1/bio/gatk/combinegvcfs"

rule concordance:
    input:
        ref=config["ref"]["genome"],
        gvcfs="output/called/{sample}.{contig}.g.vcf.gz"
    output:
        outsnp="output/concordance/{sample}.{contig}.concordance.snp.tsv",
        outindel="output/concordance/{sample}.{contig}.concordance.indel.tsv"
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "benchmarks/concordance/{sample}.{contig}.concordance.json"
    shell:"""
    gatk Concordance -R {input.ref} -eval {input.gvcfs} --truth data/ref/1000G_phase1.snps.high_confidence.b37.vcf.gz --summary {output.outsnp}
    gatk Concordance -R {input.ref} -eval {input.gvcfs} --truth data/ref/Mills_and_1000G_gold_standard.indels.b37.vcf.gz  --summary {output.outindel}
    """


rule genotype_variants:
    input:
        ref=config["ref"]["genome"],
        gvcf="output/called/all.{contig}.g.vcf.gz"
    output:
        vcf=temp("output/genotyped/all.{contig}.vcf.gz")
    params:
        extra=config["params"]["gatk"]["GenotypeGVCFs"]
    log:
        "logs/gatk/genotypegvcfs.{contig}.log"
    benchmark:
        "benchmarks/calling/genotype_variants.{contig}.json"
    wrapper:
        "0.27.1/bio/gatk/genotypegvcfs"


rule merge_variants:
    input:
        ref=get_fai(), # fai is needed to calculate aggregation over contigs below
        vcfs=lambda w: expand("output/genotyped/all.{contig}.vcf.gz", contig=get_contigs())
    output:
        vcf="output/genotyped/all.vcf.gz"
    log:
        "logs/picard/merge-genotyped.log"
    benchmark:
        "benchmarks/calling/merge_variants.json"
    wrapper:
        "0.40.2/bio/picard/mergevcfs"
