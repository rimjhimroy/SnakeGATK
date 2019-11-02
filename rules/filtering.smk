import os
import subprocess

def get_vartype_arg(wildcards):
    return "--select-type-to-include {}".format(
        "SNP" if wildcards.vartype == "snvs" else "INDEL")


rule select_calls:
    input:
        ref=config["ref"]["genome"],
        vcf="output/genotyped/all.vcf.gz"
    output:
        vcf=temp("output/filtered/all.{vartype}.vcf.gz")
    params:
        extra=get_vartype_arg
    log:
        "logs/gatk/selectvariants/{vartype}.log"
    benchmark:
        "benchmarks/filtering/select_calls.{vartype}.json"
    wrapper:
        "0.27.1/bio/gatk/selectvariants"


def get_filter(wildcards):
    return {
        "snv-hard-filter":
        config["filtering"]["hard"][wildcards.vartype]}


rule hard_filter_calls:
    input:
        ref=config["ref"]["genome"],
        vcf="output/filtered/all.{vartype}.vcf.gz"
    output:
        vcf="output/filtered/all.{vartype}.hardfiltered.vcf.gz"
    params:
        filters=get_filter
    log:
        "logs/gatk/variantfiltration/{vartype}.log"
    benchmark:
        "benchmarks/filtering/hard_filter_calls.{vartype}.json"
    wrapper:
        "0.27.1/bio/gatk/variantfiltration"

rule recalibrate_calls_snps:
    input:
        vcf="output/filtered/all.snvs.vcf.gz",
        ref=config["ref"]["genome"],
        # resources have to be given as named input files
        hapmap="data/ref/hapmap_3.3.b37.vcf.gz",
        omni="data/ref/1000G_omni2.5.b37.vcf.gz",
        g1k="data/ref/1000G_phase1.snps.high_confidence.b37.vcf.gz",
        dbsnp="data/ref/dbsnp_138.b37.vcf.gz",
        mills="data/ref/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
        # use aux to e.g. download other necessary file
    output:
            vcf="output/filtered/all.snvs.recalibrated.vcf.gz",
            tranches="output/filtered/all.snvs.tranches",
            recal="output/filtered/all.snvs.recal.vcf.gz",
            plot="output/filtered/plot.snvs.R"
    conda:
        "../envs/gatk.yaml"
    log:
        "logs/gatk/variantrecalibrator/snvs.log"
    benchmark:
        "benchmarks/filtering/recalibrate_calls_snps.json"
    params:
        # set mode, must be either SNP, INDEL or BOTH
        # resource parameter definition. Key must match named input files from above.
        extra="",  # optional
        java_opts="", # optional
    shell:"""
    gatk VariantRecalibrator -R {input.ref} -V {input.vcf} \
    --resource dbsnp,known=true,training=false,truth=false,prior=2.0:{input.dbsnp} \
    --resource 1000G,known=false,training=true,truth=false,prior=10.0:{input.g1k} \
    --resource hapmap,known=false,training=true,truth=true,prior=15.0:{input.hapmap} \
    --resource omni,known=false,training=true,truth=false,prior=12.0:{input.omni} \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode SNP --output {output.recal} --tranches-file {output.tranches} --rscript-file {output.plot}
    gatk ApplyVQSR -R {input.ref} -V {input.vcf} -O {output.vcf} --truth-sensitivity-filter-level 99.0 --tranches-file {output.tranches} --recal-file {output.recal} -mode SNP
    """

#include: rules/vqsr_indel.smk


rule merge_calls:
    input:
        vcfs=["output/filtered/all.indels.hardfiltered.vcf.gz","output/filtered/all.snvs.recalibrated.vcf.gz"]
    output:
        vcf="output/filtered/all.vcf.gz"
    log:
        "logs/picard/merge-filtered.log"
    benchmark:
        "benchmarks/filtering/merge_calls.json"
    conda:
        "../envs/picard.yaml"
    wrapper:
        "0.27.1/bio/picard/mergevcfs"
