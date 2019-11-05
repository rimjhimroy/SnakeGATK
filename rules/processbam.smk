rule add_readgroup:
    input:
        "data/bam/{sample}.bam"
    output:
        "data/bam/{sample}.RG.bam"
    log:
        "logs/picard/readgroup/{sample}.log"
    conda:
        "../envs/picard.yaml"
    benchmark:
        "benchmarks/processbam/add_readgroup.{sample}.json"
    shell:
        "picard AddOrReplaceReadGroups INPUT={input} OUTPUT={output} SORT_ORDER=coordinate RGID={wildcards.sample}.1 RGLB={wildcards.sample} RGPL=illumina RGPU=unit1 RGSM={wildcards.sample} VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true"

rule mark_duplicates:
    input:
        "data/bam/{sample}.RG.bam"
    output:
        bam="output/dedup/{sample}.bam",
        metrics="output/qc/dedup/{sample}.metrics.txt"
    log:
        "logs/picard/dedup/{sample}.log"
    params:
        config["params"]["picard"]["MarkDuplicates"]
    conda:
        "../envs/picard.yaml"
    benchmark:
        "benchmarks/processbam/mark_duplicates.{sample}.json"
    shell:
        "picard MarkDuplicates INPUT={input} OUTPUT={output.bam} METRICS_FILE={output.metrics} VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true"


rule recalibrate_base_qualities:
    input:
        bam=get_recal_input(),
        bai=get_recal_input(bai=True),
        ref=config["ref"]["genome"],
        known=config["ref"]["known-variants"]
    output:
        bam="output/recal/{sample}.bam"
    params:
        extra=get_regions_param() + config["params"]["gatk"]["BaseRecalibrator"]
    log:
        "logs/gatk/bqsr/{sample}.log"
    benchmark:
        "benchmarks/processbam/recalibrate_base_qualities.{sample}.json"
    wrapper:
        "0.27.1/bio/gatk/baserecalibrator"


rule samtools_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    benchmark:
        "benchmarks/processbam/{prefix}.samtools_index.json"
    wrapper:
        "0.27.1/bio/samtools/index"