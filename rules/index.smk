import os

configfile: "config/config.yaml"

args = {}
extensions=['.fa.amb','.fa.fai','.dict']
GENOMEBASE=os.path.splitext(config["ref"]["genome"])[0]
args['VCFS'] = glob_wildcards("data/ref/{S}.vcf.gz").S
print(args['VCFS'])

rule all:
    input: 
        [GENOMEBASE+ext for ext in extensions],
        expand('data/ref/{vcfs}.vcf.gz.tbi', zip, vcfs=args['VCFS'])

rule generate_tabix:
    input:
        expand('data/ref/{vcf}.vcf.gz', zip, vcfs=args['VCFS'])
    output:
        "data/ref/{vcf}.vcf.gz.tbi"
    benchmark:
        "benchmarks/index/generate_tabix.{vcf}.json"
    conda:
        "../envs/tabix.yaml"
    shell:"""
        gunzip < data/ref/{wildcards.vcf}.vcf.gz > data/ref/{wildcards.vcf}.temp.vcf
        bgzip data/ref/{wildcards.vcf}.temp.vcf
        mv data/ref/{wildcards.vcf}.temp.vcf.gz data/ref/{wildcards.vcf}.vcf.gz
        tabix -p vcf data/ref/{wildcards.vcf}.vcf.gz
    """

rule initialize_index:
    input:
        config["ref"]["genome"]
    output:
        GENOMEBASE+".fa.amb"
    benchmark:
        "benchmarks/index/initialize_index.json"
    wrapper:
        "0.27.1/bio/bwa/index"

rule initialize_fa_index:
    input:
        config["ref"]["genome"]
    output:
        GENOMEBASE+".fa.fai"
    conda:
        "../envs/samtools.yaml"
    benchmark:
        "benchmarks/index/initialize_fa_index.json"
    shell:
        "samtools faidx {input}"

exists = os.path.isfile(GENOMEBASE+".dict")
if exists:
    print("Dictionary file already exists! Skipping this step. If dictionary file is intended to be re-generated please delete current dictionary file")
else:
    rule generate_dictionary:
        input:
            config["ref"]["genome"]
        output:
            GENOMEBASE+".dict"
        benchmark:
            "benchmarks/index/generate_dictionary.json"
        conda:
            "../envs/gatk.yaml"
        shell:
            "gatk CreateSequenceDictionary -R {input}"




