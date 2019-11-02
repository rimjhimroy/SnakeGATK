import os

configfile: "config.yaml"

extensions=['.fa.amb','.fa.fai','.dict']
GENOMEBASE=os.path.splitext(config["ref"]["genome"])[0]
#DBSNP=config["ref"]["known-variants"]
#DBSNPBASE=os.path.splitext(config["ref"]["known-variants"])[0]

rule all:
    input: [GENOMEBASE+ext for ext in extensions]

rule initialize_index:
    input:
        config["ref"]["genome"]
    output:
        GENOMEBASE+".fa.amb"
    wrapper:
        "0.27.1/bio/bwa/index"

rule initialize_fa_index:
    input:
        config["ref"]["genome"]
    output:
        GENOMEBASE+".fa.fai"
    conda:
        "../envs/samtools.yaml"
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
         conda:
            "../envs/gatk.yaml"
         shell:
            "gatk CreateSequenceDictionary -R {input}"

exists = os.path.isfile(DBSNP)
if exists:
    index_exist=os.path.isfile(DBSNP+".tbi")
    if index_exist:
        print("Dictionary file already exists! Skipping this step. If dictionary file is intended to be re-generated please delete current dictionary file")
    else:
        rule generate_tabix:
            input:
                config["ref"]["known-variants"]
            output:
                DBSNP+".tbi"
            conda:
                "../envs/gatk.yaml"
            params:
                base=DBSNPBASE
            shell:"""
                gunzip {params.base}.gz
                bgzip {params.base}
                tabix -p vcf {params.base}.gz
            """

