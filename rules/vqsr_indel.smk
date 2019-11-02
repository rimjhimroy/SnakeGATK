rule recalibrate_calls_indel:
    input:
        vcf="output/filtered/all.indels.vcf.gz",
        ref=config["ref"]["genome"],
        # resources have to be given as named input files
        hapmap="data/ref/hapmap_3.3.b37.vcf.gz",
        omni="data/ref/1000G_omni2.5.b37.vcf.gz",
        g1k="data/ref/1000G_phase1.snps.high_confidence.b37.vcf.gz",
        dbsnp="data/ref/dbsnp_138.b37.vcf.gz",
        mills="data/ref/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
        # use aux to e.g. download other necessary file
    output:
            vcf="output/filtered/all.indels.recalibrated.vcf.gz",
            tranches="output/filtered/all.INDELs.tranches",
            recal="output/filtered/all.INDEL.recal.vcf.gz",
            plot="output/filtered/plot.INDEL.R"
    conda:
        "../envs/gatk.yaml"
    log:
        "logs/gatk/variantrecalibrator/INDEL.log"
    benchmark:
        "benchmarks/filtering/recalibrate_calls_indel.json"
    params:
        # set mode, must be either SNP, INDEL or BOTH
        # resource parameter definition. Key must match named input files from above.
        extra="",  # optional
        java_opts="", # optional
    shell:"""
    gatk VariantRecalibrator -R {input.ref} -V {input.vcf} \
    --resource dbsnp,known=true,training=false,truth=false,prior=2.0:{input.dbsnp} \
    --resource mills,known=false,training=true,truth=true,prior=12.0:{input.mills} \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -mode INDEL --output {output.recal} --tranches-file {output.tranches} --rscript-file {output.plot}
    call="gatk ApplyVQSR -R {input.ref} -V {input.vcf} -O {output.vcf} --truth-sensitivity-filter-level 99.0 --tranches-file {output.tranches} --recal-file {output.recal} -mode INDEL"
    echo $call
    eval $call
    """

