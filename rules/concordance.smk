rule concordance:
    input:
        ref=config["ref"]["genome"],
        gvcfs="output/called/{sample}.{contig}.g.vcf.gz"
    output:
        outsnp="output/concordance/{{sample}}.{{contig}}.concordance.snp.tsv",
        outindel="output/concordance/{{sample}}.{{contig}}.concordance.indel.tsv"
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "benchmarks/concordance/{sample}.{contig}.concordance.json"
    script:"""
    gatk Concordance -R {input.ref} -eval {input.gvcfs}--truth data/ref/1000G_phase1.snps.high_confidence.b37.vcf.gz --summary {output.outsnp}
    gatk Concordance -R {input.ref} -eval {input.gvcfs}--truth data/ref/Mills_and_1000G_gold_standard.indels.b37.vcf.gz  --summary {output.outindel}
    """