
rule vcf_to_tsv:
    input:
        "output/annotated/all.vcf"
    output:
        report("output/tables/calls.tsv.gz", caption="../report/calls.rst", category="Calls")
    conda:
        "../envs/rbt.yaml"
    benchmark:
        "benchmarks/stats/vcf_to_tsv.json"
    shell:
        "bcftools view --apply-filters PASS --output-type u {input} | "
        "rbt vcf-to-txt -g --fmt DP AD --info ANN | "
        "gzip > {output}"


rule plot_stats:
    input:
        "output/tables/calls.tsv.gz"
    output:
        depths=report("output/plots/depths.svg", caption="../report/depths.rst", category="Plots"),
        freqs=report("output/plots/allele-freqs.svg", caption="../report/freqs.rst", category="Plots")
    conda:
        "../envs/stats.yaml"
    benchmark:
        "benchmarks/stats/plot_stats.json"
    script:
        "../scripts/plot-depths.py"




rule concordance_all:
    input:
        ref=config["ref"]["genome"],
        gvcfs="output/annotated/all.vcf"
    output:
        outsnp="output/concordance/all.concordance.snp.tsv",
        outindel="output/concordance/all.concordance.indel.tsv"
    conda:
        "../envs/gatk.yaml"
    benchmark:
        "benchmarks/concordance/all.concordance.json"
    shell:"""
    gatk Concordance -R {input.ref} -eval {input.gvcfs}--truth data/ref/1000G_phase1.snps.high_confidence.b37.vcf.gz --summary {output.outsnp}
    gatk Concordance -R {input.ref} -eval {input.gvcfs}--truth data/ref/Mills_and_1000G_gold_standard.indels.b37.vcf.gz  --summary {output.outindel}
    """
