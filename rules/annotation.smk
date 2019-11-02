rule snpeff:
    input:
        "output/filtered/all.vcf.gz",
    output:
        vcf="output/annotated/all.vcf.gz",
        csvstats="output/snpeff/all.csv",
        stats="output/snpeff/all.html"
    log:
        "logs/snpeff/snpeff.log"
    params:
        reference=config["ref"]["name"],
        extra="-Xmx6g"
    conda:
        "../envs/snpeff.yaml"
    benchmark:
        "benchmarks/annotation/snpeff.json"
    shell:"""
        snpEff {params.extra} -stats {output.stats} -csvStats {output.csvstats} -v {params.reference} {input} > {output.vcf}
    """