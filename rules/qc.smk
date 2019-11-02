rule samtools_stats:
    input:
        "data/bam/{sample}.bam"
    output:
        "output/qc/samtools-stats/{sample}.txt"
    log:
        "logs/samtools-stats/{sample}.log"
    benchmark:
        "benchmarks/qc/samtools_stats.{sample}.json"
    wrapper:
        "0.27.1/bio/samtools/stats"



