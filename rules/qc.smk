rule samtools_stats:
    input:
        "data/bam/{sample}.bam"
    output:
        "output/qc/samtools-stats/{sample}.txt"
    log:
        "logs/samtools-stats/{sample}.log"
    wrapper:
        "0.27.1/bio/samtools/stats"



