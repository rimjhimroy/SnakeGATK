rule vcf_to_tsv:
    input:
        "output/annotated/all.vcf.gz"
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
