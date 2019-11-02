import pandas as pd
from snakemake.utils import min_version

min_version("5.7.1")

report: "../report/workflow.rst"

###### Config file and sample sheets #####
configfile: "config/config.yaml"

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
units = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit"], drop=False)



##### Wildcard constraints #####
wildcard_constraints:
    vartype="snvs|indels",
    sample="|".join(samples.index),
    unit="|".join(units["unit"])


##### Helper functions #####

def get_fai():
    return config["ref"]["genome"] + ".fai"


# contigs in reference genome
def get_contigs():
    return pd.read_table(get_fai(), header=None, usecols=[0], squeeze=True, dtype=str)


def get_sample_bams(wildcards):
    """Get all aligned reads of given sample."""
    return expand("output/recal/{sample}.bam",
                  sample=wildcards.sample)


def get_regions_param(regions=config["processing"].get("regions"), default=""):
    if regions:
        params = "--intervals '{}' ".format(regions)
        padding = config["processing"].get("region-padding")
        if padding:
            params += "--interval-padding {}".format(padding)
        return params
    return default


def get_call_variants_params(wildcards, input):
    return (get_regions_param(regions=config["processing"].get("regions"), default="--intervals {}".format(wildcards.contig)) +
            config["params"]["gatk"]["HaplotypeCaller"])

def get_recal_input(bai=False):
    # case 1: no duplicate removal
    f = "data/bam/{sample}.RG.bam"
    if config["processing"]["remove-duplicates"]:
        # case 2: remove duplicates
        f = "output/dedup/{sample}.bam"
    if bai:
        if config["processing"].get("regions"):
            # case 3: need an index because random access is required
            f += ".bai"
            return f
        else:
            # case 4: no index needed
            return []
    else:
        return f