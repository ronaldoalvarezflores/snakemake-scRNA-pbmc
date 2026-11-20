configfile: "config/config.yaml"

SAMPLES = list(config["samples"])


rule all:
    input:
        expand("results/{sample}/umap.png", sample=SAMPLES),
        expand("results/{sample}/markers.png", sample=SAMPLES),
        expand("results/{sample}/summary.txt", sample=SAMPLES)


# Shell rule: download the 10x count matrix of one sample and unpack it
rule get_data:
    output:
        directory("data/{sample}")
    params:
        url=lambda wc: config["samples"][wc.sample]
    log:
        "logs/get_data/{sample}.log"
    shell:
        """
        mkdir -p {output}
        wget -nv -O {output}/raw.tar.gz {params.url} > {log} 2>&1
        tar -xzf {output}/raw.tar.gz -C {output} --strip-components=2
        rm {output}/raw.tar.gz
        """


# Shell rule: quick summary of the raw data with standard Unix tools
rule summary:
    input:
        "data/{sample}"
    output:
        "results/{sample}/summary.txt"
    shell:
        """
        mkdir -p $(dirname {output})
        echo "sample: {wildcards.sample}" > {output}
        echo "cells (barcodes): $(wc -l < {input}/barcodes.tsv)" >> {output}
        echo "genes: $(wc -l < {input}/genes.tsv)" >> {output}
        echo "md5 matrix.mtx: $(md5sum {input}/matrix.mtx | cut -d' ' -f1)" >> {output}
        """


# Python script: import, filter and normalize
rule preprocess:
    input:
        "data/{sample}"
    output:
        "results/{sample}/preprocessed.h5ad"
    params:
        min_genes=config["min_genes"],
        min_cells=config["min_cells"],
        max_pct_mt=config["max_pct_mt"]
    resources:
        mem_mb=lambda wc, attempt: attempt * config["mem_mb"]
    retries: 3
    log:
        "logs/preprocess/{sample}.log"
    benchmark:
        "benchmarks/preprocess/{sample}.tsv"
    conda:
        "envs/scanpy.yaml"
    script:
        "scripts/preprocess.py"


# Python script: PCA, neighbours, Leiden clustering, UMAP, annotation
rule cluster:
    input:
        "results/{sample}/preprocessed.h5ad"
    output:
        "results/{sample}/clustered.h5ad"
    params:
        n_top_genes=config["n_top_genes"],
        n_pcs=config["n_pcs"],
        n_neighbors=config["n_neighbors"],
        resolution=config["resolution"],
        markers=config["markers"]
    threads: config["threads"]
    resources:
        mem_mb=lambda wc, attempt: attempt * config["mem_mb"]
    retries: 3
    log:
        "logs/cluster/{sample}.log"
    benchmark:
        "benchmarks/cluster/{sample}.tsv"
    conda:
        "envs/scanpy.yaml"
    script:
        "scripts/cluster.py"


# Python script: final figures
rule plot:
    input:
        "results/{sample}/clustered.h5ad"
    output:
        umap="results/{sample}/umap.png",
        markers="results/{sample}/markers.png"
    params:
        markers=config["markers"]
    log:
        "logs/plot/{sample}.log"
    conda:
        "envs/scanpy.yaml"
    script:
        "scripts/plot.py"


# Step 1 of the exercise
rule test:
    output:
        "test.txt"
    shell:
        "echo 'Hello from Snakemake!' > {output}"
