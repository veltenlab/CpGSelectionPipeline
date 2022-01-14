configfile: "config.yaml"

import pandas as pd
import subprocess

samples = config['samples']
sample_map = pd.read_csv(config['sample_map'])

rule all:
    input:
#        expand("data/external/Cabezas/raw/{sample}/fastq/{sample}_1.fastq",sample=samples),
#        expand("data/external/Cabezas/raw/{sample}/mapped/{sample}.bam/{sample}_1_bismark_bt2_pe.bam",sample=samples),
        expand("data/external/Cabezas/raw/{sample}/extract/",sample=samples),

#rule download:
#    input:
#        geo_acc = config['sample_map']
#    output:
#        folder = "data/external/Cabezas/raw/{sample}/fastq/",
#        first = "data/external/Cabezas/raw/{sample}/fastq/{sample}_1.fastq",
#        second = "data/external/Cabezas/raw/{sample}/fastq/{sample}_2.fastq"
#    params:
#        s = "{sample}",
#        max_size = "50G"
#    log:
#        "data/external/Cabezas/raw/{sample}/fastq/{sample}.log"
#    run:
#        print("Start " + params.s)
#        fastqs = sample_map[sample_map['GEO_Accession (exp)']==params.s]['Run']
#        merge_r1 = 'cat '
#        merge_r2 = "cat "
#        for fastq in fastqs:
#            sra_download = "/users/lvelten/mscherer/conda/envs/sra/bin/prefetch " + fastq + " --output-directory " + output.folder + " --max-size " + params.max_size
#            proc = subprocess.run(sra_download,shell=True)
#            dump_command = "/users/lvelten/mscherer/conda/envs/sra/bin/fastq-dump --split-3 " + output.folder + "/" + fastq + "/" + fastq + ".sra"+" -O " + output.folder
#            proc = subprocess.run(dump_command,shell=True)
#            merge_r1 = merge_r1 + output.folder + "/" + fastq + "_1.fastq "
#            merge_r2 = merge_r2 + output.folder + "/" + fastq + "_2.fastq "

#        merge_r1 = merge_r1 + "> " + output.folder + "/" + params.s + "_1.fastq"
#        merge_r2 = merge_r2 + "> " + output.folder + "/" + params.s + "_2.fastq"
#        proc = subprocess.run(merge_r1,shell=True)
#        proc = subprocess.run(merge_r2,shell=True)
#        cleanup = "rm -rf " + output.folder + "SRR*"
#        proc = subprocess.run(cleanup,shell=True)

#rule bismark:
#    input:
#        ref = "references/mm10/",
#        first = "data/external/Cabezas/raw/{sample}/fastq/{sample}_1.fastq",
#        second = "data/external/Cabezas/raw/{sample}/fastq/{sample}_2.fastq"
#    output:
#        "data/external/Cabezas/raw/{sample}/mapped/{sample}.bam/{sample}_1_bismark_bt2_pe.bam"
#    conda:
#        "/users/lvelten/mscherer/conda/envs/bismark.yml"
#    threads:
#        4
#    log:
#        "data/external/Cabezas/raw/{sample}/mapped/{sample}.log"
#    params:
#        threads = 4,
#        tmp_dir = "tmp"
#    shell:
#        "bismark --bowtie2 {input.ref} -1 {input.first} -2 {input.second} -o {output} --temp_dir {params.tmp_dir} -p {params.threads} ; "
#        "rm -rf {input.first}; rm -rf {input.second} ; "
#        "2> {log}"

rule bismark_extract:
    input:
        "data/external/Cabezas/raw/{sample}/mapped/{sample}.bam/{sample}_1_bismark_bt2_pe.bam/{sample}_1_bismark_bt2_pe.bam"
    output:
        "data/external/Cabezas/raw/{sample}/extract/"
    conda:
        "/users/lvelten/mscherer/conda/envs/bismark.yml"
    log:
        "data/external/Cabezas/raw/{sample}/trimmed/{sample}.log"
    params:
        threads = 10,
        genome = "references/mm10/"
    shell:
        "bismark_methylation_extractor {input} -o {output} --parallel {params.threads} --comprehensive --bedGraph --cytosine_report --no_overlap --genome_folder {params.genome} ; "
        "2> {log}"

rule sort_index:
    input:
        "data/external/Cabezas/raw/{sample}/mapped/{sample}.bam/{sample}_1_bismark_bt2_pe.bam"
    output:
        "data/external/Cabezas/raw/{sample}/mapped/{sample}.bam/{sample}_1_bismark_bt2_pe.bam.bai"
    conda:
        "/users/lvelten/mscherer/conda/envs/bismark.yml"
    shell:
        "samtools index {input} > data/external/Cabezas/raw/{sample}/mapped/{sample}.bam/{sample}_1_bismark_bt2_pe_sorted.bam ; "
        "samtools sort data/external/Cabezas/raw/{sample}/mapped/{sample}.bam/{sample}_1_bismark_bt2_pe_sorted.bam ; "
        "rm -rf data/external/Cabezas/raw/{sample}/mapped/{sample}.bam/{sample}_1_bismark_bt2_pe.bam"
