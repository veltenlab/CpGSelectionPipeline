configfile: "config_merge.yaml"

import pandas as pd
import subprocess

samples = config['samples']
sample_map = pd.read_csv(config['sample_map'])

rule all:
    input:
#        expand("../data/external/Cabezas/raw/{sample}/merged/merged.bam",sample=samples),
        expand("../data/external/Cabezas/raw/{sample}/extract/",sample=samples)
#        expand("../data/external/Cabezas/raw/{sample}/merged/merged.bam.bai",sample=samples),

#rule merge:
#    input:
#        geo_acc = config['sample_map'],
#        folder = "../data/external/Cabezas/raw/"        
#    output:
#        merged = "../data/external/Cabezas/raw/{sample}/merged/merged.bam"
#    params:
#        s = "{sample}",
#        max_size = "50G"
#    log:
#        "../data/external/Cabezas/raw/{sample}/fastq/{sample}.log"
#    run:
#        srrs = sample_map[sample_map['GEO_Accession (exp)']==params.s]['Run']
#        all_bams = " "
#        for srr in srrs:
#            all_bams = all_bams + input.folder + "/" + srr + "/mapped/" + srr + ".bam/" + srr + "_1_bismark_bt2_pe.bam/" + srr + "_1_bismark_bt2_pe.bam "

#        merge_cmd = "samtools merge " + output.merged + all_bams + " -n"
#        proc = subprocess.run(merge_cmd,shell=True)
#        #cleanup = "rm -rf " + input.folder + "SRR*"
#        #proc = subprocess.run(cleanup,shell=True)

#rule sort_index:
#    input:
#        "../data/external/Cabezas/raw/{sample}/merged/merged.bam"
#    output:
#        "../data/external/Cabezas/raw/{sample}/merged/merged.bam.bai"
#    conda:
#        "/users/lvelten/mscherer/conda/envs/bismark.yml"
#    params:
#        s = "{sample}"
#    shell:
#        "samtools sort {input} -o ../data/external/Cabezas/raw/{params.s}/merged/merged_sorted.bam ; "
#        "samtools index ../data/external/Cabezas/raw/{params.s}/merged/merged_sorted.bam"

rule bismark_extract:
    input:
        "../data/external/Cabezas/raw/{sample}/merged/merged.bam"
    output:
        "../data/external/Cabezas/raw/{sample}/extract/"
    conda:
        "/users/lvelten/mscherer/conda/envs/bismark.yml"
    log:
        "../data/external/Cabezas/raw/{sample}/trimmed/{sample}.log"
    params:
        threads = 10,
        genome = "../references/mm10/",
        s = "{sample}"
    shell:
        "bismark_methylation_extractor {input} -o {output} --parallel {params.threads} --comprehensive --bedGraph --cytosine_report --no_overlap --genome_folder {params.genome} ;  "
        "rm -rf ../data/external/Cabezas/raw/{params.s}/extract/*_context_merged.txt ; "
        "2> {log}"

