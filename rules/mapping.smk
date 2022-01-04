# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: mapping.smk
# time: 2021/04/16
import os


rule mapping:
    output:
        temp(config['workspace'] + '/samples/{sample}/mapping/{library}.sorted.bam')
    input:
        fq1=rules.trim.output.fq1,
        fq2=rules.trim.output.fq2
    priority: 10
    log:
       bwa=config['workspace'] + '/log/mapping/{sample}/{library}_bwa.log',
       samblaster=config['workspace'] + '/log/mapping/{sample}/{library}_samblaster.log'
    params:
        rg='\'@RG\\tID:{library}\\tSM:{sample}\\tLB:{library}\\tPL:ILLUMINA\'',
        index=config['genome']['bwa_index'],
        tmp=config['workspace'] + '/samples/{sample}/mapping/{library}.tmp',
    threads: 8 if workflow.cores >= 8 else workflow.cores
    shell:
        'bwa mem -t {threads} -R {params.rg} {params.index} {input.fq1} {input.fq2} 2>{log.bwa}'
        ' | samblaster 2>{log.samblaster}'
        ' | samtools sort -T {params.tmp} -o {output} -'


def get_all_bam(wildcards):
    return [
        config['workspace'] + f'/samples/{wildcards.sample}/mapping/{library}.sorted.bam'
        for library in config['samples'][wildcards.sample]['fastq']
    ]


rule merge_bam:
    output:
        config['workspace'] + '/samples/{sample}/mapping/{sample}.sorted.bam'
    input:
        get_all_bam
    shell:
        'samtools merge {output} {input}'


rule index:
    output:
        config['workspace'] + '/samples/{sample}/mapping/{sample}.sorted.bam.bai'
    input:
        rules.merge_bam.output
    priority: 5
    shell:
        'samtools index {input}'


rule chrom_sizes:
    output:
        config['workspace'] + '/samples/{sample}/mapping/{sample}.chrom_sizes'
    input:
        bam=rules.merge_bam.output,
        index=rules.index.output
    shell:
        'samtools idxstats {input.bam}'
        ' | grep \'^chr[0-9XY]\{{1,2\}}[[:space:]]\''
        ' | awk \'BEGIN{{OFS="\\t"}}{{print $1, $2}}\' > {output}'


rule clean_bed:
    output:
        config['workspace'] + '/samples/{sample}/mapping/{sample}_clean.bed'
    input:
        rules.chrom_sizes.output
    params:
        blacklist=config['genome']['blacklist']
    shell:
        'bedtools complement -g {input} -i {params.blacklist} > {output}'