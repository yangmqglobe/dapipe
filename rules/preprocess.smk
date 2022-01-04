# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: preprocess.smk
# time: 2021/04/15
import os


def get_adapter(wildcards):
    if 'adapter' in config:
        return config['adapter']
    else:
        return  f'{BASE_DIR}/resources/adapter.fa'


rule trim:
    output:
        fq1=temp(config['workspace'] + '/samples/{sample}/preprocess/{library}_r1_trimed.fq.gz'),
        fq2=temp(config['workspace'] + '/samples/{sample}/preprocess/{library}_r2_trimed.fq.gz'),
        json=config['workspace'] + '/samples/{sample}/qc/{library}_fastp.json',
        html=config['workspace'] + '/samples/{sample}/qc/{library}_fastp.html'
    input:
        fq1=lambda wildcards: config['samples'][wildcards.sample]['fastq'][wildcards.library]['fq1'],
        fq2=lambda wildcards: config['samples'][wildcards.sample]['fastq'][wildcards.library]['fq2']
    params:
        adapter=get_adapter
    log:
        config['workspace'] + '/log/preprocess/{sample}/{library}_fastp.log'
    threads: 8 if workflow.cores > 8 else workflow.cores
    shell:
        'fastp -w {threads} --adapter_fasta {params.adapter}'
        ' -i {input.fq1} -I {input.fq2} -o {output.fq1} -O {output.fq2}'
        ' -j {output.json} -h {output.html} 2>{log}'