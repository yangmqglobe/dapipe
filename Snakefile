# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: Snakefile
# time: 2021/04/15
include: 'rules/preprocess.smk'
include: 'rules/mapping.smk'
include: 'rules/callpeak.smk'
include: 'rules/aggregate.smk'
include: 'rules/comparison.smk'
include: 'rules/qc.smk'


BASE_DIR = os.path.dirname(workflow.snakefile)


rule basic_all:
    input:
        expand(
            config['workspace'] + '/samples/{sample}/qc/align_summary.txt',
            sample=config['samples']
        ),
        expand(
            config['workspace'] + '/samples/{sample}/qc/{sample}_fragment_sizes_counts.txt',
            sample=config['samples']
        ),
        expand(
            config['workspace'] + '/samples/{sample}/callpeak/{sample}_peaks.bed',
            sample=config['samples']
        ),
        expand(
            config['workspace'] + '/samples/{sample}/callpeak/{sample}_treat_pileup.bw',
            sample=config['samples']
        )


rule aggregate_all:
    input:
        config['workspace'] + '/aggregate/all_sample_rpkm_qnorm.txt',
        expand(
            config['workspace'] + '/aggregate/all_sample_pcaplot.{fmt}',
            fmt=config['plot_formats']
        ),
        config['workspace'] + '/aggregate/all_uniq_peaks_annotated.bed'


rule compare_all:
    input:
        expand(
            config['workspace'] + '/comparisons/{comparison}/{comparison}_result.txt',
            comparison=config['comparisons'] if 'comparisons' in config else list()
        )


rule all:
    input:
        rules.basic_all.input,
        rules.aggregate_all.input,
        rules.comparisons.input
        