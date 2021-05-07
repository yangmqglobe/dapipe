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


rule all:
    input:
        expand(
            config['workspace'] + '/samples/{sample}/callpeak/{sample}_treat_pileup.bw',
            sample=config['samples']
        ),
        config['workspace'] + '/aggregate/all_sample_rpkm_qnorm.txt',
        expand(
            config['workspace'] + '/aggregate/all_sample_pcaplot.{fmt}',
            fmt=config['plot_formats']
        ),
        expand(
            config['workspace'] + '/comparisons/{comparison}/{comparison}_result.txt',
            comparison=config['comparisons']
        )