# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: qc.smk
# time: 2022/01/04
import json

import pysam
import pandas as pd


rule total_reads:
    output:
        temp(config['workspace'] + '/samples/{sample}/qc/{sample}_total_reads.json')
    input:
        bam=rules.merge_bam.output,
        index=rules.index.output
    run:
        total = pysam.view('-c', f'{input.bam}')
        with open(output[0], 'w') as f:
            json.dump({'total_reads': int(total)}, f)


rule mapped_reads:
    output:
        temp(config['workspace'] + '/samples/{sample}/qc/{sample}_mapped_reads.json')
    input:
        bam=rules.merge_bam.output,
        index=rules.index.output
    run:
        mapped = pysam.view('-c', '-F', '4', f'{input.bam}')
        with open(output[0], 'w') as f:
            json.dump({'mapped_reads': int(mapped)}, f)


rule dup_reads:
    output:
        temp(config['workspace'] + '/samples/{sample}/qc/{sample}_dup_reads.json')
    input:
        bam=rules.merge_bam.output,
        index=rules.index.output
    run:
        dup = pysam.view('-c', '-f', '1024', f'{input.bam}')
        with open(output[0], 'w') as f:
            json.dump({'dup_reads': int(dup)}, f)


rule chrM_reads:
    output:
        temp(config['workspace'] + '/samples/{sample}/qc/{sample}_chrM_reads.json')
    input:
        bam=rules.merge_bam.output,
        index=rules.index.output
    run:
        chrm = pysam.view('-c', f'{input.bam}', 'chrM')
        with open(output[0], 'w') as f:
            json.dump({'chrM_reads': int(chrm)}, f)


rule clean_reads:
    output:
        temp(config['workspace'] + '/samples/{sample}/qc/{sample}_clean_reads.json')
    input:
        bam=rules.merge_bam.output,
        index=rules.index.output,
        bed=rules.clean_bed.output
    params:
        include=lambda wildcards: config['filter']['include'],
        exclude=lambda wildcards: config['filter']['exclude'],
        mapq=lambda wildcards: config['filter']['mapq']
    run:
        clean = pysam.view(
            '-c', '-f', f'{params.include}', '-F', f'{params.exclude}', '-q', f'{params.mapq}',
            '-L', f'{input.bed}', f'{input.bam}'
        )
        with open(output[0], 'w') as f:
            json.dump({'clean_reads': int(clean)}, f)


rule promoters:
    output:
        config['workspace'] + '/share/promoters.bed'
    input:
        config['genome']['txdb']
    params:
        script=lambda wildcards: BASE_DIR + '/tools/promoters.R'
    shell:
        'Rscript {params.script} {input} {output}'


rule promoter_reads:
    output:
        temp(config['workspace'] + '/samples/{sample}/{sample}_promoter_reads.json')
    input:
        bam=rules.merge_bam.output,
        index=rules.index.output,
        bed=rules.clean_bed.output,
        promoters=rules.promoters.output
    params:
        include=lambda wildcards: config['filter']['include'],
        exclude=lambda wildcards: config['filter']['exclude'],
        mapq=lambda wildcards: config['filter']['mapq']
    run:
        promoter = pysam.view(
            '-c', '-f', f'{params.include}', '-F', f'{params.exclude}', '-q', f'{params.mapq}',
            '-L', f'{input.promoters}', f'{input.bam}'
        )
        with open(output[0], 'w') as f:
            json.dump({'promoter_reads': int(promoter)}, f)


def load_json(filename):
    with open(filename, 'r') as f:
        return json.load(f)


rule qc:
    output:
        config['workspace'] + '/samples/{sample}/qc/align_summary.txt'
    input:
        total=rules.total_reads.output,
        mapped=rules.mapped_reads.output,
        dup=rules.dup_reads.output,
        chrm=rules.chrM_reads.output,
        clean=rules.clean_reads.output,
        promoter=rules.promoter_reads.output
    run:
        # collect results
        total = load_json(input[0])
        mapped = load_json(input[1])
        dup = load_json(input[2])
        chrm = load_json(input[3])
        clean = load_json(input[4])
        promoter = load_json(input[5])

        qc = total | mapped | dup | chrm | clean | promoter

        qc = pd.Series(qc, name='counts').to_frame()
        qc['frac'] = qc['counts'] / total['total_reads']

        qc.to_csv(output[0], sep='\t', header=False)


rule fragment_sizes:
    output:
        config['workspace'] + '/samples/{sample}/qc/{sample}_fragment_sizes_counts.txt'
    input:
        bam=rules.merge_bam.output,
        index=rules.index.output,
        bed=rules.clean_bed.output
    params:
        script=lambda wildcards: BASE_DIR + '/tools/fragmentsizes.py',
        include=lambda wildcards: config['filter']['include'],
        exclude=lambda wildcards: config['filter']['exclude'],
        mapq=lambda wildcards: config['filter']['mapq']
    shell:
        'python {params.script} -f {params.include} -F {params.exclude} -q {params.mapq}'
        ' -L {input.bed} {input.bam} {output}'
