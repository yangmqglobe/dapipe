# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: qc.smk
# time: 2022/01/04
import json
import re
import os

from jinja2 import Environment, FileSystemLoader
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
        temp(config['workspace'] + '/samples/{sample}/qc/{sample}_promoter_reads.json')
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


rule peak_reads:
    output:
        temp(config['workspace'] + '/samples/{sample}/qc/{sample}_peak_reads.json')
    input:
        cutsites=rules.sort_cutsites.output,
        peaks=rules.narrowPeak2bed.output
    run:
        stdout = shell(
            f'bedtools intersect -wa -u -a {input.cutsites} -b {input.peaks} | wc -l',
            iterable=True
        )
        peak = ''.join(stdout)
        with open(output[0], 'w') as f:
            json.dump({'peak_reads': int(peak)}, f)


rule fragment_sizes:
    output:
        config['workspace'] + '/samples/{sample}/qc/{sample}_fragment_sizes_counts.json'
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


def load_json(filename):
    filename = str(filename)
    with open(filename, 'r') as f:
        return json.load(f)


def grouped_int(value):
    return f'{value:,}'


def percentage(value):
    return f'{value:.2%}'


def get_template():
    env = Environment(loader=FileSystemLoader(f'{BASE_DIR}/tools'))
    env.filters['grouped_int'] = grouped_int
    env.filters['percentage'] = percentage
    return env.get_template('qctemplate.html')


def get_all_library_qc(wildcards):
    return [
        config['workspace'] + f'/samples/{wildcards.sample}/qc/{library}_fastp.json'
        for library in config['samples'][wildcards.sample]['fastq']
    ]


def read_library_qc(filename):
    filename = str(filename)
    library = re.match(r'^(.+?)_fastp.json$', os.path.basename(filename)).group(1)
    with open(filename) as f:
        data = json.load(f)
    return {
        'library' : library,
        'total_reads': data['summary']['before_filtering']['total_reads'],
        'pass_reads': data['summary']['after_filtering']['total_reads'],
        'total_bases': data['summary']['before_filtering']['total_bases'],
        'pass_bases': data['summary']['after_filtering']['total_bases'],
        'total_q20_bases': data['summary']['before_filtering']['q20_bases'],
        'pass_q20_bases': data['summary']['after_filtering']['q20_bases'],
        'total_q30_bases': data['summary']['before_filtering']['q30_bases'],
        'pass_q30_bases': data['summary']['after_filtering']['q30_bases'],
        'pass_gc_content': data['summary']['after_filtering']['gc_content']
    }


rule qc:
    output:
        config['workspace'] + '/samples/{sample}/qc/{sample}_qc.josn',
        config['workspace'] + '/samples/{sample}/qc/{sample}_qc.html'
    input:
        fastp=get_all_library_qc,
        total=rules.total_reads.output,
        mapped=rules.mapped_reads.output,
        dup=rules.dup_reads.output,
        chrm=rules.chrM_reads.output,
        clean=rules.clean_reads.output,
        promoter=rules.promoter_reads.output,
        peak=rules.peak_reads.output,
        fragment_sizes=rules.fragment_sizes.output
    params:
        template_dir=lambda wildcards: BASE_DIR + '/tools'
    run:
        # collect sequencing results
        sequencing = [read_library_qc(library) for library in input.fastp]

        # collect align results
        total = load_json(input.total)
        mapped = load_json(input.mapped)
        dup = load_json(input.dup)
        chrm = load_json(input.chrm)
        clean = load_json(input.clean)
        promoter = load_json(input.promoter)
        peak = load_json(input.peak)
        fragment_sizes = load_json(input.fragment_sizes)

        align = total | mapped | dup | chrm | clean | promoter | peak | fragment_sizes
        
        qc = {
            'sample': wildcards.sample,
            'sequencing': sequencing,
            'align': align
        }

        # write json
        with open(output[0], 'w') as f:
            json.dump(qc, f)
        # render report
        template = get_template()
        report = template.render(qc=qc)

        with open(output[1], 'w') as f:
            f.write(report)
