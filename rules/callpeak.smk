# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: callpeak.smk
# time: 2021/04/16
from psutil import virtual_memory


rule cutsites:
    output:
        temp(config['workspace'] + '/samples/{sample}/callpeak/{sample}_cutsites.bed')
    input:
        bam=rules.merge_bam.output,
        index=rules.index.output,
        bed=rules.clean_bed.output
    params:
        script=os.path.dirname(workflow.snakefile) + '/tools/cutsites.py'
    shell:
        'python {params.script} -L {input.bed} {input.bam} {output}'


rule sort_cutsites:
    output:
        config['workspace'] + '/samples/{sample}/callpeak/{sample}_cutsites.sorted.bed'
    input:
        rules.cutsites.output
    shell:
        'sort -k1,1V -k2,2n -k3,3n {input} > {output}'


rule callpeak:
    output:
        narrowPeak=config['workspace'] + '/samples/{sample}/callpeak/{sample}_peaks.narrowPeak',
        treat_bdg=temp(config['workspace'] + '/samples/{sample}/callpeak/{sample}_treat_pileup.bdg')
    input:
        rules.sort_cutsites.output
    log:
        config['workspace'] + '/log/callpeak/{sample}_macs2.log'
    params:
        genome_size=config['genome']['size'],
        outdir=config['workspace'] + '/samples/{sample}/callpeak',
        name='{sample}'
    shell:
        'macs2 callpeak -t {input} -f BED -g {params.genome_size} --keep-dup all --outdir {params.outdir}'
        ' -n {params.name} -B --SPMR --nomodel --shift -75 --extsize 150 --nolambda --call-summits --max-gap 250 2>{log}'


rule narrowPeak2bed:
    output:
        narrowPeak=config['workspace'] + '/samples/{sample}/callpeak/{sample}_peaks.bed'
    input:
        rules.callpeak.output.narrowPeak
    shell:
        'awk \'BEGIN {{OFS="\\t"}}{{print $1,$2,$3,$4,$5}}\' {input} > {output}'
        ' && bedSort {output} {output}'


rule bdg_clip:
    output:
        temp(config['workspace'] + '/samples/{sample}/callpeak/{sample}_treat_pileup_clip.bdg')
    input:
        rules.callpeak.output.treat_bdg,
        rules.chrom_sizes.output
    shell:
        'bedClip {input} {output}'


rule bdg_sort:
    output:
        temp(config['workspace'] + '/samples/{sample}/callpeak/{sample}_treat_pileup_sort.bdg')
    input:
        rules.bdg_clip.output
    shell:
        'bedSort {input} {output}'


rule bdg2bw:
    output:
        config['workspace'] + '/samples/{sample}/callpeak/{sample}_treat_pileup.bw'
    input:
        rules.bdg_sort.output,
        rules.chrom_sizes.output
    shell:
        'bedGraphToBigWig {input} {output}'