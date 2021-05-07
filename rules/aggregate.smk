# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: aggregate.smk
# time: 2021/04/16
rule merge_peaks:
    output:
        config['workspace'] + '/aggregate/all_uniq_peaks.bed'
    input:
        expand(
            config['workspace'] + '/samples/{sample}/callpeak/{sample}_peaks.narrowPeak',
            sample=config['samples']
        )
    shell:
        'cat {input} | sort -k1,1V -k2,2n -k3,3n'
        ' | bedtools merge -i stdin'
        ' | awk \'BEGIN{{OFS="\\t"}}{{print $1, $2, $3, "PEAK_"NR, ".", "."}}\''
        ' > {output}'


rule annotate:
    output:
        config['workspace'] + '/aggregate/all_uniq_peaks_annotated.bed'
    input:
        rules.merge_peaks.output,
        config['genome']['txdb']
    params:
        script=os.path.dirname(workflow.snakefile) + '/tools/annotate.R'
    shell:
        'Rscript {params.script} {input} {output}'


rule count:
    output:
        temp(config['workspace'] + '/aggregate/{sample}_cutsites_counts.bed')
    input:
        peaks=rules.merge_peaks.output,
        cutsites=rules.sort_cutsites.output
    params:
        names='{sample}'
    shell:
        'bedtools intersect -sorted -wa -C -a {input.peaks} -b {input.cutsites} > {output}'


rule merge_counts:
    output:
        raw=config['workspace'] + '/aggregate/all_sample_raw_counts.txt',
        rpkm=config['workspace'] + '/aggregate/all_sample_rpkm.txt',
        qnorm=config['workspace'] + '/aggregate/all_sample_rpkm_qnorm.txt'
    input:
        counts=expand(
            config['workspace'] + '/aggregate/{sample}_cutsites_counts.bed',
            sample=config['samples']
        ),
        peaks=rules.merge_peaks.output,
        cutsites=expand(
            config['workspace'] + '/samples/{sample}/callpeak/{sample}_cutsites.sorted.bed',
            sample=config['samples']
        )
    params:
        names=expand('{sample}', sample=config['samples'])
    run:
        from snakemake.utils import linecount
        import pandas as pd
        import qnorm

        counts = [
            pd.read_table(file, usecols=[3, 6], names=['peak', name], index_col=0, squeeze=True)
            for name, file in zip(params.names, input.counts)
        ]
        counts = pd.concat(counts, axis=1)
        counts.to_csv(output.raw, sep='\t', index_label='')

        peaks = pd.read_table(
            input.peaks[0], usecols=[1, 2, 3], names=['start', 'end', 'peak'], index_col=2
        )
        length = peaks['end'] - peaks['start']
        rpk = counts.divide(length / 1000, axis=0)

        total = pd.Series({
            name: linecount(file)
            for name, file in zip(params.names, input.cutsites)
        })
        rpkm = rpk.divide(total / 1000000, axis=1)
        rpkm.to_csv(output.rpkm, sep='\t', index_label='')

        norm = qnorm.quantile_normalize(rpkm)
        norm.to_csv(output.qnorm, sep='\t', index_label='')


rule metadata:
    output:
        config['workspace'] + '/aggregate/metadata.txt',
    input:
        rules.merge_counts.output
    run:
        import pandas as pd

        meta = {sample: meta['meta'] for sample, meta in config['samples'].items()}
        meta = pd.DataFrame.from_dict(meta, orient='index')

        meta_cols = [
            col for col in meta.columns if meta[col].unique().shape[0] != 1
        ]
        if len(meta_cols) == 0:
            raise ValueError('there is not difference between samples!')
        meta = meta[meta_cols]

        meta.to_csv(output[0], sep='\t')


rule comparisons:
    output:
        config['workspace'] + '/aggregate/comparisons.txt'
    input:
        rules.merge_counts.output
    run:
        import pandas as pd

        meta = {sample: meta['meta'] for sample, meta in config['samples'].items()}
        meta = pd.DataFrame.from_dict(meta, orient='index')

        meta_cols = list({
            comparison['condition'] for comparison in config['comparisons'].values()
        })
        if len(meta_cols) == 0:
            raise ValueError('comparisons not defind!')
        meta = meta[meta_cols]

        meta.to_csv(output[0], sep='\t')


rule pcaplot:
    output:
        config['workspace'] + '/aggregate/all_sample_pcaplot.{fmt}',
    input:
        rules.comparisons.output,
        rules.merge_counts.output.qnorm
    params:
        script=os.path.dirname(workflow.snakefile) + '/tools/pcaplot.R'
    shell:
        'Rscript {params.script} {input} {output}'
