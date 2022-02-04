# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: qc.smk
# time: 2022/01/04


rule promoters:
    output:
        config['workspace'] + '/share/promoters.bed'
    input:
        config['genome']['txdb']
    params:
        script=lambda wildcards: BASE_DIR + '/tools/promoters.R'
    shell:
        'Rscript {params.script} {input} {output}'


rule fragmentsizes:
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


rule qc:
    output:
        config['workspace'] + '/samples/{sample}/qc/align_summary.txt'
    input:
        bam=rules.merge_bam.output,
        index=rules.index.output,
        bed=rules.clean_bed.output,
        promoters=rules.promoters.output
    threads: 8 if workflow.cores >= 8 else workflow.cores
    run:
        from multiprocessing import Pool
        import pysam

        pool = Pool(processes=8)
        total = pool.apply_async(pysam.view, ('-c', f'{input.bam}'))
        mapped = pool.apply_async(pysam.view, ('-c', '-F', '4', f'{input.bam}'))
        dup = pool.apply_async(pysam.view, ('-c', '-f', '1024', f'{input.bam}'))
        chrm = pool.apply_async(pysam.view, ('-c', f'{input.bam}', 'chrM'))
        clean = pool.apply_async(
            pysam.view, (
                '-c', '-f', '2', '-F', '1028', '-q', '10', '-L', f'{input.bed}', f'{input.bam}'
            )
        )
        promoters = pool.apply_async(
            pysam.view, (
                '-c', '-f', '2', '-F', '1028', '-q', '10', '-L', f'{input.promoters}', f'{input.bam}'
            )
        )
        pool.close()
        pool.join()
        # collect results
        total = int(total.get().strip())
        mapped = int(mapped.get().strip())
        dup = int(dup.get().strip())
        chrm = int(chrm.get().strip())
        clean = int(clean.get().strip())
        promoters = int(promoters.get().strip())

        with open(output[0], 'w') as f:
            f.write(f'total_reads\t{total}\t{total/total:.2%}\n')
            f.write(f'mapped_reads\t{mapped}\t{mapped/total:.2%}\n')
            f.write(f'dup_reads\t{dup}\t{dup/total:.2%}\n')
            f.write(f'chrM_reads\t{chrm}\t{chrm/total:.2%}\n')
            f.write(f'clean_reads\t{clean}\t{clean/total:.2%}\n')
            f.write(f'promoters_reads\t{promoters}\t{promoters/total:.2%}\n')
