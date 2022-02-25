# -*- coding:utf-8 -*-
#
# author: yangmqglobe
# file: fragmentsizes.py
# time: 2021/02/16
from argparse import ArgumentParser
from collections import Counter
import json

import pysam


def read_regions(filename):
    with open(filename) as f:
        for line in f:
            chrom, start, end = line.strip('\n').split('\t')[:3]
            start, end = int(start), int(end)
            yield chrom, start, end


def fetch_reads(bam, regions=None):
    if regions is None:
        yield from bam
    else:
        regions = read_regions(regions)
        for region in regions:
            for read in bam.fetch(*region):
                yield read


def fetch_fragment_lens(in_bam,
                        regions=None,
                        include=2,
                        exclude=1024,
                        mapq=10):
    with pysam.AlignmentFile(in_bam, 'rb') as bam:
        for read in fetch_reads(bam, regions):
            if not read.flag & include:
                continue
            if read.flag & exclude:
                continue
            if read.mapq < mapq:
                continue
            yield abs(read.template_length)


def run(in_bam, out, regions=None, include=2, exclude=1024, mapq=10):
    counts = Counter()
    for fl in fetch_fragment_lens(in_bam, regions, include, exclude, mapq):
        counts[fl] += 1
    counts = [{'size': size, 'count': count} for size, count in counts.items()]
    counts = sorted(counts, key=lambda x: x['size'])

    # output
    with open(out, 'w') as f:
        json.dump({'fragment_sizes': counts}, f)


def main():
    parser = ArgumentParser(description='Count fragment size')
    parser.add_argument('-L',
                        dest='regions',
                        help='only include reads overlapping this BED FILE',
                        default=None)
    parser.add_argument(
        '-f',
        dest='include',
        help='only include reads with all  of the FLAGs in INT present',
        default=2,
        type=int)
    parser.add_argument(
        '-F',
        dest='exclude',
        help='only include reads with none of the FLAGS in INT present',
        default=1028,
        type=int)
    parser.add_argument('-q',
                        dest='mapq',
                        help='only include reads with mapping quality >= INT',
                        default=10,
                        type=int)
    parser.add_argument('<in_bam>', help='input bam file')
    parser.add_argument('<out>', help='output json file')
    args = vars(parser.parse_args())
    run(args['<in_bam>'], args['<out>'], args['regions'], args['include'],
        args['exclude'], args['mapq'])


if __name__ == '__main__':
    main()
