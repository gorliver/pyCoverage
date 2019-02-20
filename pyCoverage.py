#!/usr/bin/env python3

from __future__ import print_function
import pysam
import random
import multiprocessing
import sys


def check_read(read, min_qual=20):
    """ filter reads, default is mapq >=20. Used in pysam callback """
    if read.mapq >= min_qual:
        return True
    return False


def count_wrapper(args):
    """ wrapper for multiprocessing """
    return countReads.count_reads(*args)


def get_handle(file):
    """ check input bed file is stdin or a file """
    d = []
    if file == "-":
        b = sys.stdin
    else:
        b = open(file, 'r')

    d = read_data(b)
    b.close()
    return d


def read_data(file_handle):
    """ read data from bed file """
    data = []
    for line in file_handle:
        arr = line.strip().split("\t")
        data.append(arr)
    return data


class countReads(object):
    def __init__(self, bam_file=None, bed_file=None, threads=8):
        self.bam_file = bam_file
        self.bed_file = bed_file
        self.bed_regions = None
        self.threads = threads

    def count_reads(self, chrom, st, ed):
        """ count reads in the interval using AlignmentFile.count """
        bam_file = self.bam_file
        bam = pysam.AlignmentFile(bam_file, mode='rb')

        count = bam.count(
            contig=str(chrom),
            start=int(st),
            end=int(ed),
            read_callback=check_read)
        #        print("got ", count)
        return (chrom, st, ed, count)

    def get_bed_regions(self):
        """ get intervals from bed file or stdin """
        if self.bed_regions is not None:
            return self.bed_regions
        bed_regions = []
        #        bed_file_handle = get_handle(self.bed_file)
        #        bed_regions = read_data(bed_file_handle)
        bed_regions = get_handle(self.bed_file)
        self.bed_regions = bed_regions
        #        bed_file_handle.close()
        #        print("file name is ",self.bed_file)

        #        print(bed_regions[0:5])
        return bed_regions

    def work_splitter(self, func, regions=None, threads=1):
        """ split count works """
        args = [[self, x[0], int(x[1]), int(x[2])] for x in regions]
        random.shuffle(regions)
        pool = multiprocessing.Pool(threads)
        out = pool.map_async(func, args).get(999999)
        pool.close()
        pool.join()
        return out

    def get_counts(self):
        """ count number of reads in each intervals, store the results in a
        dict with key is the position and value is the counts """
        bed_regions = self.get_bed_regions()
        result = self.work_splitter(count_wrapper, bed_regions, self.threads)

        count_dict = dict()
        for i in result:
            pos = i[0] + ":" + str(i[1]) + "-" + str(i[2])
            count_dict[pos] = i[3]
        return count_dict


if __name__ == "__main__":

    bam_file = ''
    bed_file = ''
    threads = 40
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        print("python3 {} <bam_file> <bed_file> <threads:40>".format(
            sys.argv[0]))
        exit()

    if len(sys.argv) == 4:
        threads = sys.argv[3]

    bam_file = sys.argv[1]
    bed_file = sys.argv[2]

    cr = countReads(bam_file, bed_file, threads)

    results = cr.get_counts()

    intervals = cr.get_bed_regions()
    intervals.sort(key=lambda x: (str(x[0]), int(x[1])))
    for i in intervals:
        pos = i[0] + ":" + str(i[1]) + "-" + str(i[2])
        count = results[pos]
        i.append(count)
        print("\t".join(str(x) for x in i))
