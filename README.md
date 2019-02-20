# pyCoverage
Count reads in each intervals from a bed file. Reads are stored in bam file. pyCoverage is similar to coverageBed with `-counts` and `-b bam_file`.
But pyCoverage can use multiple process and the run very fast.

Run pyCoverage using 8 threads:

```{bash}
pyCoverage.py  bam_sorted.bam interval.bed 8 > interval.cov
```

The counts are added to the end of each line in `interval.bed`.

pyCoverage.py accept stdin:
```{bash}
awk '$1==4' interval.bed | pyCoverage.py  bam_sorted.bam - 16 > interval_chr4.cov
```

Bam file must be sorted and indexed.

Which means pyCoverage will use 40 threads if the threads number is not specified:

```{bash}
pyCoverage.py  bam_sorted.bam interval.bed > interval.cov
```
will use 40 threads.
