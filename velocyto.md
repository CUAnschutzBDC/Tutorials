# Testing velocyto

## The working run
Here is what my directory looked like before running velocyto
```{bash}
(veloctyo) [wellskri@compute21 test]$ ls
```

```
barcodes.tsv  cellsorted_possorted_genome_bam.bam  possorted_genome_bam.bam
```

Here is my command
```{bash}
(veloctyo) [wellskri@compute21 test]$ velocyto run -b barcodes.tsv -o velocyto -vv possorted_genome_bam.bam /beevol/home/rbilab/ref/cellranger/mouse/refdata-gex-mm10-2020-A/genes/genes.gtf
```

The `-vv` is just to show me the info and debugging output from velocyto

From the debugging, one line tells me that sorting will be skipped
```
2021-04-01 18:34:39,262 - WARNING - The file /beevol/home/wellskri/Analysis/test/cellsorted_possorted_genome_bam.bam already exists. The sorting step will be skipped and the existing file will be used.
```

## Trouble shooting
I also tried with your directory structure:

```{bash}
(veloctyo) [wellskri@compute21 test]$ ls
```

```
barcodes.tsv   WT_cellsorted_possorted_genome_bam.bam  WT_possorted_genome_bam.bam
```

Here I have `WT_cellsorted_possorted_genome_bam.bam` and `WT_possorted_genome_bam.bam`

When I run this the file tries to sort because it is looking for `cellsorted_WT_possorted_genome_bam.bam`

```{bash}
(veloctyo) [wellskri@compute21 test]$ velocyto run -b barcodes.tsv -o velocyto -vv WT_possorted_genome_bam.bam /beevol/home/rbilab/ref/cellranger/mouse/refdata-gex-mm10-2020-A/genes/genes.gtf
```

Here is the line that tells me the file is being sorted
```
2021-04-01 18:45:31,383 - INFO - Starting the sorting process of /beevol/home/wellskri/Analysis/test/WT_possorted_genome_bam.bam the output will be at: /beevol/home/wellskri/Analysis/test/cellsorted_WT_possorted_genome_bam.bam
2021-04-01 18:45:31,384 - INFO - Command being run is: samtools sort -l 7 -m 2048M -t CB -O BAM -@ 16 -o /beevol/home/wellskri/Analysis/test/cellsorted_WT_possorted_genome_bam.bam /beevol/home/wellskri/Analysis/test/WT_possorted_genome_bam.bam
2021-04-01 18:45:31,384 - INFO - While the bam sorting happens do other things...
```

## Fixing the problem
To solve this issue, you will need to rename your file. The file that was successful in samtools sort you called `WT_cellsorted_possorted_genome_bam.bam`. You need to rename it to `cellsorted_WT_possorted_genome_bam.bam`

You can do that by running the following command
```{bash}
(veloctyo) [wellskri@compute21 test]$ mv WT_cellsorted_possorted_genome_bam.bam cellsorted_WT_possorted_genome_bam.bam
(veloctyo) [wellskri@compute21 test]$ ls
```

```
barcodes.tsv  cellsorted_WT_possorted_genome_bam.bam  WT_possorted_genome_bam.bam
```

Now running the same command as before
```{bash}
(veloctyo) [wellskri@compute21 test]$ velocyto run -b barcodes.tsv -o velocyto -vv WT_possorted_genome_bam.bam /beevol/home/rbilab/ref/cellranger/mouse/refdata-gex-mm10-2020-A/genes/genes.gtf
```

I again see the sorting step is skipped
```
2021-04-01 18:51:03,801 - WARNING - The file /beevol/home/wellskri/Analysis/test/cellsorted_WT_possorted_genome_bam.bam already exists. The sorting step will be skipped and the existing file will be used.
```