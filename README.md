#SplitThreader

Make a graph out of your highly rearranged genome and use it to check for evidence of complex gene fusions, reconstruct the history of variants in a region (such as why an oncogene has become amplified), and check for cycles that could indicate extrachromosomal double minutes. 


##Quick reference guide:
```bash
Example:

VARIANTS=tests/sniffles.bedpe
GENOME=tests/human_hg19.genome
SplitThreader.py check --variants $VARIANTS --genome $GENOME


VARIANTS=tests/sniffles.bedpe
GENOME=tests/human_hg19.genome
COVERAGE=copy_numbers.segmented.tab
SplitThreader.py Flow --variants $VARIANTS --genome $GENOME --coverage $COVERAGE --out flow_test


VARIANTS=tests/sniffles.bedpe
GENOME=tests/human_hg19.genome
SplitThreader.py Evolution --variants $VARIANTS --genome $GENOME --out evolution_test --chrom 17 --start 37000000 --end 41000000


VARIANTS=tests/sniffles.bedpe
GENOME=tests/human_hg19.genome
LIST=~/Desktop/SplitThreader_testcases/commandline/all_sizes.quivered_hq.fusion_finder.mtc99.mlcbp100.mdbl100kb.bedpe.with_abundance.pair.genes.summary.2_fl_reads.bedpe.list
GENES=~/Desktop/SplitThreader_testcases/commandline/gencode.v19.annotation.gtf.genes.bed
SplitThreader.py Fusions --variants $VARIANTS --genome $GENOME --annotation $GENES --list $LIST --out test

