# This is a very old version of SplitThreader that was initially coded in Python

I made this repo public only as a fun historical artifact. Please see [splitthreader.com](splitthreader.com) and the [actual SplitThreader repo](https://github.com/MariaNattestad/SplitThreader) for the version of SplitThreader you would ever want to use.


## Old README usage information.
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
LIST=~/Desktop/SplitThreader_testcases/commandline/all_sizes.quivered_hq.fusion_finder.mtc99.mlcbp100.mdbl100kb.bedpe.with_abundance.pair.genes.summary.2_fl_reads.bedpe.list
GENES=~/Desktop/SplitThreader_testcases/commandline/gencode.v19.annotation.gtf.genes.bed
SplitThreader.py Fusions --variants $VARIANTS --genome $GENOME --annotation $GENES --list $LIST --out test


VARIANTS=tests/sniffles.bedpe
GENOME=tests/human_hg19.genome
SplitThreader.py Evolution --variants $VARIANTS --genome $GENOME --out evolution_test --chrom 17 --start 37000000 --end 41000000



################################    Creating copy number input for SplitThreader    ################################

# Just do this filtering once:
BAM=analysis/bwamem.hg19.position_sorted.bam
samtools view -b -q 60 $BAM > ${BAM%.bam}.mq60.bam

BAM=${BAM%.bam}.mq60.bam
BASE=coverage/mq60/SKBR3_bwamem_hg19.position_sorted.mq60

# getting coverage for every basepair in the genome
## bedtools genomecov -d -ibam $BAM -g /seq/schatz/mnattest/reference_genomes/human/hg19.genome > $BASE.coverage

############### 1 kb bins ################
# making 1kb bins on the coverage
### awk 'BEGIN{num=1000;sum=0;possum=0}{if(NR!=1){if(possum==1000 || chrom!=$1){print chrom,pos,sum/possum,possum;sum=0;possum=0;}}

# Create a .csv file for the SplitThreader web app:
### awk 'BEGIN{start=0;print "chromosome,start,end,unsegmented_coverage"}{if(chrom!=$1){start=0}; print $1,start,$2,$3;start=$2;chr


############## 10 kb bins ################
# making 10kb bins on the coverage
awk 'BEGIN{num=10000;sum=0;possum=0}{if(NR!=1){if(possum==10000 || chrom!=$1){print chrom,pos,sum/possum,possum;sum=0;possum=0;}}{s

# Create a .csv file for the SplitThreader web app:
awk 'BEGIN{start=0;print "chromosome,start,end,unsegmented_coverage"}{if(chrom!=$1){start=0}; print $1,start,$2,$3;start=$2;chrom=$

####################################################################################################################
```
