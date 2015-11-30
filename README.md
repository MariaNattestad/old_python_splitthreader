#SplitThreader

Make a graph out of your highly rearranged genome and use it to check for evidence of complex gene fusions, reconstruct the history of variants in a region (such as why an oncogene has become amplified), and check for cycles that could indicate extrachromosomal double minutes. 


##Example:
```bash
VARIANTS=Sniffles_SKBR3_Oct28_with_inv_dups.spansplit.bedpe
NODES=Sniffles_SKBR3_Oct28_with_inv_dups.spansplit.nodes.bed
SplitThreader.py check --variants $VARIANTS --nodes $NODES
SplitThreader.py fusions --variants $VARIANTS --nodes $NODES --annotation gencode.v19.annotation.gtf.genes.bed --list IsoSeq_fusions_with_5_reads_but_no_direct_Sniffles --out test
SplitThreader.py parsimony --variants $VARIANTS --nodes $NODES --out test --chrom 17 --start 37000000 --end 41000000
SplitThreader.py cycles --variants $VARIANTS --nodes $NODES --out test --annotation gencode.v19.annotation.gtf.genes.bed --depth 20
```


