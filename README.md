#SplitThreader

Make a graph out of your highly rearranged genome and use it to check for evidence of complex gene fusions, reconstruct the history of variants in a region (such as why an oncogene has become amplified), and check for cycles that could indicate extrachromosomal double minutes. 


##Quick reference guide:
```bash
VARIANTS=Sniffles_SKBR3_Oct28_with_inv_dups.spansplit.bedpe
NODES=Sniffles_SKBR3_Oct28_with_inv_dups.spansplit.nodes.bed

SplitThreader.py check --variants $VARIANTS --nodes $NODES

SplitThreader.py fusions --variants $VARIANTS --nodes $NODES --annotation gencode.v19.annotation.gtf.genes.bed --gene_name_column 8 --list file_with_list_fusions --out test

SplitThreader.py parsimony --variants $VARIANTS --nodes $NODES --out test --chrom 17 --start 37000000 --end 41000000

SplitThreader.py cycles --variants $VARIANTS --nodes $NODES --out test --annotation gencode.v19.annotation.gtf.genes.bed --depth 20
```



##SplitThreader.py check
This is a quick sanity check that will help you answer the following questions:
1. Is the input data in the right format?
    If you encounter issues, make sure you are only using the output of SpanSplit. The *.spansplit.bedpe file is the variants file, and *.spansplit.nodes.bed is the nodes file. 
    - Variants should look like this:
        `X       155032876       155032912       X       155216895       155216915       210     -1      -       +       DUP     24      25      23      21      34      19      25`
    - Nodes should look like this:
        `20      0       213149  1`
2. How big is the graph?
    - Number of nodes and edges is proportional to the number of variants. The size of the graph greatly affects the runtime of cycles, which is currently the only genome-wide program. Parsimony can also be very slow on large graphs if the region selected contains many variants. Check the respective programs for how to deal with this runtime issue. 
3. Does the annotation fit onto the graph?
    Provide annotation using --annotation and give it a bed file with the annotation you wish to use. 
    The annotation file can look something like this:
    `1       11869   14412   ENSG00000223972.4       .       +       pseudogene      DDX11L1`
    The important columns are 1=chromosome, 2=start location, 3=end location, 6=strand, and you can indicate the column with the gene name using --gene_name_column. To pick DDX11L1 as the gene name, do --gene_name_column 8
    You will see the placement rate reported as a percentage. If this is low, check to make sure the chromosome names are consistent. Often this is an issue of one file having a "chr" prefix when the other does not. 



##SplitThreader.py fusions
If you have some evidence that genes may be connected, such as RNA-seq or IsoSeq evidence (PacBio long-read RNA-sequencing), then Fusions will evaluate whether these have genomic evidence in the form of variants that bring the genes close together. 
You provide a list of gene fusions as a file where each line contains two gene names (these must match the column of the annotation set using the --gene_name_column parameter), separated by spaces or tabs.
SplitThreader fusions then evaluates each of these fusions by finding all the paths connecting the two genes in any direction, then reporting the best scoring fusion path according to the following rules:
1. The total length of the fusion gene must end up being less than 1 Mb
2. Paths creating a fusion that reads through both genes in the same direction are favored (so it is possible for both to be transcribed in the correct direction)
3. Paths with a shorter total length are favored. 
4. The path between genes must thread through at most 2 variants (novel adjacencies), 1 is preferred

After finding all the fusions and outputting the result, Fusions also reports which of the gene fusions share the same threading path. This may be useful if your annotation includes many overlapping genes, such as pseudogenes or several isoforms of the same gene, so you can choose not to count the same causal variant in multiple gene fusions. 


##SplitThreader.py parsimony
Oncogene amplification can be caused by a series of variants that copy the oncogene around to different regions of the genome. Parsimony reconstructs some of this history by looking at the likely order in which these variants occurred. Select a region around your favorite oncoge, if possible choose the start and end locations of this region as quiet, non-amplified regions that can act as the baseline.
The output is a bed file showing the genomic coordinates of the sequences in each path. Parsimony also draws a picture showing which nodes are included and excluded in each path through this region. 


##SplitThreader.py cycles
Cycles will identify every loop in the graph, in the hopes that some of them may be actual double minutes. For each cycle, we report the total length, number of variants to thread through, the chromosomes involved, and the minimum read depth on the variants. If you provide annotation through the --annotation parameter, Cycles will also report the genes contained in each cycle. 
A bed file is created with the genomic coordinates of the sequence in each of these cycles. 


