
# Install bedtools

# Make spansplit.py executable and copy it into a directory in your PATH
# chmod +x bin/spansplit.py

# From the command-line
SplitThreader_prep.sh  variants.bedpe reads.bed output_prefix

python
# In Python:

from SplitThreader import *


g = Graph() 

# From now on access the graph using g.

prefix = "~/Desktop/SplitThreader_results/my_sample"

# Nodes and edges files are output from running SplitThreader_prep.sh
nodes_filename = prefix + ".spansplit.coverage_on_nodes.bed"
edges_filename = prefix + ".spansplit.bedpe"

g.read_spansplit(nodes_filename,edges_filename)


# Points in the genome are tuples of (chromosome, position)
point1 = ("17",37278042)
point2 = ("8",121560937)

# You can find a node by its genomic position
g.find_nodename_by_position(point1)
g.find_nodename_by_position(point2)

# You can find the shortest distance between any 2 points and that shortest path between them
path, distance = self.g.calculate_distance(point1,point2)
print path
print distance


# You can search for where genes are located in the graph by loading in some annotation in bed format
# Format (white-space separated)
# chromosome  start  stop  feature_name  strand  more_columns_are_ignored
annot_filename = "gencode.v19.annotation.gtf.genes.bed"
g.read_annotation(annot_filename,name_field=8)


# You can check what the most likely paths are between any two genes:
report = g.gene_fusion_report("CYTH1","EIF3H",verbose=True)
print report

# You can also create a bed file showing the coordinates of the sequences in a given path
self.g.franken_path(report["path"],output_filename)

## At this point, we can use the bed file using bedtools (outside python) to frankenstein the sequences together into a DNA sequence matching that path in the graph:
# bedtools getfasta -s -fi $REFERENCE -bed $BED_FRANKENSTEIN_PATH -fo $BED_FRANKENSTEIN_PATH.fasta
## Then just collapse those sequences together into a single fasta entry
# cat <(echo ">franken_junction" ) <(cat $BED_FRANKENSTEIN_PATH.fasta | grep -v ">" | tr -d "\n" ) <(echo "") > $BED_FRANKENSTEIN_PATH.unfolded.fasta
## Now you can realign the RNA-seq or IsoSeq data, or even the DNA sequences near those breakpoints, onto the franken-junction to see how well the reads support the junction. 


