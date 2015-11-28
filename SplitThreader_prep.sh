


# Variant calls in bedpe format, from Sniffles
BEDPE=

# Reads must be long reads because the spansplit.py program only counts reads that map all the way across an interval of size SLOP centered at the breakpoint
READ_BED=

# Output prefix:
NAME=

# File where each line has a chromosome name (column 1) and the length of the chromosome (column 2)
GENOME_FILE=

# Example inputs:
# READ_BED=/seq/schatz/mnattest/hg19_skbr3/analysis/bwamem.hg19.position_sorted.bed
# BEDPE=/seq/schatz/SKBR3/Sniffles/SKBR3_Sniffles_10.19.primary_chrom.over10kb.merge_within_1kb.darkzone_filtered.bedpe
# NAME=/seq/schatz/SKBR3/SpanSplit/SKBR3_Sniffles_10.19
# GENOME_FILE=/seq/schatz/SKBR3/reference/genome_files_for_bedtools/human_g1k_v37.fasta.genome

SLOP=1000
MQ=60


#############################################################################

# Standardize variant calls and expand intervals so next we can grab all the reads aligning at the breakpoints
awk -v slop=$SLOP '{mid1=int(($3+$2)/2); mid2=int(($5+$6)/2); $2=mid1-slop*1.5; $3=mid1+slop*1.5; $5=mid2-slop*1.5; $6=mid2+slop*1.5; if($2<0){$2=0}; if($3<0){$3=0}; if($5<0){$5=0}; if($6<0){$6=0}; print}' OFS="\t" $BEDPE > $NAME.calls.expanded.bedpe 
awk '{print $1,$2,$3,$7,1; print $4,$5,$6,$7,2}' OFS="\t" $NAME.calls.expanded.bedpe > $NAME.calls.encoding.bed

# Collect reads and encode them with which breakpoint they are located nearby
time bedtools intersect -a $READ_BED -b $NAME.calls.encoding.bed -wb > $NAME.reads_near_breakpoints.bed

# Count reads spanning and split
spansplit.py -bedpe $BEDPE -reads $NAME.reads_near_breakpoints.bed -out $NAME -mq $MQ -slop $SLOP -genome-file $GENOME_FILE


# Get coverage on each node
time bedtools intersect -wb -a $READ_BED -b $NAME.spansplit.nodes.bed > $NAME.spansplit.reads_on_nodes.bed 

time sort -k7,7 -k8,8n $NAME.spansplit.reads_on_nodes.bed > $NAME.spansplit.reads_on_nodes.sorted.bed

time awk '{if(name!=$10 && length(name)!=0){print chrom,pos,pos2,name,sum/len;sum=0};name=$10;chrom=$7;pos=$8;pos2=$9;name=$10;len=($9-$8);sum+=($3-$2)}END{print chrom,pos,pos2,name,sum/len}' OFS="\t" $NAME.spansplit.reads_on_nodes.sorted.bed > $NAME.spansplit.coverage_on_nodes.bed


#output:

#chrom1 start1 stop1   chrom2 start2 stop2  variantID score strand1 strand2 typeOfSV numreads left1 span1 right1 left2 span2 right2


