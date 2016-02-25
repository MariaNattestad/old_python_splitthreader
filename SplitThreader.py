#!/usr/bin/env python

from SplitThreaderLib import *
import argparse
import time

# Example:
# VARIANTS=tests/sniffles.bedpe
# GENOME=tests/human_hg19.genome
# SplitThreader.py check --variants $VARIANTS --genome $GENOME

# VARIANTS=tests/sniffles.bedpe
# GENOME=tests/human_hg19.genome
# SplitThreader.py Evolution --variants $VARIANTS --genome $GENOME --out evolution_test --chrom 17 --start 37000000 --end 41000000

# VARIANTS=tests/sniffles.bedpe
# GENOME=tests/human_hg19.genome
# LIST=~/Desktop/SplitThreader_testcases/commandline/all_sizes.quivered_hq.fusion_finder.mtc99.mlcbp100.mdbl100kb.bedpe.with_abundance.pair.genes.summary.2_fl_reads.bedpe.list
# GENES=~/Desktop/SplitThreader_testcases/commandline/gencode.v19.annotation.gtf.genes.bed
# SplitThreader.py Fusions --variants $VARIANTS --genome $GENOME --annotation $GENES --list $LIST --out test


def initialize(args):
    variants_filename = args.variants
    genome_file = args.genome_file

    g = Graph()
    g.read_sniffles(variants_filename,genome_file=genome_file)
    g.count_span_vs_split_edges()
    
    return g

def check(args):
    g = initialize(args)
    print "Number of nodes:", len(g.nodes)
    print "Number of edges:", len(g.edges)
    annotation_filename = args.annotation_file
    gene_name_column = args.gene_name_column

    if annotation_filename != None:
        g.read_annotation(annotation_filename,name_field=gene_name_column,by_node_access = True)


def fusions(args):
    annotation_filename = args.annotation_file
    gene_pair_list_file = args.gene_pair_list_file
    gene_name_column = args.gene_name_column
    output_prefix = args.output_prefix

    g = initialize(args)
    g.read_annotation(annotation_filename,name_field=gene_name_column,by_node_access = False)

    g.search_for_fusions(gene_pair_list_file,output_prefix=output_prefix)


def evolution(args):
    chromosome = args.zoom_region_chrom
    start = args.zoom_region_start
    end = args.zoom_region_end
    output_prefix = args.output_prefix
    depth_limit = args.depth_limit

    g = initialize(args)

    g.local_evolution(chromosome,start,end,output_prefix,degree=3,min_weight_required=10,depth_limit=depth_limit)


def flow(args):

    coverage_file = args.coverage
    sniffles_filename = args.variants
    genome_file = args.genome_file
    output_prefix = args.output_prefix
    max_variant_CNV_distance = 100000

    g = initialize(args)

    ############################################################################################################################################################
    # Step 1:  Score split read variants independently of each other and with each breakpoint scored independently of the other
        # Categories:
        # 0: No CNV
        # 1: CNV present only in the opposite direction
        # 2: Multiple CNVs in the correct direction
        # 3: Exactly 1 matching CNV in the correct direction

    SRVs = g.score_SRVs(coverage_file=coverage_file,max_variant_CNV_distance=max_variant_CNV_distance)

    SRVs = g.summarize_SRVs(SRVs)
    print "\n_____________________________"
    print "SRVs by CNV evidence:"
    g.count_by_category(SRVs)

    g.annotate_sniffles_file_with_variant_flow_evaluations(reports = SRVs, sniffles_filename = sniffles_filename, output_filename = output_prefix + ".SRVs.csv") 


    ############################################################################################################################################################
    # Step 2:  Score CNVs
        # Categories:
        # 1: Good variant evidence
        # Many: Messy variant evidence (multiple variants could explain it)
        # 0: No variant evidence

    CNVs_by_SRV_evidence = g.score_CNVs(coverage_file=coverage_file,max_variant_CNV_distance=max_variant_CNV_distance)

    print "\n_____________________________"
    print "CNVs by SRV evidence:"
    g.count_by_category(CNVs_by_SRV_evidence)


    print "\n_____________________________"
    print "CNV features:"

    CNV_graph = Graph()
    CNV_features = CNV_graph.investigate_CNVs_for_quality(coverage_file=coverage_file,genome_file=genome_file,threshold_for_long_CN_segments=100000)

    CNV_graph.count_by_category(CNV_features)


    g.output_CNVs_with_quality_and_SRV_concordance_analysis(coverage_file=coverage_file, CNV_features=CNV_features, CNVs_by_SRV_evidence=CNVs_by_SRV_evidence,output_filename = output_prefix + ".CNVs.tab")


    



def main():
    # Software description
    parser=argparse.ArgumentParser(description="SplitThreader makes a graph out of the genome with sequences represented  \
        by nodes that are split at variant breakpoints and reconnected with the number of reads spanning and split as edge weights")
   
    subparsers = parser.add_subparsers(title="Available commands")
    parser_check = subparsers.add_parser('check', help='Run sanity checks on graph (DO THIS FIRST)')
    parser_fusions = subparsers.add_parser('Fusions', help='Find gene fusions')
    parser_evolution = subparsers.add_parser('Evolution', help='Reconstruct history of a region by finding the most parsimonious set of paths through the graph')
    parser_flow = subparsers.add_parser('Flow', help='Analyze variants for concordance with copy numbers')

    # Parameters for basic checking of the graph
    parser_check.add_argument("--variants",help="bedpe file from running spansplit.py",dest="variants",required=True)
    parser_check.add_argument("--genome",help="genome file where column 1 is chromosome and column 2 is the length of that chromosome",dest="genome_file",required=True)
    parser_check.add_argument("--annotation",help="annotation file",dest="annotation_file",default=None)
    parser_check.add_argument("--gene_name_column",help="which column (1-indexed) of annotation file contains the gene names you want to report?",type=int,dest="gene_name_column",default=8)
    parser_check.set_defaults(func=check)


    #  Parameters for "fusions" program: for finding evidence for gene fusions given annotation and a list of gene pairs
    parser_fusions.add_argument("--variants",help="bedpe file from running spansplit.py",dest="variants",required=True)
    parser_fusions.add_argument("--genome",help="genome file where column 1 is chromosome and column 2 is the length of that chromosome",dest="genome_file",required=True)
    parser_fusions.add_argument("--out",help="prefix for all output files",dest="output_prefix",required=True)
    parser_fusions.add_argument("--annotation",help="annotation file",dest="annotation_file",required=True)
    parser_fusions.add_argument("--list",help="file with list of gene pairs, one pair per line separated by whitespace",dest="gene_pair_list_file",required=True)
    parser_fusions.add_argument("--gene_name_column",help="which column (1-indexed) of annotation file contains the gene names you wish to search by?",type=int,dest="gene_name_column",default=8)
    parser_fusions.set_defaults(func=fusions)
    
    #  Parameters for "evolution" program: Reconstructing history of mutations in a region of the genome
    parser_evolution.add_argument("--variants",help="bedpe file from running spansplit.py",dest="variants",required=True)
    parser_evolution.add_argument("--genome",help="genome file where column 1 is chromosome and column 2 is the length of that chromosome",dest="genome_file",required=True)
    parser_evolution.add_argument("--out",help="prefix for all output files",dest="output_prefix",required=True)
    parser_evolution.add_argument("--chrom",help="which chromosome to focus historical reconstruction on",dest="zoom_region_chrom",type=str,required=True)
    parser_evolution.add_argument("--start",help="which start position within the chromosome to focus historical reconstruction on",dest="zoom_region_start",type=int,required=True)
    parser_evolution.add_argument("--end",help="which end position within the chromosome to focus historical reconstruction on",dest="zoom_region_end",type=int,required=True)
    parser_evolution.add_argument("--depth",help="Number of edges deep to search in the graph. Influences runtime. Default = 30",type=int,dest="depth_limit",default=30)
    parser_evolution.set_defaults(func=evolution)

    parser_flow.add_argument("--variants",help="bedpe file from running spansplit.py",dest="variants",required=True)
    parser_flow.add_argument("--genome",help="genome file where column 1 is chromosome and column 2 is the length of that chromosome",dest="genome_file",required=True)
    parser_flow.add_argument("--coverage",help="coverage bed file with columns: chr\tstart\tend\tunsegmented_coverage\tsegmented_coverage. MUST BE SORTED BY CHROMOSOME THEN BY START POSITION, sort -k1,1 -k2,2n will do this",dest="coverage",required=True)
    parser_flow.add_argument("--out",help="prefix for all output files",dest="output_prefix",required=True)
    parser_flow.set_defaults(func=flow)



    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
    main()














