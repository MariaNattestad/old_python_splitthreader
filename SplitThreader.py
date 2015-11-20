#!/usr/bin/env python
from SplitThreaderLib import *
import argparse
import time


def initialize(args):
    nodes_filename = args.nodes_from_spansplit
    edges_filename = args.variants_from_spansplit


    g = Graph()
    g.read_spansplit(nodes_filename,edges_filename)

    return g


def check(args):
    
    g = initialize(args)
    print "Number of nodes:", len(g.nodes)
    print "Number of edges:", len(g.edges)
    

def fusions(args):
    print "Running fusions"

    annotation_filename = args.annotation_file
    gene_pair_list_file = args.gene_pair_list_file
    gene_name_column = args.gene_name_column

    g = initialize(args)
    g.read_annotation(annotation_filename,name_field=gene_name_column)

    f=open(gene_pair_list_file)
    print "gene_1\tgene2\tnumber_of_variants_to_thread_through\tbp_distance\tstrand_gene_1\tstrand_gene_2"

    for line in f:
        fields = line.strip().split()
        report =  g.gene_fusion_report(fields[0],fields[1])
        if report == None:

            print "%s\t%s\t%s\t%s\t%s\t%s" % (fields[0],fields[1], "none","none","none","none")
        else:
            print "%s\t%s\t%d\t%d\t%s\t%s" % (fields[0],fields[1], report["number_of_splits"], report["distance"],report["Gene1_direction"],report["Gene2_direction"])
            # self.g.franken_path(report["path"],"/Users/mnattest/Desktop/SplitThreader_testcases/Extra_Sniffles_franken_%s_%s.txt" % (fields[0],fields[1]))
        
    f.close()




def parsimony(args):
    region_chromosome = args.zoom_region_chrom
    print region_chromosome

    g = initialize(args)


def main():

    # Software description
    parser=argparse.ArgumentParser(description="SplitThreader makes a graph out of the genome with sequences represented  \
        by nodes that are split at variant breakpoints and reconnected with the number of reads spanning and split as edge weights")
    
   
    subparsers = parser.add_subparsers(title="Available commands")
    parser_check = subparsers.add_parser('check', help='Run sanity checks on graph (DO THIS FIRST)')
    parser_fusions = subparsers.add_parser('fusions', help='Find gene fusions')
    parser_parsimony = subparsers.add_parser('parsimony', help='Reconstruct history of a region by finding the most parsimonious set of paths through the graph')

    # Parameters needed for all programs
    # parser.add_argument("--variants",help="bedpe file from running spansplit.py",dest="variants_from_spansplit",required=True) # Variant calls in bedpe format (like from Sniffles)
    # parser.add_argument("--nodes",help="nodes file from running spansplit.py",dest="nodes_from_spansplit",required=True)


    # Parameters for basic checking of the graph
    
    parser_check.set_defaults(func=check)
    parser_check.add_argument("--variants",help="bedpe file from running spansplit.py",dest="variants_from_spansplit",required=True) # Variant calls in bedpe format (like from Sniffles)
    parser_check.add_argument("--nodes",help="nodes file from running spansplit.py",dest="nodes_from_spansplit",required=True)


    #  Parameters for "fusions" program: for finding evidence for gene fusions given annotation and a list of gene pairs
    parser_fusions.add_argument("--variants",help="bedpe file from running spansplit.py",dest="variants_from_spansplit",required=True) # Variant calls in bedpe format (like from Sniffles)
    parser_fusions.add_argument("--nodes",help="nodes file from running spansplit.py",dest="nodes_from_spansplit",required=True)
    parser_fusions.add_argument("--annotation",help="annotation file",dest="annotation_file",required=True)
    parser_fusions.add_argument("--list",help="file with list of gene pairs, one pair per line separated by whitespace",dest="gene_pair_list_file",required=True)
    parser_fusions.add_argument("--gene_name_column",help="which column (1-indexed) of annotation file contains the gene names you wish to search by?",type=int,dest="gene_name_column",required=True)
    parser_fusions.set_defaults(func=fusions)
    
    #  Parameters for "parsimony" program: Reconstructing history of mutations in a region of the genome
    parser_parsimony.add_argument("--variants",help="bedpe file from running spansplit.py",dest="variants_from_spansplit",required=True) # Variant calls in bedpe format (like from Sniffles)
    parser_parsimony.add_argument("--nodes",help="nodes file from running spansplit.py",dest="nodes_from_spansplit",required=True)
    parser_parsimony.add_argument("--chrom",help="which chromosome to focus historical reconstruction on",dest="zoom_region_chrom",required=True)
    parser_parsimony.set_defaults(func=parsimony)
    

    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
    main()

