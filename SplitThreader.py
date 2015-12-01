#!/usr/bin/env python
from SplitThreaderLib import *
import argparse
import time

# Example:
# VARIANTS=Sniffles_SKBR3_Oct28_with_inv_dups.spansplit.bedpe
# NODES=Sniffles_SKBR3_Oct28_with_inv_dups.spansplit.nodes.bed
# SplitThreader.py check --variants $VARIANTS --nodes $NODES
# SplitThreader.py fusions --variants $VARIANTS --nodes $NODES --annotation gencode.v19.annotation.gtf.genes.bed --list IsoSeq_fusions_with_5_reads_but_no_direct_Sniffles --out test
# SplitThreader.py parsimony --variants $VARIANTS --nodes $NODES --out test --chrom 17 --start 37000000 --end 41000000
# SplitThreader.py cycles --variants $VARIANTS --nodes $NODES --out test --annotation gencode.v19.annotation.gtf.genes.bed --depth 20

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
    annotation_filename = args.annotation_file
    gene_name_column = args.gene_name_column

    if annotation_filename != None:
        g.read_annotation(annotation_filename,name_field=gene_name_column,by_node_access = True)

    

def fusions(args):
    annotation_filename = args.annotation_file
    gene_pair_list_file = args.gene_pair_list_file
    gene_name_column = args.gene_name_column
    output_fusion_report = args.output_prefix + ".fusion_reports.txt"

    g = initialize(args)
    g.read_annotation(annotation_filename,name_field=gene_name_column,by_node_access = False)

    f=open(gene_pair_list_file)
    
    f_output_fusion_report = open(output_fusion_report,"w")
    f_output_fusion_report.write("gene_1\tgene2\tnumber_of_variants_to_thread_through\ttotal_fusion_gene_length\tstrand_gene_1\tstrand_gene_2\tminimum_read_depth_on_breakpoints\n")

    gene_fusion_reports = []
    gene_names = []
    for line in f:
        fields = line.strip().split()
        gene_names.append((fields[0],fields[1]))
        gene_fusion_reports.append(g.gene_fusion_report(fields[0],fields[1]))

    total_putative_fusions = len(gene_fusion_reports)
    total_supported_fusions = 0
    
    all_supported_fusion_paths = []
    all_supported_fusion_names = []

    for i in xrange(len(gene_fusion_reports)):
        report = gene_fusion_reports[i]
        if report == None:
            f_output_fusion_report.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (gene_names[i][0],gene_names[i][1], "none","none","none","none","none"))
        else:
            total_supported_fusions += 1
            all_supported_fusion_paths.append(report["path"])
            all_supported_fusion_names.append(report["Gene1"] + "-" + report["Gene2"])
            f_output_fusion_report.write("%s\t%s\t%d\t%d\t%s\t%s\t%d\n" % (gene_names[i][0],gene_names[i][1], report["number_of_splits"], report["distance"],report["Gene1_direction"],report["Gene2_direction"],g.min_weight(report["path"])))
            # self.g.franken_path(report["path"],"/Users/mnattest/Desktop/SplitThreader_testcases/Extra_Sniffles_franken_%s_%s.txt" % (fields[0],fields[1]))
    f_output_fusion_report.close()
    f.close()
    print "Found evidence of %d gene fusions out of the list of %d given:" % (total_supported_fusions,total_putative_fusions)


    print "%d of these are through unique paths" % (len(g.collapse_redundant_paths(all_supported_fusion_paths,split_edges_only = True,names=all_supported_fusion_names)))

    print "Wrote output to %s" % (output_fusion_report)


def parsimony(args):
    chromosome = args.zoom_region_chrom
    start = args.zoom_region_start
    end = args.zoom_region_end
    output_prefix = args.output_prefix

    g = initialize(args)

    g.local_parsimony(chromosome,start,end,output_prefix)

    # s = g.subgraph_from_genome_interval(chromosome,start,end)
    # print "Number of nodes:", len(s.nodes)
    # print "Number of edges:", len(s.edges)
    # # print s.nodes
    # # for edge in s.edges:
    # #     print edge,edge.weight
    # # s.to_json(output_prefix+".subgraph.json")

    # recordings = s.parsimony(depth_limit = 30)
    # print "Parsimonious set of", len(recordings), "paths found"
    # for report in recordings:
    #     print report
    # s.karyotype_from_parsimony(recordings,output_filename=output_prefix)

def cycles(args):
    g = initialize(args)
    output_cycles = args.output_prefix + ".cycles_summary.txt"
    # output_json = args.output_prefix + ".cycles.json"
    output_bed = args.output_prefix + ".cycles.bed"
    depth_limit = args.depth_limit
    annotation_filename = args.annotation_file
    gene_name_column = args.gene_name_column

    if annotation_filename != None:
        g.read_annotation(annotation_filename,name_field=gene_name_column,by_node_access = True)

    cycles = g.find_cycles(depth_limit=depth_limit)
    print "%d cycles found at depth limit %d" % (len(cycles),depth_limit)

    f=open(output_cycles,"w")
    f.write("path_ID\tnumber_of_splits\tminimum_read_depth_on_breakpoints\tcycle_length_in_bp\tchromosomes_involved\tgenes_in_cycle\n")
    path_counter = 0
    for cycle in cycles:
        path = cycle[:-1] # remove the last element because it is the same as the first
        chromosomes = set()
        for item in path:
            chromosomes.add(g.port_from_path_item(item).node.attributes["chrom"])
        path_counter += 1
        f.write("path_%03d\t%d\t%d\t%d\t" % (path_counter,g.count_splits_in_path(cycle), g.min_weight(path), g.find_total_length(path)) + ",".join(chromosomes) + "\t" + ",".join(g.genes_on_path(path)) + "\n")

    f.close()
    g.franken_paths(cycles,output_bed)
    print "Wrote output to %s and %s:" % (output_cycles,output_bed)
    print "\tSummary information on each cycle: %s" % (output_cycles)
    print "\tGenomic coordinates of the sequences involved, in bed format: %s" % (output_bed)

    # g.paths_to_json(cycles,output_json)

def main():

    # Software description
    parser=argparse.ArgumentParser(description="SplitThreader makes a graph out of the genome with sequences represented  \
        by nodes that are split at variant breakpoints and reconnected with the number of reads spanning and split as edge weights")
    
   
    subparsers = parser.add_subparsers(title="Available commands")
    parser_check = subparsers.add_parser('check', help='Run sanity checks on graph (DO THIS FIRST)')
    parser_fusions = subparsers.add_parser('fusions', help='Find gene fusions')
    parser_parsimony = subparsers.add_parser('parsimony', help='Reconstruct history of a region by finding the most parsimonious set of paths through the graph')
    parser_cycles = subparsers.add_parser('cycles', help='Detects cycles in the graph')

    # Parameters needed for all programs
    # parser.add_argument("--variants",help="bedpe file from running spansplit.py",dest="variants_from_spansplit",required=True) # Variant calls in bedpe format (like from Sniffles)
    # parser.add_argument("--nodes",help="nodes file from running spansplit.py",dest="nodes_from_spansplit",required=True)


    # Parameters for basic checking of the graph
    parser_check.set_defaults(func=check)
    parser_check.add_argument("--variants",help="bedpe file from running spansplit.py",dest="variants_from_spansplit",required=True) # Variant calls in bedpe format (like from Sniffles)
    parser_check.add_argument("--nodes",help="nodes file from running spansplit.py",dest="nodes_from_spansplit",required=True)
    parser_check.add_argument("--annotation",help="annotation file",dest="annotation_file",default=None)
    parser_check.add_argument("--gene_name_column",help="which column (1-indexed) of annotation file contains the gene names you want to report?",type=int,dest="gene_name_column",default=8)


    #  Parameters for "fusions" program: for finding evidence for gene fusions given annotation and a list of gene pairs
    parser_fusions.add_argument("--variants",help="bedpe file from running spansplit.py",dest="variants_from_spansplit",required=True) # Variant calls in bedpe format (like from Sniffles)
    parser_fusions.add_argument("--nodes",help="nodes file from running spansplit.py",dest="nodes_from_spansplit",required=True)
    parser_fusions.add_argument("--out",help="prefix for all output files",dest="output_prefix",required=True)
    parser_fusions.add_argument("--annotation",help="annotation file",dest="annotation_file",required=True)
    parser_fusions.add_argument("--list",help="file with list of gene pairs, one pair per line separated by whitespace",dest="gene_pair_list_file",required=True)
    parser_fusions.add_argument("--gene_name_column",help="which column (1-indexed) of annotation file contains the gene names you wish to search by?",type=int,dest="gene_name_column",default=8)
    parser_fusions.set_defaults(func=fusions)
    
    #  Parameters for "parsimony" program: Reconstructing history of mutations in a region of the genome
    parser_parsimony.add_argument("--variants",help="bedpe file from running spansplit.py",dest="variants_from_spansplit",required=True) # Variant calls in bedpe format (like from Sniffles)
    parser_parsimony.add_argument("--nodes",help="nodes file from running spansplit.py",dest="nodes_from_spansplit",required=True)
    parser_parsimony.add_argument("--out",help="prefix for all output files",dest="output_prefix",required=True)
    parser_parsimony.add_argument("--chrom",help="which chromosome to focus historical reconstruction on",dest="zoom_region_chrom",type=str,required=True)
    parser_parsimony.add_argument("--start",help="which start position within the chromosome to focus historical reconstruction on",dest="zoom_region_start",type=int,required=True)
    parser_parsimony.add_argument("--end",help="which end position within the chromosome to focus historical reconstruction on",dest="zoom_region_end",type=int,required=True)
    parser_parsimony.set_defaults(func=parsimony)

    #  Parameters for "cycles" program: cycle detection outputting paths
    parser_cycles.add_argument("--variants",help="bedpe file from running spansplit.py",dest="variants_from_spansplit",required=True) # Variant calls in bedpe format (like from Sniffles)
    parser_cycles.add_argument("--nodes",help="nodes file from running spansplit.py",dest="nodes_from_spansplit",required=True)
    parser_cycles.add_argument("--out",help="prefix for all output files",dest="output_prefix",required=True)
    parser_cycles.add_argument("--annotation",help="annotation file",dest="annotation_file",default=None)
    parser_cycles.add_argument("--gene_name_column",help="which column (1-indexed) of annotation file contains the gene names you want to report?",type=int,dest="gene_name_column",default=8)
    parser_cycles.add_argument("--depth",help="Number of edges deep to search in the graph. Influences runtime. Default = 20",type=int,dest="depth_limit",default=20)
    parser_cycles.set_defaults(func=cycles)
    

    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
    main()

