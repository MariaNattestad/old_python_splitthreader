#! /usr/bin/env python
import argparse
import numpy as np

def run(args):
    print "running"
    interval_size = args.slop

    bedpe_file = args.bedpe  
    reads_file = args.reads  
    output_file = args.out + ".spansplit.bedpe"
    slop = args.slop

    mapping_quality = args.mq
    number_of_bed_file_columns = args.offset

    # format of $NAME.variant.expanded.bedpe
    # ####  side 1 ####  #### side 2 ####
    # chrom1 start1 stop1   chrom2 start2 stop2  variantID score strand1 strand2 typeOfSV IDS:*,numreads


    # format of $NAME.reads_near_breakpoints.bed
    # ######## read info #######  ######### variant info ##### 1 for first in bedpe file, 2 for second in bedpe file
    # chrom start stop readname   chrom start stop variantID 1/2

    left_dict = {}
    span_dict = {}
    right_dict = {}
    count_pass_mapping_quality = 0
    count_fail_mapping_quality = 0
    count_on_edge = 0
    readnames_1 = {}
    readnames_2 = {}

    f=open(reads_file)
    for line in f:
        neat = line.strip().split()
        readchrom = neat[0]
        readstart = int(neat[1])
        readstop = int(neat[2])
        readname = neat[3]
        read_mapping_quality = float(neat[4])
        if read_mapping_quality >= mapping_quality:
            count_pass_mapping_quality += 1
            variantchrom = neat[number_of_bed_file_columns]
            variantstart = int(neat[number_of_bed_file_columns+1])
            variantstop = int(neat[number_of_bed_file_columns+2])
            variantID = neat[number_of_bed_file_columns+3]
            side = int(neat[number_of_bed_file_columns+4])
            if side == 1:
                readnames_1[variantID] = readnames_1.get(variantID,set()).union(set([readname]))
            elif side == 2:
                readnames_2[variantID] = readnames_2.get(variantID,set()).union(set([readname]))

            key = (variantID,side)
            if variantstop-variantstart != 3*interval_size:
                #print variantstart,variantstop,variantstop-variantstart
        		count_on_edge += 1
        		#continue
            # check if read spans left side
            if readstart <= variantstart and readstop >= variantstart+interval_size:
                left_dict[key] = left_dict.get(key,0) + 1
            # check if read spans breakpoint interval
            if readstart <= variantstart+interval_size and readstop >= variantstop-interval_size:
                span_dict[key] = span_dict.get(key,0) + 1
            # check if read spans right side
            if readstart <= variantstop-interval_size and readstop >= variantstop:
                right_dict[key] = right_dict.get(key,0) + 1
        else:
            count_fail_mapping_quality += 1

    f.close()
    print "len(left_dict): %d" % (len(left_dict))
    print "len(span_dict): %d" % (len(span_dict))
    print "len(right_dict): %d" % (len(right_dict))
    print "Reads passing mq test: %d" % count_pass_mapping_quality
    print "Reads failing mq test: %d" % count_fail_mapping_quality
    print "intervals on edge: %d" % count_on_edge

    
    # ####  side 1 ####      #### side 2 ####
    # chrom1 start1 stop1   chrom2 start2 stop2  variantID score strand1 strand2 typeOfSV IDS:*,numreads
    fout = open(output_file,'w')
    f=open(bedpe_file)
    for line in f:
        neat = line.strip().split()
        variantID = neat[6]
        # readnames on the left:
        reads1 = readnames_1.get(variantID,set())
        # readnames on the right:
        reads2 = readnames_2.get(variantID,set())
        # readnames occurring on both sides
        num_split = len(reads1.intersection(reads2)) 
        ##  first side of the translocation
        key1 = (variantID,1)
        left1 = left_dict.get(key1,0)
        span1 = span_dict.get(key1,0)
        right1 = right_dict.get(key1,0)
        ## second side of the translocation
        key2 = (variantID,2)
        left2 = left_dict.get(key2,0)
        span2 = span_dict.get(key2,0)
        right2 = right_dict.get(key2,0)
        original_line = "\t".join(line.strip().split()[0:11]) # excluding number of split reads

        # Check if inverted duplication (reads loop back onto themselves at this breakpoint):
        chrom1 = neat[0]
        pos1 = int(neat[1]) # beginning of the interval 
        other_pos1 = int(neat[2]) # end of the first breakpoint interval

        chrom2 = neat[3]
        pos2 = int(neat[4]) # beginning of the interval 
        other_pos2 = int(neat[5]) # end of the second breakpoint interval


        if chrom1==chrom2 and (abs(other_pos1-pos2) < slop or abs(other_pos2-pos1) < slop):
            pos1 = min([pos1,pos2])
            pos2 = pos1
            other_pos1 = max([other_pos1,other_pos2])
            other_pos2 = other_pos1
            original_line = "%s\t%d\t%d\t%s\t%d\t%d\t%s" % (chrom1,pos1,other_pos1,chrom2,pos2,other_pos2,"\t".join(neat[6:11]))
            num_split = int(neat[11])
            # print "Inverted duplication:",original_line
        elif chrom1==chrom2 and (abs(other_pos1-pos2) < 100000 or abs(other_pos2-pos1) < 100000): # if positions are within 100 kb they may share reads naturally
            num_split = int(neat[11])
        
        fout.write(original_line+"\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n" % (num_split,left1,span1,right1,left2,span2,right2))
    #   chrom1 start1 stop1   chrom2 start2 stop2  variantID score strand1 strand2 typeOfSV numreads left1 span1 right1 left2 span2 right2

    f.close()
    fout.close()

    nodes(args)

    
def nodes(args):
    bedpe_file = args.out + ".spansplit.bedpe"
    output_file = args.out + ".spansplit.nodes.bed"

    small_chromosomes_only = args.filter
    chromosome_lengths_file = args.genome_file
    slop = args.slop
    chromosome_lengths = {}

    f=open(chromosome_lengths_file)
    for line in f:
        fields = line.strip().split()
        chromosome_lengths[fields[0]] = int(fields[1])

    
    ############################# cut the chromosomes at each breakpoint to make nodes ################################
    f=open(bedpe_file)
    breakpoints = {}
    # dictionary with an entry for each chromosome, each is a list of positions where the breakpoints are
    for line in f:
        neat = line.strip().split()
        # variantID = neat[6]
        chrom1 = neat[0]
        pos1 = int(neat[1]) # beginning of the interval becomes the official breakpoint (for cutting nodes apart)
        chrom2 = neat[3]
        pos2 = int(neat[4]) # beginning of the interval becomes the official breakpoint (for cutting nodes apart)
        # filter out any breakpoints to or from the alternate chromosomes, patches, etc.
        if (len(chrom1) < 6 and len(chrom2) < 6) or not small_chromosomes_only:  
            # check if this is an in-place variant like an inverted duplication where the breakpoints are within SLOP of each other 
            # and should therefore be considered a single breakpoint
            # other_pos1 = int(neat[2]) # end of the first breakpoint interval
            # other_pos2 = int(neat[5]) # end of the second breakpoint interval
            # if chrom1==chrom2 and (abs(other_pos1-pos2) < slop or abs(other_pos2-pos1) < slop):
            #     # record only one breakpoint in the dictionary
            #     breakpoints[chrom1] = breakpoints.get(chrom1,[])+[pos1]
            # else:

            # record both breakpoints in the dictionary
            breakpoints[chrom1] = breakpoints.get(chrom1,[])+[pos1]
            breakpoints[chrom2] = breakpoints.get(chrom2,[])+[pos2]
    # output summary statistics on the breakpoints in each chromosome
    print "Chromosomes: %d" % len(breakpoints)
    print "Number of breakpoints in each chromosome:"
    for key in breakpoints:
        print "%s: %d" % (key,len(breakpoints[key]))
    f.close()

    distances = []
    fout = open(output_file,'w')
    for key in breakpoints:
        print "key:",key
        # sort positions within each chromosome
        array = np.array([0]+list(np.unique(breakpoints[key]))+[chromosome_lengths[key]]) 
        array.sort()
        # output file format:
        # chrom start stop name 
        for i in xrange(1,len(array)):
            distances.append(array[i]-array[i-1])
            fout.write("%s\t%d\t%d\t%d\n" % (key,array[i-1],array[i],i))
    fout.close()

def main():
    parser=argparse.ArgumentParser(description="Counts reads spanning and split at breakpoints. Useful for inspection and preparation for SplitThreader.")
    parser.add_argument("-bedpe",help="Variant calls in bedpe format (like from Sniffles)",dest="bedpe",required=True)
    parser.add_argument("-out",help="output prefix",dest="out",required=True)
    parser.add_argument("-reads",help="reads_near_breakpoints.bed",dest="reads",required=True)
    parser.add_argument("-genome-file",help="genome file: column 1 is chromosome name, column 2 is length of chromosome. Tab-separated.",dest="genome_file",required=True)
    parser.add_argument("-mq",help="mapping quality threshold (inclusive)",dest="mq",type=int,default=0)
    parser.add_argument("-bedcolumns",help="number of columns in read bedfile before intersect (default=6)",dest="offset",default=6)
    parser.add_argument("-filter",help="filter to primary chromosomes only (names shorter than 6 characters). Default=True",dest="filter",type=bool,default=True)
    parser.add_argument("-slop",help="interval size (bp) that the reads must span across to be counted (default=1000)",dest="slop",type=int,default=1000)
    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
    main()

