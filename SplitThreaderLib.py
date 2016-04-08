#! /usr/bin/env python

#######################################################################################################################################
################################################        SplitThreader.py         ######################################################
#######################################################################################################################################

# Node has 2 ports
    # Each port has an unrestricted number of edges to other ports
# Graph has nodes and edges


import random
import operator
import numpy as np
import time

class FileFormatError(Exception):
    pass

class Node(object):
    def __init__(self,name,attributes = {}):

        self.name = name
        self.ports = {}
        self.attributes = attributes
        self.weight = 0

        # make two ports always
        self.add_port("start")
        self.add_port("stop")

    @property
    def x(self):
        # import random
        return self.attributes.get("x",random.random())

    @property
    def y(self):
        # import random
        return self.attributes.get("y",random.random())

    def add_port(self,port_name):
        self.ports[port_name] = Port(port_name,self)

    # human-readable outputs:
    def __str__(self):
        return "Node: %s" % self.name
    def __repr__(self):
        return self.name


class Port(object):
    def __init__(self,name,node):
        self.name = name
        self.node = node
        # self.edges = [] # keep it a list instead of a dictionary so we can have multiple edges between the same ports
        self.edges = {}

    def add_edge(self,edge_instance):
        # self.edges.append(edge_instance)
        self.edges[edge_instance.glide(self)] = edge_instance

    # human-readable outputs:
    def __str__(self):
        return "%s:%s" % (self.node.name, self.name)
    def __repr__(self):
        return "%s:%s" % (self.node.name, self.name)

    def jump(self):
        ports = self.node.ports
        if self == ports["start"]:
            return ports["stop"]
        elif self == ports["stop"]:
            return ports["start"]
        else:
            assert("jump(): current port not in edge")


class Edge(object):
    def __init__(self,port1,port2):
        self.ports = (port1,port2)
        self.weight = 0
        self.spansplit = "split" # "span" or "split"
        self.variant_name = ""

    @property
    def n1(self): return self.ports[0].node.name
    @property
    def n2(self): return self.ports[1].node.name
    @property
    def p1(self): return self.ports[0].name
    @property
    def p2(self): return self.ports[1].name

    def __str__(self):
        names = map(str,self.ports)
        return "--".join(names)

    def __repr__(self):
        names = map(str,self.ports)
        return "--".join(names)

    def glide(self,current_port):
        if self.ports[0] == current_port:
            return self.ports[1]
        elif self.ports[1] == current_port:
            return self.ports[0]
        else:
            assert("glide(): current port not in edge")

def reverse(side):
    if side == "+":
        return "-"
    elif side == "-":
        return "+"
    elif side == "start":
        return "stop"
    elif side == "stop":
        return "start"
    else:
        raise Exception("reverse only takes +, -, start, stop")

class Graph(object):
    FileFormatError=FileFormatError

    def __init__(self):
        self.nodes = {} # key: node_name; value: Node object
        self.edges = [] # list of Edge objects
        self.node_lookup = {} # key: (chrom, pos, strand) tuple; #value: node_name
        self.annotation = {}
        self.annotation_by_node = {}

    def __str__(self):
        return ", ".join(self.nodes.keys())

    def __repr__(self):
        return "Graph: [" +", ".join(self.nodes.keys()) + "]"

    def print_nodes(self):
        for node in self.nodes:
            print node, self.nodes[node].attributes
    def print_edges(self):        
        for edge in self.edges:
            print edge, edge.weight, edge.spansplit

    def create_nodes_with_attributes(self,node_attributes):
        for node_name in node_attributes:
            if node_name in self.nodes.keys():
                print "Node already exists in graph, not overwriting"
            else:
                self.nodes[node_name] = Node(node_name,node_attributes[node_name])
                attr = node_attributes[node_name]
                if "chrom" in attr and "start" in attr and "stop" in attr:
                    self.node_lookup[(attr["chrom"],attr["start"],"-")] = node_name
                    self.node_lookup[(attr["chrom"],attr["stop"],"+")] = node_name

    def create_nodes(self,node_names):
        for node_name in node_names:
            if node_name in self.nodes.keys():
                print "Node already exists in graph, not overwriting"
            else:
                self.nodes[node_name] = Node(node_name)

    def create_edges(self,edge_list,weight_list = [],spansplit="split"):
        # edge_list contains node:port -- node:port pairs

        for i in xrange(len(edge_list)):
            edge = edge_list[i]
            side1 = edge[0]
            node1 = side1[0]
            port1 = side1[1]

            side2 = edge[1]
            node2 = side2[0]
            port2 = side2[1]
            # Create new edge between: self.nodes[node1].ports[port1]  ----  self.nodes[node2].ports[port2])
            e = Edge(self.nodes[node1].ports[port1],self.nodes[node2].ports[port2])
            if len(weight_list) == len(edge_list):
                e.weight = weight_list[i]
            e.spansplit = spansplit

            # Add edge to Graph
            self.edges.append(e)
            # Add edge to both ports
            self.nodes[node1].ports[port1].edges[self.nodes[node2].ports[port2]] = e
            self.nodes[node2].ports[port2].edges[self.nodes[node1].ports[port1]] = e


    def create_3_edges(self,key1,key2,weight_split,weight_span1,weight_span2,verbose=False,split_variant_name=""):
    
        port1 = ""
        if key1[2] == "+":
            port1 = "stop"
        elif key1[2] == "-":
            port1 = "start"
        else:
            raise FileFormatError("strands must be + or - in columns 9 and 10")

        port2 = ""
        if key2[2] == "+":
            port2 = "stop"
        elif key2[2] == "-":
            port2 = "start"
        else:
            raise FileFormatError("strands must be + or - in columns 9 and 10")


        # split
        node1 = self.node_lookup.get(key1, "NA")
        node2 = self.node_lookup.get(key2, "NA")
        if node1 == "NA" or node2 == "NA":
            print "Keys attempted"
            print key1
            print key2
            print "___________________________________"
            print "Dictionary:"
            print self.node_lookup
            raise FileFormatError("node in spansplit file doesn't exist in spansplit.nodes file")
        self.create_edge(node1,port1,node2,port2,weight_split,spansplit="split",split_variant_name=split_variant_name)
        if verbose==True:
            print "nodes:", node1, node2
            
        # span 1
        rev_key1 = (key1[0],key1[1],reverse(key1[2]))
        rev_node1 = self.node_lookup.get(rev_key1, "NA")
        rev_port1 = reverse(port1)
        self.create_edge(node1,port1,rev_node1,rev_port1,weight_span1,spansplit="span")

        # span 2 if not the same as span 1
        if key1 != key2:
            rev_key2 = (key2[0],key2[1],reverse(key2[2]))
            rev_node2 = self.node_lookup.get(rev_key2, "NA")
            rev_port2 = reverse(port2)
            
            self.create_edge(node2,port2,rev_node2,rev_port2,weight_span2,spansplit="span")
        elif verbose == True:
            print "Inverted duplication detected at ", key1


    def create_edge(self,node1,port1,node2,port2,weight,spansplit="split",split_variant_name=""):
        
        e = Edge(self.nodes[node1].ports[port1], self.nodes[node2].ports[port2])
        e.weight = weight
        e.spansplit=spansplit
        e.variant_name = split_variant_name
        self.edges.append(e)
        self.nodes[node1].ports[port1].edges[self.nodes[node2].ports[port2]] = e
        self.nodes[node2].ports[port2].edges[self.nodes[node1].ports[port1]] = e



    def add_portal(self,portal_name="Portal",all_ports=False,fill_to_flow = True):
        # Set up portals at the unconnected nodes
        portals = set()
        portal_edge_weights = {}
        for node_name in self.nodes:
            no_start = False
            no_stop = False
            if len(self.nodes[node_name].ports["start"].edges.values()) == 0:
                no_start = True
            if len(self.nodes[node_name].ports["stop"].edges.values()) == 0:
                no_stop = True
            if no_start and no_stop:
                pass # Don't connect a node to a portal if it has no edges
            elif all_ports == True:
                portals.add(self.nodes[node_name].ports["start"])
                portals.add(self.nodes[node_name].ports["stop"])
            elif fill_to_flow == True:
                sum_port1, sum_port2 = self.calculate_flow_on_node(node_name)
                node_weight = self.nodes[node_name].weight
                if sum_port1 < node_weight:
                    portals.add(self.nodes[node_name].ports["start"])
                    portal_edge_weights[self.nodes[node_name].ports["start"]] = node_weight - sum_port1
                if sum_port2 < node_weight:
                    portals.add(self.nodes[node_name].ports["stop"])
                    portal_edge_weights[self.nodes[node_name].ports["stop"]] = node_weight - sum_port2
            elif no_start:
                portals.add(self.nodes[node_name].ports["start"])
            elif no_stop:
                portals.add(self.nodes[node_name].ports["stop"])

        self.create_nodes_with_attributes({portal_name:{"chrom":"0","start":0,"stop":0,"x":0,"y":0}})
        for dead_end_port in portals:
            weight = portal_edge_weights.get(dead_end_port,float("inf"))
            self.create_edge(dead_end_port.node.name, dead_end_port.name, portal_name, "stop", weight,spansplit="portal")
            # self.create_edge(node1,port1,node2,port2,weight):
        
        # self.print_edges()

    def edges_from_path(self,path,split_edges_only = False):
        second_port = self.port_from_path_item(path[0])
        edge_list = []
        for item in path[1:]:
            first_port = self.opposite_port_from_path_item(item)
            this_edge = second_port.edges[first_port]
            if split_edges_only == True:
                if this_edge.spansplit == "split":
                    edge_list.append(this_edge)
            else:
                edge_list.append(this_edge)
            second_port = self.port_from_path_item(item)
        return edge_list

    def port_from_path_item(self,item):
        node,port = item.split(":")
        return self.nodes[node].ports[port]

    def opposite_port_from_path_item(self,item):
        node,port = item.split(":")
        return self.nodes[node].ports[reverse(port)]

    def min_weight(self,path, by_nodes = False):
        
        if not by_nodes: # use edge weights
            edge_list = self.edges_from_path(path)
            edge_dict = {}
            for edge in edge_list:
                if edge.weight == float("inf"):
                   pass
                else:
                    edge_dict[edge] = edge_dict.get(edge,0) + 1

            minweight = -1
            for edge in edge_dict:
                corrected_weight = edge.weight*1.0/edge_dict[edge]
                if minweight == -1 or corrected_weight < minweight:
                    minweight = corrected_weight
            if minweight == -1:
                print "no edges except to portal, using node weights instead"
                return self.min_weight(path,by_nodes=True)
            else:
                return minweight

        else: # use node weights
            node_names = self.node_names_from_path(path)

            node_dict = {}
            for node_name in node_names:
                if node_name != "Portal":
                    node_dict[node_name] = node_dict.get(node_name,0) + 1

            minweight = -1
            for node_name in node_dict:
                node = self.nodes[node_name]
                corrected_weight = node.weight*1.0/node_dict[node_name]
                if minweight == -1 or corrected_weight < minweight:
                    minweight = corrected_weight

            return minweight

    def find_nodename_by_position(self,point1): # point1 = ("chrom",position)
        chrom1 = point1[0]
        pos1 = point1[1]

        matching_nodes = []
        for node_name in self.nodes:
            if chrom1 == self.nodes[node_name].attributes["chrom"]:
                # if a point is on the breakpoint between nodes, assign it to the left node
                if pos1 > self.nodes[node_name].attributes["start"] and pos1 <= self.nodes[node_name].attributes["stop"]: 
                    matching_nodes.append(node_name)
        if len(matching_nodes) > 1:
            print point1
            raise FileFormatError("More than one node contains the point")
        if len(matching_nodes) == 0:
            print point1
            raise FileFormatError("No nodes match the point")
        return matching_nodes[0]

    def find_nodenames_by_gene(self,gene_name):
        error_message = "Gene not in annotation"
        annot = self.annotation.get(gene_name,error_message)
        if annot == error_message:
            return None
        node_name_start = self.find_nodename_by_position((annot["chrom"],annot["start"]))
        node_name_stop = self.find_nodename_by_position((annot["chrom"],annot["stop"]))
        if node_name_start != node_name_stop:
            return [node_name_start,node_name_stop]
        else:
            return [node_name_start]

    def calculate_distance(self,point1,point2,depth_limit=5):
        chrom1 = point1[0]
        chrom2 = point2[0]
        pos1 = point1[1]
        pos2 = point2[1]

        node_name1 = self.find_nodename_by_position(point1)
        node_name2 = self.find_nodename_by_position(point2)
        if node_name1 == node_name2:
            if pos2 > pos1: # forward direction, return stop port as path
                return [str(self.nodes[node_name1].ports["stop"])], abs(pos2-pos1) 
            if pos1 > pos2: # reverse direction, return start port as path
                return [str(self.nodes[node_name1].ports["start"])], abs(pos2-pos1) 
        else:
            node1 = self.nodes[node_name1]
            node2 = self.nodes[node_name2]

        
        all_connections = []
        all_distances = []
        
        for strand1 in ["+","-"]:
            for strand2 in ["+","-"]:
                paths,distances = self.find_paths_by_directions(node1,node2,pos1,pos2,strand1,strand2,depth_limit=depth_limit)        
                all_connections += paths
                all_distances += distances

        min_distance = -1
        shortest_path = []
        for i in xrange(len(all_distances)):
            if all_distances[i] < min_distance or min_distance == -1:
                min_distance = all_distances[i]
                shortest_path = all_connections[i]
        return shortest_path,min_distance

    def find_path_total_lengths(self,allpaths,count_only_nodes_in_set=[]):
        total_lengths = []
        for path in allpaths:
            total_lengths.append(self.find_total_length(path,count_only_nodes_in_set))
        return total_lengths

    def find_total_length(self,path,count_only_nodes_in_set=[]):
        distance = 0
        for item in path:
            intermediate_node = self.port_from_path_item(item).node
            if len(count_only_nodes_in_set)==0 or intermediate_node.name in count_only_nodes_in_set:
                distance += intermediate_node.attributes["stop"] - intermediate_node.attributes["start"]
        return distance

    def find_paths_by_directions(self,node1,node2,pos1,pos2,strand1,strand2,depth_limit=5,split_depth_limit=3):
        node1_start_port = "start" if strand1=="+" else "stop"
        node2_end_port = "stop" if strand2=="+" else "start"

        allpaths = self.breadth_first_search(node1.ports[node1_start_port], node2.ports[node2_end_port], depth_limit = depth_limit, stop_when_found = False,split_depth_limit=split_depth_limit)
        distances = []
        for path in allpaths:
            distance = 0
            for item in path[1:-1]:
                intermediate_node = self.port_from_path_item(item).node
                distance += intermediate_node.attributes["stop"] - intermediate_node.attributes["start"]

            distance += abs(node1.attributes[reverse(node1_start_port)]-pos1)
            distance += abs(pos2 - node2.attributes[reverse(node2_end_port)])
            distances.append(distance)
            # print node1.attributes["stop"]
            # total length of intermediate nodes + partial lengths of first and last nodes
        return allpaths,distances


    def read_annotation(self,annotation_file,name_field=4,strand_field=6,by_node_access = False):
        # It takes much longer to index all the genes by which nodes contain them, so by_node_access is turned off by default. It is useful to turn on when you want to see all the genes in a given path, region, or cycle
        f = open(annotation_file)
        counter = 0
        successes = 0
        # print "sample gene names used:"
        for line in f:
            fields = line.strip().split()
            chrom = fields[0]
            start = float(fields[1])
            stop = float(fields[2])
            gene_name = fields[name_field-1]
            strand = fields[strand_field-1]
            self.annotation[gene_name] = {"chrom":chrom,"start":start,"stop":stop,"strand":strand}
            if by_node_access == True:
                # get names of all nodes that are between these two points
                all_nodes_on_gene = self.nodes_within_genome_interval(chrom,start,stop)
                if len(all_nodes_on_gene) > 0:
                    successes += 1
                for node_name in all_nodes_on_gene:
                    self.annotation_by_node[node_name] = self.annotation_by_node.get(node_name,[]) + [gene_name]
            counter += 1
            # if counter < 5:
            #     print gene_name
        f.close()


        if by_node_access == True:
            print "Reading annotation: %.1f %% of genes placed successfully" % (successes*100./counter)
            
            if successes*100./counter < 20:
                print "Check whether the chromosome names in the annotation file matches the variants and nodes files."


    def node_names_from_path(self,path):
        node_names = []
        for item in path:
            node_names.append(self.port_from_path_item(item).node.name)
        return node_names

    def genes_on_path(self,path):
        genes = []
        for item in path:
            node_name = self.port_from_path_item(item).node.name
            genes = genes + self.annotation_by_node.get(node_name,[])
        return list(set(genes))

    def gene_fusion_distance(self,gene_name1,gene_name2,additional_info = None, depth_limit=20, verbose=False,split_depth_limit = 3):
        annot1 = self.annotation.get(gene_name1)
        annot2 = self.annotation.get(gene_name2)
        
        reports = []
        for gene_end1 in ["start","stop"]:
            pos1 = annot1[gene_end1]
            node_name1 = self.find_nodename_by_position((annot1["chrom"],pos1))
            node1 = self.nodes[node_name1]
            if verbose:
                print "1:",node1
                
            for gene_end2 in ["start","stop"]:
                pos2 = annot2[gene_end2]
                node_name2 = self.find_nodename_by_position((annot2["chrom"],pos2))
                node2 = self.nodes[node_name2]
                if verbose:
                    print "2:",node2
                for strand1 in ["+","-"]:
                    for strand2 in ["+","-"]:
                        paths,distances = self.find_paths_by_directions(node1,node2,pos1,pos2,strand1,strand2,depth_limit=depth_limit,split_depth_limit=split_depth_limit)  
                        if len(paths)>0:
                            direction_gene1 = "sense"
                            if (annot1["strand"] == "-" and gene_end1 == "stop") or (annot1["strand"]=="+" and gene_end1 == "start"):
                                pass
                            else:
                                direction_gene1 = "antisense"
                            direction_gene2 = "sense"
                            if (annot2["strand"] == "-" and gene_end2 == "start") or (annot2["strand"]=="+" and gene_end2 == "stop"):
                                pass
                            else:
                                direction_gene2 = "antisense"
                            # min_distance = distances[0]
                            # path = paths[0]
                            # if len(distances) > 1:
                            #     for i in xrange(len(distances)):
                            #         if distances[i] < min_distance:
                            #             min_distance = distances[i]
                            #             path = paths[i]
                            for i in xrange(len(paths)):
                                reports.append({"Gene1":gene_name1,"Gene2":gene_name2, "Gene1_direction":direction_gene1, "Gene2_direction":direction_gene2, "path":paths[i],"distance":distances[i],"info":additional_info})
        return reports

    def count_splits_in_path(self,path):
        num_splits = 0
        edges = self.edges_from_path(path)
        for edge in edges:
            if edge.spansplit == "split":
                num_splits += 1

        return num_splits

    def split_weights_in_path(self,path):
        split_weights = []
        edges = self.edges_from_path(path)
        for edge in edges:
            if edge.spansplit == "split":
                if edge.weight == int(edge.weight):
                    split_weights.append(int(edge.weight))
                else:
                    split_weights.append(edge.weight)

        
        return split_weights

    def gene_fusion_report(self,gene_name1,gene_name2,additional_info = None, depth_limit=15,verbose = False,output_file_all_candidates=None):
        if verbose:
            print gene_name1,"-",gene_name2

        reports = self.gene_fusion_distance(gene_name1,gene_name2,additional_info = additional_info, depth_limit=depth_limit,verbose=False)
            
        if len(reports) == 0:
            if verbose:
                print "No gene fusion detected"
        else:
            if verbose:
                print len(reports), "gene fusion(s) detected"

            scores = []
            for report in reports:

                # Contents of report:
                #{"Gene1":gene_name1,"Gene2":gene_name2, "Gene1_direction":direction_gene1, "Gene2_direction":direction_gene2, "path":paths[i],"distance":distances[i],"info":additional_info}

                all_split_weights = self.split_weights_in_path(report["path"])
                score = 0
                num_splits = self.count_splits_in_path(report["path"])

                report["number_of_splits"] = num_splits
                report["split_weights"] = all_split_weights
                report["variant_names"] = self.variants_from_path(report["path"])
                report["rna_split_read_count"] = int(report["info"][0])
                
                if num_splits==1 and report["distance"]<100000 and report["Gene1_direction"] == report["Gene2_direction"]:
                    score = 150
                elif num_splits==2 and report["distance"]<100000 and report["Gene1_direction"] == report["Gene2_direction"]: 
                    score = 120
                elif num_splits==1 and report["distance"]<1000000 and report["Gene1_direction"] == report["Gene2_direction"]:
                    score = 100
                elif num_splits==1 and report["distance"]<1000000: 
                    score = 70
                elif num_splits==2 and report["distance"]<1000000 and report["Gene1_direction"] == report["Gene2_direction"]:
                    score = 50
                elif num_splits==2 and report["distance"]<1000000:
                    score = 40
                else:
                    score = 20
                if num_splits == 0:
                    score = 0
                else:
                    score += np.min([np.min(all_split_weights)/10,20])
                
                scores.append(score)

            scores = np.array(scores)
            indices = np.argsort(scores)[::-1]

            to_return = None
            for index in indices:
                report = reports[index]
                if output_file_all_candidates != None:
                    output_file_all_candidates.write("%s,%s,%s,%s,%s,%s,%d,%d,%s,%s,%d,%.2f\n" % (report["Gene1"],report["Gene1_direction"], self.annotation.get(report["Gene1"])["chrom"], report["Gene2"], report["Gene2_direction"], self.annotation.get(report["Gene2"])["chrom"], report["rna_split_read_count"], report["number_of_splits"], "|".join(map(str,report["split_weights"])),  "|".join(report["variant_names"]), report["distance"], scores[index]   )) 
                    
                if scores[index] == max(scores):
                    # Flip the genes around if both are antisense
                    if report["Gene1_direction"]=="antisense" and report["Gene2_direction"]=="antisense":
                        report["Gene1_direction"]="sense"
                        report["Gene2_direction"]="sense"
                        tmp = report["Gene1"]
                        report["Gene1"] = report["Gene2"]
                        report["Gene2"] = tmp

                    to_return = report
            if verbose:
                print '__________________________'
            if max(scores) <= 40:
                return None

            return to_return

    def find_intact_length(self,path,count_only_nodes_in_set=[]):
        longest_uninterrupted_length_so_far = 0
        current_uninterrupted_length = 0
        current_chromosome = ""
        position_where_we_left_off = 0
        for item in path:
            node,port = item.split(":")
            
            # Attributes: # {"chrom":fields[0],"start":int(fields[1]),"stop":int(fields[2]),"x":int(fields[1]),"y":y}
            seq_length = self.nodes[node].attributes["stop"]-self.nodes[node].attributes["start"]
            this_chromosome = self.nodes[node].attributes["chrom"]
            this_position = self.nodes[node].attributes[reverse(port)] # port refers to after the jump, so we reverse that to get the entry point into this node
            if (len(count_only_nodes_in_set) == 0 or node in count_only_nodes_in_set) and (this_chromosome == current_chromosome and this_position == position_where_we_left_off): 
                # count the sequence length only if this node is within the special set (for regional parsimony)
                current_uninterrupted_length += seq_length
            else:
                # Save if this path is the best so far
                if current_uninterrupted_length > longest_uninterrupted_length_so_far:
                    longest_uninterrupted_length_so_far = current_uninterrupted_length
                    # longest_uninterrupted_path_so_far = path
                # Reset chromosome and length
                current_uninterrupted_length = seq_length
                current_chromosome = this_chromosome
            position_where_we_left_off = self.nodes[node].attributes[port] # port refers to after the jump, so that reflects the exit port out of this node
        return longest_uninterrupted_length_so_far

    def find_path_intact_lengths(self,allpaths,count_only_nodes_in_set=[]):
        intact_lengths = []
        for path in allpaths:
            intact_lengths.append(self.find_intact_length(path,count_only_nodes_in_set))
        return intact_lengths

    def subtract(self,path,weight, by_nodes = False):
        # print "Subtracting", weight
        if by_nodes == True:
            nodes = self.node_names_from_path(path)
            for node_name in nodes:
                node = self.nodes[node_name]
                node.weight = node.weight - weight
        else:
            edges = self.edges_from_path(path)
            for edge in edges:
                edge.weight = edge.weight - weight
                # if edge.weight == 0:
                    # print "zeroed out", edge

    def breadth_first_search(self, current_port, destination_port, depth_limit = -1, stop_when_found = False, split_depth_limit=None):
        
        ###############  Basic steps:  #####################################
        # jump
        # make queue with current_port in it
        # while something is in the queue:
            # find all edges of current_port
                # glide
                # jump
                # if this new port matches
                    # return
                # else 
                    # add it to the queue
            # let the while loop repeat to keep exploring the queue
        #############################################################

        allpaths = []

        current_port = current_port.jump()
        # if match, then return:
        if str(current_port) == str(destination_port):
            if stop_when_found:
                return [[str(destination_port)]]
            else:
                allpaths.append([str(destination_port)])

        queue = [(current_port,[current_port])]
        while queue:
            (port, path) = queue.pop(0)
            if split_depth_limit != None and len(path)>0 and self.count_splits_in_path(map(str,path)) > split_depth_limit:
                    continue
                # return allpaths
            if depth_limit != -1 and len(path) > depth_limit:
                return allpaths
            edges = port.edges
            for edge in edges.values():
                # ignore if already in path
                # glide
                current_port = edge.glide(port)
                # jump
                current_port = current_port.jump()
                # if match, then return
                if str(current_port) == str(destination_port):
                    if stop_when_found:
                        return [map(str,path) + [str(destination_port)]]
                    else:
                        allpaths.append(map(str,path) + [str(destination_port)])
                # else append to queue
                else:
                    queue.append((current_port, path+[current_port]))
        return allpaths



    def depth_first_search_recurse(self,current_port,destination_port,allpaths,depth_limit,cycle_limit,path_so_far=[],stop_when_found = False,depth=0):
        # saving the ports after jumping, so if the path contains A:start, it means you went through A in the reverse direction. A:stop means forward direction. 

        ############# Basic steps: ################
        # jump (read sequence)
        # if current_port == destination_port:
        #     return allpaths + [path_so_far]
        # find edges of port
        # for edge in edges:
            # glide
            # recurse
        ###########################################
        
        jumped_port = current_port.jump()
        
        saveport = str(jumped_port)
        if str(jumped_port) == str(destination_port):
            # print "MATCH"
            allpaths.append(path_so_far + [saveport]) # new
        else:
            edges = jumped_port.edges.values()
            for edge in edges:
                glide_port = edge.glide(jumped_port)
                if stop_when_found and len(allpaths)>0 or depth > depth_limit or path_so_far.count(str(jumped_port)) > cycle_limit:
                    return
                else: # keep recursing
                    self.depth_first_search_recurse(current_port=glide_port, destination_port=destination_port, allpaths=allpaths, path_so_far=path_so_far+[saveport],stop_when_found=stop_when_found,depth_limit=depth_limit,cycle_limit=cycle_limit,depth=depth+1) # new

    def depth_first_search(self,current_port,destination_port,stop_when_found = False,depth_limit=1000,cycle_limit=2):
        allpaths = []
        self.depth_first_search_recurse(current_port=current_port,destination_port=destination_port,allpaths=allpaths,stop_when_found=stop_when_found,depth_limit=depth_limit,cycle_limit=cycle_limit)
        return allpaths

    def group_redundant_gene_fusions(self,paths_dict_by_fusion_name,split_edges_only = False):
        paths_used = []
        edge_lists_used = []
        names_used = []

        groups = {}

        for fusion_name in paths_dict_by_fusion_name:
            path = paths_dict_by_fusion_name[fusion_name]
            edge_list = set(self.edges_from_path(path,split_edges_only=split_edges_only))
            if paths_used == []:
                paths_used.append(path)
                edge_lists_used.append(edge_list)
                names_used.append(fusion_name)
                groups[fusion_name] = [fusion_name]
            else:
                redundant = False
                for j in xrange(len(edge_lists_used)):
                    previous_edge_list = edge_lists_used[j]
                    if edge_list == previous_edge_list:
                        redundant = True
                        groups[names_used[j]].append(fusion_name)
                        break
                if redundant == False:
                    paths_used.append(path)
                    edge_lists_used.append(edge_list)
                    names_used.append(fusion_name)
                    groups[fusion_name] = [fusion_name]

        return groups.values()

    def collapse_redundant_paths(self,paths,split_edges_only = False,names = None):
        paths_used = []
        edge_lists_used = []
        names_used = []
        if names != None and len(names)!=len(paths):
            print "To label the paths with names, the length of names (list) must be the same as the length of paths (list)"
            return None

        i = 0
        for path in paths:
            edge_list = set(self.edges_from_path(path,split_edges_only=split_edges_only))
            if paths_used == []:
                paths_used.append(path)
                edge_lists_used.append(edge_list)
                if names != None:
                    names_used.append(names[i])
                    print names[i]
            else:
                redundant = False
                for j in xrange(len(edge_lists_used)):
                    previous_edge_list = edge_lists_used[j]
                    if edge_list == previous_edge_list:
                        redundant = True
                        if names != None:
                            # print names_used[j], "<----", names[i], "\t(same path, possibly overlapping genes)"
                            print "\t%s\t(same fusion path as %s)" % (names[i], names_used[j])
                        break
                if redundant == False:
                    paths_used.append(path)
                    edge_lists_used.append(edge_list)
                    if names != None:
                        names_used.append(names[i])
                        print names[i]
            i += 1

        return paths_used

    def depth_first_search_for_group(self,current_port,destination_ports,allpaths,depth_limit,cycle_limit,path_so_far=[],stop_when_found = False,depth=0):
        jumped_port = current_port.jump()   
        
        saveport = str(jumped_port)
        if depth > 0 and (str(jumped_port.node.name) in map(str,destination_ports)):
            # print "MATCH"
            allpaths.append(path_so_far + [saveport]) # new
        else:
            edges = jumped_port.edges.values()
            for edge in edges:
                glide_port = edge.glide(jumped_port)
                if stop_when_found and len(allpaths)>0 or depth > depth_limit or path_so_far.count(str(jumped_port)) > cycle_limit:
                    return
                else: # keep recursing
                    self.depth_first_search_for_group(current_port=glide_port, destination_ports=destination_ports, allpaths=allpaths, path_so_far=path_so_far+[saveport],stop_when_found=stop_when_found,depth_limit=depth_limit,cycle_limit=cycle_limit,depth=depth+1)


    def special_subgraph_from_genome_interval(self,chromosome,start,end,degree,verbose=True):
        primary_nodes = self.nodes_within_genome_interval(chromosome,start,end)
        connection_path_nodes = set()

        # Search recursively for the primary nodes and if found, add all nodes in the path to the final list of nodes for this subgraph
        for primary_node in primary_nodes:
            node = self.nodes[primary_node]
            for port in node.ports:
                current_port = node.ports[port]
                allpaths = []
                self.depth_first_search_for_group(current_port=current_port,destination_ports=primary_nodes,allpaths=allpaths,stop_when_found=False,depth_limit=degree,cycle_limit=2)
                path_nodes = set()
                for path in allpaths:
                    for path_item in path:
                        path_nodes.add(path_item.split(":")[0])
                connection_path_nodes = connection_path_nodes.union(path_nodes)

        s = self.enclosed_subgraph_from_nodes(primary_nodes = primary_nodes, intermediate_nodes=list(connection_path_nodes))
        # print "SPECIAL NODES:", s.nodes
        return s

    def enclosed_subgraph_from_nodes(self,primary_nodes,intermediate_nodes):
        # primary_nodes and intermediate_nodes are lists of strings
        
        # print "enclosed_subgraph_from_nodes()"

        s = Graph() # s as in subgraph

        region_nodes = primary_nodes + intermediate_nodes

        # print "Region nodes:"
        # print region_nodes

        node_attributes = {}
        edges_to_add = {}

        portals = set()

        for node_name in region_nodes:
            # Copy the node itself
            node_attributes[node_name] = self.nodes[node_name].attributes
            # Find all the first-degree nodes and their edges
            for port_name in self.nodes[node_name].ports:
                edges = self.nodes[node_name].ports[port_name].edges
                if len(edges) == 0:
                    # print "NO EDGES:", port_name
                    portals.add(self.nodes[node_name].ports[port_name])
                for other_port in edges:
                    edge = edges[other_port]

                    # If other side is also in region or this node gets an exit
                    within_region = str(other_port.node.name) in map(str,region_nodes)
                    # if within_region: 
                    #     print "within_region", edge

                    split_edge_outside_region = (node_name in primary_nodes and edge.spansplit=="split")
                    # if split_edge_outside_region:
                    #     print "split_edge_outside_region", edge

                    if within_region or split_edge_outside_region:
                        # print node_name, "connecting to", other_port.node
                        # Grab the node on the other side
                        node_attributes[other_port.node.name] = other_port.node.attributes
                        # Create a dictionary of edges to avoid using the same multiple times
                        edges_to_add[str(edge)] = edge

                        if split_edge_outside_region:
                            # print other_port.jump()
                            portals.add(other_port.jump())
                    else:
                        # print node_name, "NOT connecting to", other_port.node
                        # print str(other_port.node.name), "NOT IN", map(str,region_nodes)
                        portals.add(self.nodes[node_name].ports[port_name])

        
        # print "edges_to_add:", edges_to_add
        # print "all_nodes:", node_attributes.keys()
        # Add nodes to the graph
        s.create_nodes_with_attributes(node_attributes)

        for node_name in node_attributes:
            s.nodes[node_name].weight = self.nodes[node_name].weight

        edge_counter = 0
        # Add edges to the graph
        for edge_name in edges_to_add:
            edge = edges_to_add[edge_name]
            node1 = edge.ports[0].node.name
            port1 = edge.ports[0].name
            node2 = edge.ports[1].node.name
            port2 = edge.ports[1].name
            weight = edge.weight
            spansplit = edge.spansplit
            s.create_edge(node1,port1,node2,port2,weight=weight,spansplit=spansplit)
            edge_counter += 1
        # print "total edges:", edge_counter


        # Add portal 
        s.create_nodes_with_attributes({"Portal":{"chrom":"0","start":0,"stop":0,"x":0,"y":0}})
        for dead_end_port in portals:
            weight = float("inf") #portal_edge_weights.get(dead_end_port,float("inf"))
            s.create_edge(dead_end_port.node.name, dead_end_port.name, "Portal", "stop", weight,spansplit="portal")
            # self.create_edge(node1,port1,node2,port2,weight):



        return s

    def subgraph_from_genome_interval(self,chromosome,start,end,degree,verbose=True):
        nodes = self.nodes_within_genome_interval(chromosome,start,end)

        if verbose: print "subgraph nodes:", nodes
        
        s = self.subgraph_from_nodes(nodes,degree_given=degree)
        return s

    def nodes_within_genome_interval(self,chromosome,start,end):
        matching_nodes = []
        for node_name in self.nodes:
            if chromosome == self.nodes[node_name].attributes["chrom"]:
                # If at least one end of the node is within the interval (1 and 2) or the node encompasses the interval entirely (3)
                if (self.nodes[node_name].attributes["start"] > start and self.nodes[node_name].attributes["start"] < end) or \
                (self.nodes[node_name].attributes["stop"] > start and self.nodes[node_name].attributes["stop"] < end) or \
                (self.nodes[node_name].attributes["start"] < start and self.nodes[node_name].attributes["stop"] > end): 
                    matching_nodes.append(node_name)
        
        return matching_nodes

    def node_at_genome_position(self,chromosome,pos):
        matching_nodes = []
        for node_name in self.nodes:
            if chromosome == self.nodes[node_name].attributes["chrom"]:
                # If at least one end of the node is within the interval (1 and 2) or the node encompasses the interval entirely (3)
                if self.nodes[node_name].attributes["start"] <= pos and self.nodes[node_name].attributes["stop"] > pos: 
                    return node_name
        
        # if we don't find anything:
        return None

    def first_degree_nodes(self,node_names):
        first_degree_nodes = node_names + []
        for node_name in node_names:
            node = self.nodes[node_name]
            for port_name in self.nodes[node_name].ports:
                edges = self.nodes[node_name].ports[port_name].edges
                for other_port in edges:
                    first_degree_nodes.append(other_port.node.name)

        return list(set(first_degree_nodes))


    def subgraph_from_nodes(self,node_names_given,degree_given=1): 
        print "subgraph_from_nodes()"
        # Adds the nodes given along with all their first-degree nodes and the edges in between
        node_names = node_names_given + []
        degree = degree_given - 1 # the process later goes one degree out already, so only add degrees here if above 1

        for i in xrange(degree):
            node_names = self.first_degree_nodes(node_names)

        s = Graph() # s as in subgraph

        node_attributes = {}
        edges_to_add = {}
        for node_name in node_names:
            # Copy the node itself
            node_attributes[node_name] = self.nodes[node_name].attributes
            # Find all the first-degree nodes and their edges
            for port_name in self.nodes[node_name].ports:
                edges = self.nodes[node_name].ports[port_name].edges
                for other_port in edges:
                    edge = edges[other_port]
                    # Grab the node on the other side
                    node_attributes[other_port.node.name] = other_port.node.attributes
                    # Create a dictionary of edges to avoid using the same multiple times
                    edges_to_add[str(edge)] = edge

        
        print "all_nodes:", node_attributes.keys()
        # Add nodes to the graph
        s.create_nodes_with_attributes(node_attributes)

        edge_counter = 0
        # Add edges to the graph
        for edge_name in edges_to_add:
            edge = edges_to_add[edge_name]
            node1 = edge.ports[0].node.name
            port1 = edge.ports[0].name
            node2 = edge.ports[1].node.name
            port2 = edge.ports[1].name
            weight = edge.weight
            spansplit = edge.spansplit
            s.create_edge(node1,port1,node2,port2,weight=weight,spansplit=spansplit)
            edge_counter += 1
        print "total edges:", edge_counter

        return s

    def parsimony(self,use_breadth_first_search=False,portal_name="Portal",verbose=False,depth_limit=20,cycle_limit=0,min_weight_required = 10,chop_end_nodes=0,count_only_nodes_in_set=[]):
        # import time
        # self.add_portal(fill_to_flow = True ) # ????????????????????????????
        # recording = []

        before = time.time()
        allpaths = self.depth_first_search(self.nodes[portal_name].ports["start"],self.nodes[portal_name].ports["start"],cycle_limit=cycle_limit,depth_limit=depth_limit)
        if verbose:
            print "DFS:  %.2f seconds" % (time.time()-before)
            print "Number of paths", len(allpaths)
        
        # Chop off the end nodes on each path for the subgraph/local region case as these are inflated and can differentiate paths that are actually the same within the specific region
        if chop_end_nodes != 0:
            newpaths = []
            for path in allpaths:
                if len(path) > chop_end_nodes*2:
                    newpaths.append(path[chop_end_nodes:-chop_end_nodes])
            allpaths = newpaths

        allpaths = self.collapse_redundant_paths(allpaths)

        if verbose:
            print "Number of unique paths", len(allpaths)

        intact_lengths = self.find_path_intact_lengths(allpaths,count_only_nodes_in_set)
        total_lengths = self.find_path_total_lengths(allpaths,count_only_nodes_in_set)
        sorting_parameters = []

        # Sorting by 1) Longest intact path within region, 2) Longest total path within region, 3) Fewest splits
        for i in xrange(len(allpaths)):
            sorting_parameters.append((i,intact_lengths[i],total_lengths[i],-1*self.count_splits_in_path(allpaths[i])))
        
        ordering = sorted(sorting_parameters,key=operator.itemgetter(1,2,3))[::-1]

        recordings = []
        for item in ordering:
            index = item[0]
            path = allpaths[index]

            weight = self.min_weight(path,by_nodes=False)
            if weight < min_weight_required:
                continue
            else:
                recordings.append([path,intact_lengths[index],weight])
                self.subtract(path=path,weight=weight,by_nodes = False)

        return recordings

    def franken_paths(self,paths,output_filename,header = False):
        # print "Franken-path:"
        # print path
        path_counter = 1
        f=open(output_filename,"w")
        if header == True:
            f.write("chromosome\tstart\tend\tpath_ID\ttotal_length_of_path\tstrand\n")
        else:
            print "Header for bed file:\nchromosome\tstart\tend\tpath_ID\ttotal_length_of_path\tstrand\n"
        for path in paths:
            total_length = self.find_total_length(path)
            # intact_length = self.find_intact_length(path)
            for item in path:
                node_name,port = item.split(":")
                node = self.nodes[node_name]
                strand = "+"
                if port == "start":
                    strand = "-"
                f.write("%s\t%d\t%d\tpath_%03d\t%d\t%s\n" % (node.attributes["chrom"],node.attributes["start"], node.attributes["stop"],path_counter, total_length, strand))
            path_counter += 1      
        f.close()

    def local_evolution(self,chromosome,start,end,output_prefix,degree=3,min_weight_required=10,depth_limit=30):
    
        s = self.special_subgraph_from_genome_interval(chromosome,start,end,degree=degree)   

        print s

        s.count_span_vs_split_edges()
        
        print s.nodes
        print s.print_edges()

        print "Created subgraph of this region with %d nodes and %d edges\n" % (len(s.nodes),len(s.edges))
        # s.to_json(output_prefix+".subgraph.json")
        region_nodes = set(s.nodes_within_genome_interval(chromosome,start,end))

        print region_nodes

        reports = s.parsimony(depth_limit = depth_limit,chop_end_nodes=0,count_only_nodes_in_set=region_nodes,min_weight_required=min_weight_required) 

        filtered_reports = []
        output = ""
        for report in reports[::-1]:
            path = report[0]
            path_node_names = s.node_names_from_path(path)
            if len(set(path_node_names).intersection(region_nodes)):
                print report
                filtered_reports.append(report)
                for node in region_nodes:
                    if node in path_node_names:
                        output += "======"
                    else:
                        output += "      "
                output += "\n"
        print output
        filtered_reports = filtered_reports[::-1] # flip back around, to undo the ordering that is best for drawing, and return to the order that the paths were found in: longest intact sequence first

        print "Found a parsimonious set of", len(filtered_reports), "paths."
        
        s.boxes_from_parsimony(filtered_reports,output_filename=output_prefix+".boxes.csv")

    def boxes_from_parsimony(self,recordings,output_filename):
        path_counter = 1

        total_span_counts_by_node = {}
        f = open(output_filename,'w')
        f.write("chromosome,start,end,y_start,height,path_ID\n")
        for record in recordings:
            path = record[0][1:-1] # Cut off Portal nodes on either side of the path
            minweight = record[2]
            longest_intact_length = record[1]
            for item in path:
                node_name,port = item.split(":")
                node = self.nodes[node_name]
                y_start = total_span_counts_by_node.get(node_name,0) # current height of coverage for this node
                f.write("%s,%d,%d,%d,%d,path_%d\n" % (node.attributes["chrom"],node.attributes["start"], node.attributes["stop"],y_start,minweight,path_counter))
                total_span_counts_by_node[node_name] = total_span_counts_by_node.get(node_name,0) + minweight
            path_counter += 1  
        f.close()

    def check_overlapping_genes(self,gene1,gene2):
        
        if gene1 == gene2:
            return True
        annot1 = self.annotation[gene1]
        annot2 = self.annotation[gene2]

        if annot1["chrom"] == annot2["chrom"] and not (annot1["stop"] < annot2["start"] or annot2["stop"] < annot1["start"]):
            return True
        else:
            return False

    def group_redundant_gene_fusions_by_overlapping_genes(self,list_of_fusion_name_tuples):

        print "all fusions:", list_of_fusion_name_tuples
        
        # loop through all fusions in list_of_fusion_name_tuples
        # for each, check if it matches an existing group

        # if it matches multiple groups, combine those groups
        # add the new fusion to the matching group or combined group
        # if no matching group, make it a new group

        groups = []
        for fusion1 in list_of_fusion_name_tuples:
            # print "groups before:", groups
            matching_groups = []
            for index in xrange(len(groups)):
                group = groups[index]
                for fusion2 in group:
                    if (self.check_overlapping_genes(fusion1[0],fusion2[0]) and self.check_overlapping_genes(fusion1[1],fusion2[1])) or (self.check_overlapping_genes(fusion1[0],fusion2[1]) and self.check_overlapping_genes(fusion1[1],fusion2[0])):
                        matching_groups.append(index)
                        break

            if len(matching_groups) > 0:
                if len(matching_groups) > 1:
                    for index in matching_groups[1:]:
                        groups[matching_groups[0]] = groups[matching_groups[0]] + groups[index]
                        groups[index] = []
                groups[matching_groups[0]].append(fusion1)
            else:
                groups.append([fusion1])

        final_groups = []
        for group in groups:
            if len(group) > 0:
                final_groups.append(group)

        return final_groups
        
        # Return a list of lists with the gene fusions

    def variants_from_path(self,path):
        variant_names = []

        edges = self.edges_from_path(path,split_edges_only = True)

        for edge in edges:
            variant_names.append(edge.variant_name)
        return variant_names


    def search_for_fusions(self,gene_pair_list_file,output_prefix):
        f=open(gene_pair_list_file)
        
        gene_fusion_reports = []
        gene_names = []

        filename_all_candidates = output_prefix + ".all_candidates.csv"
        output_file_all_candidates = open(filename_all_candidates,"w")
        output_file_all_candidates.write("gene1,strand1,chrom1,gene2,strand2,chrom2,RNA_split_read_count,variant_count,DNA_split_reads_at_variants,variant_names,transcript_length,score\n")
            
        for line in f:
            fields = line.strip().split()
            gene_names.append((fields[0],fields[1]))
            gene_fusion_reports.append(self.gene_fusion_report(fields[0],fields[1],fields[2:],output_file_all_candidates=output_file_all_candidates))
        output_file_all_candidates.close()


        total_putative_fusions = len(gene_fusion_reports)
        total_supported_fusions = 0
        
        all_supported_fusion_paths = []
        all_supported_fusion_names = []
        paths_dict_by_fusion_name = {}
        reports_dict_by_fusion_name = {}

        for i in xrange(len(gene_fusion_reports)):
            report = gene_fusion_reports[i]
            if report == None:
                pass
            else:
                total_supported_fusions += 1
                key = (report["Gene1"],report["Gene2"])
                paths_dict_by_fusion_name[key] = report["path"]
                reports_dict_by_fusion_name[key] = report
                
        f.close()
        print "Found evidence of %d gene fusions out of the list of %d given:" % (total_supported_fusions,total_putative_fusions)

        list_of_fusion_name_tuples = paths_dict_by_fusion_name.keys()
        groups = self.group_redundant_gene_fusions_by_overlapping_genes(list_of_fusion_name_tuples)

        final_reports = []
        num_split_RNA_reads = []
        redundant_fusions = {}
        for group in groups:
            gene_lengths = {}
            RNA_counts = {}
            for fusion_name in group:
                annot1 = self.annotation[fusion_name[0]]
                # print annot1
                annot2 = self.annotation[fusion_name[1]]
                # print annot2
                length1 = annot1["stop"]-annot1["start"]
                length2 = annot2["stop"]-annot2["start"]
                gene_lengths[fusion_name] = length1 + length2

                RNA_counts[fusion_name] = int(reports_dict_by_fusion_name[fusion_name]["info"][0])
           
            max_fusion_name = ""
            max_RNA_count = 0
            for fusion_name in gene_lengths:
                if RNA_counts[fusion_name] > max_RNA_count or (RNA_counts[fusion_name] == max_RNA_count and gene_lengths[fusion_name] > gene_lengths[max_fusion_name]):
                    max_fusion_name = fusion_name
                    max_RNA_count = RNA_counts[fusion_name]
            redundant_fusions[max_fusion_name] = set(group)-set([max_fusion_name])

            final_reports.append(reports_dict_by_fusion_name[max_fusion_name])
            num_split_RNA_reads.append(int(reports_dict_by_fusion_name[max_fusion_name]["info"][0]))

        ###################################################################################################
        ##################     Table Output  and .csv Output for Visualizer    ############################
        ###################################################################################################

        f_output_fusion_report_csv = open(output_prefix + ".fusion_report.csv","w")
        f_output_fusion_report_csv.write("gene1,strand1,chrom1,gene2,strand2,chrom2,RNA_split_read_count,variant_count,DNA_split_reads_at_variants,variant_names,overlapping_fusions,transcript_length\n")

        # import numpy as np
        sorted_indices = np.argsort(num_split_RNA_reads)[::-1]
        for index in sorted_indices:
            report = final_reports[index]
            alternate_names = []
            for item in redundant_fusions[(report["Gene1"], report["Gene2"])]:
                alternate_names.append(item[0] + "=" + item[1])

            gene1_direction = report["Gene1_direction"]
            if report["Gene1_direction"] == "antisense":
                gene1_direction = "a"
            elif report["Gene1_direction"] == "sense":
                gene1_direction = "s"
            else:
                "ERROR: report['Gene1_direction'] =", report["Gene1_direction"]
            
            gene2_direction = report["Gene2_direction"]
            if report["Gene2_direction"] == "antisense":
                gene2_direction = "a"
            elif report["Gene2_direction"] == "sense":
                gene2_direction = "s"
            else:
                "ERROR: report['Gene2_direction'] =", report["Gene2_direction"]

            variant_names = self.variants_from_path(report["path"])
            f_output_fusion_report_csv.write("%s,%s,%s,%s,%s,%s,%d,%d,%s,%s,%s,%d\n" % (report["Gene1"],gene1_direction, self.annotation.get(report["Gene1"])["chrom"], report["Gene2"], gene2_direction, self.annotation.get(report["Gene2"])["chrom"], report["rna_split_read_count"], report["number_of_splits"], "|".join(map(str,report["split_weights"])),  "|".join(report["variant_names"]), ",".join(alternate_names),report["distance"]   )) 

        print "s = sense, a = anti-sense"
        f_output_fusion_report_csv.close()


    def breakpoints_from_sniffles(self,sniffles_filename):
        # Read file once to note all the breakpoints
        breakpoints_by_chromosome = {}
        f = open(sniffles_filename)
        for line in f:
            fields = line.strip().split()
            chrom1 = fields[0]
            pos1 = (int(fields[1]) + int(fields[2]))/2
            breakpoints_by_chromosome[chrom1] = breakpoints_by_chromosome.get(chrom1,[]) + [pos1]
            chrom2 = fields[3]
            pos2 = (int(fields[4]) + int(fields[5]))/2
            breakpoints_by_chromosome[chrom2] = breakpoints_by_chromosome.get(chrom2,[]) + [pos2]

        f.close()
        return breakpoints_by_chromosome


    def create_edges_from_sniffles(self,sniffles_filename):

        f=open(sniffles_filename)
        for line in f:
            fields = line.strip().split()
            chrom1 = fields[0]
            chrom2 = fields[3]
            pos1 = (int(fields[1]) + int(fields[2]))/2
            pos2 = (int(fields[4]) + int(fields[5]))/2
            strand1 = fields[8]
            strand2 = fields[9]
            weight_split = float(fields[11])
            weight_span1 = float(fields[12])
            weight_span2 = float(fields[13])

            variant_name = fields[6]

            key1 = (chrom1,pos1,strand1)
            key2 = (chrom2,pos2,strand2)

            self.create_3_edges(key1, key2, weight_split, weight_span1, weight_span2,split_variant_name = variant_name)

        f.close()

    def read_sniffles(self,sniffles_filename,genome_file):
        
        ###############     Create Nodes    #################
        # Find breakpoints
        breakpoints_by_chromosome = self.breakpoints_from_sniffles(sniffles_filename)
        
        # Split genome up into nodes using these breakpoints
        nodes = self.nodes_from_breakpoints(breakpoints_by_chromosome,genome_file)

        self.create_nodes_with_attributes(nodes)

        ###############     Create Edges    #################
        # Read file again to connect the nodes with edges
        self.create_edges_from_sniffles(sniffles_filename)


    def count_span_vs_split_edges(self):
        num_split = 0
        num_span = 0
        num_CNV = 0
        for edge in self.edges:
            if edge.spansplit=="split":
                num_split += 1
            elif edge.spansplit=="span":
                num_span += 1
            elif edge.spansplit=="CNV":
                num_CNV += 1
            else:
                print edge.spansplit
        print "split:",num_split
        print "span:",num_span
        if num_CNV > 0:
            print "CNVs:",num_CNV


    def create_span_edges_at_CNVs(self,breakpoints_from_copy_number):

        for chrom in breakpoints_from_copy_number:
            pos_list = breakpoints_from_copy_number[chrom]
            for pos in pos_list:
                key1 = (chrom,pos,"+")
                node1 = self.node_lookup.get(key1, "NA")
                port1 = "stop"
                if node1 == "NA":
                    print "Edge cannot be created because no node has this start or end position as a port:"
                    print key1
                    raise FileFormatError("trying to create edge at position that is not a node start or end port")
                rev_key1 = (key1[0],key1[1],reverse(key1[2]))
                rev_node1 = self.node_lookup.get(rev_key1, "NA")
                rev_port1 = reverse(port1)

                weight_span = min([self.nodes[node1].weight, self.nodes[rev_node1].weight])
                self.create_edge(node1,port1,rev_node1,rev_port1,weight_span,spansplit="CNV")
                
    def read_copy_number_profile_and_sniffles(self, sniffles_filename, genome_file, coverage_file, min_distance_to_variants_before_cutting_at_CNV = 100000, cov_diff_threshold_to_split=None):


        if cov_diff_threshold_to_split == None:
            cov_diff_threshold_to_split = self.average_coverage(coverage_file)/2


        ###############     Create Nodes    #################
        # Find breakpoints
        breakpoints_from_sniffles = self.breakpoints_from_sniffles(sniffles_filename)
        breakpoints_from_copy_number = self.segmented_coverage_to_CNV_calls(coverage_file,cov_diff_threshold_to_split=cov_diff_threshold_to_split)

        ####################################################
        ############   Magical merging step   ##############
        
        all_chromosomes = list(set(breakpoints_from_sniffles.keys() + breakpoints_from_copy_number.keys()))

        final_breakpoints_from_copy_number = {}

        NUM_CNV_ADDED = 0
        TOTAL_CNV = 0
        all_breakpoints_by_chromosome = {}
        for chrom in all_chromosomes:
            
            # We are going to check for CNVs that are not already at variant cut sites and only add the new ones
            
            # Then for each CNV we check if it is in a location not already close to a variant, if so we add it
            for CNV in breakpoints_from_copy_number.get(chrom,[]): # Copy number variant
                close_to_SRV = False
                TOTAL_CNV += 1
                for SRV in breakpoints_from_sniffles.get(chrom,[]): # Split read variant
                    if abs(CNV-SRV) < min_distance_to_variants_before_cutting_at_CNV:
                        close_to_SRV = True
                if close_to_SRV == False:
                    final_breakpoints_from_copy_number[chrom] = final_breakpoints_from_copy_number.get(chrom,[]) + [CNV]
                    NUM_CNV_ADDED += 1
            all_breakpoints_by_chromosome[chrom] = breakpoints_from_sniffles.get(chrom,[]) + final_breakpoints_from_copy_number.get(chrom,[])

        print "NUM_CNV_ADDED:", NUM_CNV_ADDED
        print "TOTAL_CNV:", TOTAL_CNV


        ####################################################

        # Split genome up into nodes using these breakpoints
        nodes = self.nodes_from_breakpoints(all_breakpoints_by_chromosome,genome_file)

        self.create_nodes_with_attributes(nodes)

        # Add node weights from the segmented copy number profile (these inform spanning edge weights for the CNV breakpoints)
        self.add_node_weights_from_seg_copy_number(coverage_file)

        ###############     Create Edges    #################
        # Read file again to connect the nodes with edges
        self.create_edges_from_sniffles(sniffles_filename)

        self.create_span_edges_at_CNVs(final_breakpoints_from_copy_number)

        return final_breakpoints_from_copy_number


    def add_node_weights_from_seg_copy_number(self,coverage_file,get_average_coverage=False):
        
        # dictionary of coverage (read all into memory), loop through each node to determine which coverage bins are on it, average them all

        f = open(coverage_file)
      
        interval_lengths = []

        segmented_coverage_by_chromosome = {}
        for line in f:
            fields = line.strip().split()
            if fields[1].isdigit(): # Ignore header 
                chrom = fields[0]
                start = float(fields[1])
                stop = float(fields[2])
                interval_lengths.append(stop-start)
                segmented_coverage = float(fields[4])
                segmented_coverage_by_chromosome[chrom] = segmented_coverage_by_chromosome.get(chrom, {})
                segmented_coverage_by_chromosome[chrom][start] = segmented_coverage                

        f.close()
        
        ################# Calculate weighted average coverage on each node ##############
        interval_length = int(interval_lengths[0])

        for node_name in self.nodes:
            node = self.nodes[node_name]
            
            chrom = node.attributes["chrom"]
            start = node.attributes["start"]
            stop = node.attributes["stop"]
            
            # Round to interval length
            seg_start = int(start - (start % interval_length))
            seg_stop = int(stop - (stop % interval_length))
            
            # For nodes so small they are contained within a single coverage bin, extend seg_stop just so the node can get that bin's coverage
            if seg_stop == seg_start:
                seg_stop = seg_start + interval_length

            segments_on_node = []
            for pos in xrange(seg_start,seg_stop,interval_length):
                segments_on_node.append(segmented_coverage_by_chromosome.get(chrom,{}).get(pos,0))

            node.weight = np.mean(segments_on_node)
            
        if get_average_coverage:
            coverage_on_nodes = []
            total_node_length = 0
            for node_name in self.nodes:
                node = self.nodes[node_name]
                length = (node.attributes["stop"]-node.attributes["start"])
                total_node_length += length
                coverage_on_nodes.append(node.weight*length)

            print "Average coverage on nodes:", sum(coverage_on_nodes)/total_node_length
            return sum(coverage_on_nodes)/total_node_length


    def select_spanning_port(self,port_object):
        span_ports = []
        for other_port in port_object.edges:
            edge = port_object.edges[other_port]
            if edge.spansplit == "span":
                span_ports.append(other_port)
        if len(span_ports) == 0:
            print "NO SPANNING EDGES FROM PORT", port_object, "(can still be CNV)"
            return None
        elif len(span_ports) > 1:
            print "MULTIPLE SPANNING EDGES FROM PORT", port_object, ":", span_ports
        return span_ports[0]



    def annotate_sniffles_file_with_variant_flow_evaluations(self,reports, sniffles_filename, output_filename):

        f = open(sniffles_filename)
        fout = open(output_filename, 'w')

        fout.write("chrom1,pos1,stop1,chrom2,pos2,stop2,variant_name,score,strand1,strand2,type,split,span1,span2,flow_category,description\n")
        for line in f:
            fields = line.strip().split()
            variant_name = fields[6]
            fout.write(",".join(fields) + "," + reports.get(variant_name, "VARIANT_NOT_IN_GRAPH") + "\n" )
        f.close()
        fout.close()
        
    def filter_sniffles_file_by_CNV_presence(self,reports, sniffles_filename, output_filename):
        f = open(sniffles_filename)
        fout = open(output_filename, 'w')

        num_variants_included = 0

        for line in f:
            fields = line.strip().split()
            variant_name = fields[6]
            score = reports.get(variant_name,"Missing")
            if score == "Missing":
                print "WARNING: Variant not found:", variant_name
            elif score.find("Perfect") >= 0 or score.find("Great") >= 0:
                num_variants_included += 1
                fout.write("\t".join(fields) + "\t" + reports.get(variant_name, "VARIANT_NOT_IN_GRAPH") + "\n" )

        f.close()
        fout.close()

    def segmented_coverage_to_CNV_calls_with_diff(self,coverage_file,cov_diff_threshold_to_split=0):
        # MUST BE SORTED BY CHROMOSOME THEN BY START POSITION

        # sample:
        # chromosome      start   end     unsegmented_coverage    coverage
        # 1       0       10000   1e-05   7.83570386725362
        # 1       10000   20000   1e-05   7.83570386725362

        f = open(coverage_file)

        breakpoints_by_chromosome = {}

        current_chromosome = ""
        current_coverage = -1

        for line in f:
            fields = line.strip().split()
            chrom = fields[0]
            if not fields[1].isdigit(): # Header 
                # print line.strip()
                pass
            else:
                start = int(fields[1])
                end = int(fields[2])
                # unsegmented_coverage = float(fields[3])
                segmented_coverage = float(fields[4])
                if chrom != current_chromosome:
                    current_coverage = segmented_coverage
                    current_chromosome = chrom
                    # You can't have breakpoints on the first bin, so continue to the next bin
                elif abs(segmented_coverage - current_coverage) > cov_diff_threshold_to_split:
                    breakpoints_by_chromosome[chrom] = breakpoints_by_chromosome.get(chrom, {})
                    breakpoints_by_chromosome[chrom][start] = (current_coverage,segmented_coverage)
                    current_coverage = segmented_coverage
        f.close()
                    
        return breakpoints_by_chromosome



    def average_coverage(self,coverage_file):
        f = open(coverage_file)

        total_coverage = 0
        total_bases = 0

        for line in f:
            fields = line.strip().split()
            chrom = fields[0]
            if not fields[1].isdigit(): # Header 
                # print line.strip()
                pass
            else:
                start = int(fields[1])
                end = int(fields[2])
                # unsegmented_coverage = float(fields[3])
                segmented_coverage = float(fields[4])
                total_coverage += segmented_coverage*(end-start)
                total_bases += end-start

        f.close()

        return total_coverage/total_bases


    def investigate_CNVs_for_quality(self,coverage_file,genome_file,threshold_for_long_CN_segments,coverage_threshold=None):
        
        breakpoints_from_copy_number = self.segmented_coverage_to_CNV_calls(coverage_file,cov_diff_threshold_to_split=0)


        if coverage_threshold == None:
            average_coverage = self.average_coverage(coverage_file)
            coverage_threshold = average_coverage/2

        nodes = self.nodes_from_breakpoints(breakpoints_from_copy_number,genome_file)

        self.create_nodes_with_attributes(nodes)

        self.add_node_weights_from_seg_copy_number(coverage_file)

        self.create_span_edges_at_CNVs(breakpoints_from_copy_number)

        CNV_quality = {}


        for edge in self.edges:
            if edge.spansplit == "CNV":
                
                port1 = edge.ports[0]
                port2 = edge.ports[1]

                node1 = port1.node
                node2 = port2.node

                port_name_1 = str(port1).split(":")[1]
                port_name_2 = str(port2).split(":")[1]
                if port_name_1 == "stop" and port_name_2 == "start":
                    pass
                else:
                    print "ERROR: CNV not left to right"

                if node1.attributes["stop"] != node2.attributes["start"]:
                    print "ERROR: node start and stop don't have the same positions across CNV edge"

                node1_length = node1.attributes["stop"] - node1.attributes["start"]
                node2_length = node2.attributes["stop"] - node2.attributes["start"]


                chrom = node1.attributes["chrom"]
                pos = node1.attributes["stop"]
                key = (chrom,pos)
                CNV_quality[key] = ""

                if node1_length > threshold_for_long_CN_segments and node2_length > threshold_for_long_CN_segments:
                    CNV_quality[key] += "Wide"
                else:
                    CNV_quality[key] += "Thin"

                if abs(node1.weight-node2.weight) > coverage_threshold:
                    CNV_quality[key] += "\tHigh"
                else:
                    CNV_quality[key] += "\tLow"


        return CNV_quality


    def output_CNVs_with_quality_and_SRV_concordance_analysis(self, coverage_file, CNV_features, CNVs_by_SRV_evidence, output_filename):

        breakpoints_from_copy_number = self.segmented_coverage_to_CNV_calls_with_diff(coverage_file,cov_diff_threshold_to_split=0)

        all_CNVs = set(CNV_features.keys() + CNVs_by_SRV_evidence.keys())
        
        f = open(output_filename,"w")
        f.write("chromosome\tposition\tleft_copy_number\tright_copy_number\tCN_change\tsegment_lengths\tCN_change_categorical\tSRV_evidence\n")

        for CNV in all_CNVs:
            chrom = CNV[0]
            pos = CNV[1]
            features = CNV_features.get(CNV,"No_info")
            SRV_evidence = CNVs_by_SRV_evidence.get(CNV,"No_info")

            left_CN = breakpoints_from_copy_number[chrom][pos][0]
            right_CN = breakpoints_from_copy_number[chrom][pos][1]
            
            f.write("%s\t%d\t%.2f\t%.2f\t%.2f\t%s\t%s\n" % (chrom,pos,left_CN,right_CN,abs(left_CN-right_CN),features,SRV_evidence))

        f.close()


    def score_CNVs(self,coverage_file, max_variant_CNV_distance):
        breakpoints_from_copy_number = self.segmented_coverage_to_CNV_calls_with_diff(coverage_file,cov_diff_threshold_to_split=0)

        CNV_category = {} # by (chrom,pos)

        ###  Needs to split up SRV into its two breakpoints and treat them separately 
        

        for chrom in breakpoints_from_copy_number:
            for CNV_pos in breakpoints_from_copy_number[chrom]:
                CNV_before_and_after = breakpoints_from_copy_number[chrom][CNV_pos]
                CNV_diff = CNV_before_and_after[1] - CNV_before_and_after[0]
                split_edges = self.split_edges_within_genome_interval(chrom, start=CNV_pos-max_variant_CNV_distance, end=CNV_pos+max_variant_CNV_distance)
                category = "Unknown"

                if len(split_edges) == 0:
                    category = "None"
                else:
                    num_SRVs_with_correct_direction = 0
                    for key in split_edges:
                        SRV_port = key[2]
                        # Check if direction of CN change matches direction of split reads
                        if SRV_port == "start" and CNV_diff > 0:
                            num_SRVs_with_correct_direction += 1
                        if SRV_port == "stop" and CNV_diff < 0:
                            num_SRVs_with_correct_direction += 1
                    if num_SRVs_with_correct_direction == 1:
                        category = "1"
                    elif num_SRVs_with_correct_direction == 0:
                        category = "None"
                    else:
                        category = "Many"
                # Save the score
                CNV_category[(chrom,CNV_pos)] = category
    
        return CNV_category

    def count_by_category(self,reports):
        counts_by_category = {}
        
        for variant in reports:
            category = reports[variant]
            counts_by_category[category] = counts_by_category.get(category, 0) + 1

        for category in counts_by_category:
            print category, ":\t", counts_by_category[category]
            
    def output_CNVs_with_flow_evaluations(self,CNV_category):
        for CNV in CNV_category:
            chrom = CNV[0]
            pos = CNV[1]
            flow_category = CNV_category[CNV]

            print "%s\t%d\t%s\n" % (chrom, pos, flow_category)

    def split_edges_within_genome_interval(self,chrom,start,end):
        nodes_near_CNV = self.nodes_within_genome_interval(chrom,start,end)
        split_edges = {}
        for node_name in nodes_near_CNV:
            node = self.nodes[node_name]
            if node.attributes["start"] > start:
                for port in node.ports["start"].edges:
                    edge = node.ports["start"].edges[port]
                    if edge.spansplit == "split":
                        split_edges[(node.attributes["chrom"],node.attributes["start"],"start")] = edge
            if node.attributes["stop"] < end:
                for port in node.ports["stop"].edges:
                    edge = node.ports["stop"].edges[port]
                    if edge.spansplit == "split":
                        split_edges[(node.attributes["chrom"],node.attributes["stop"],"stop")] = edge
        

        # split_edges is a dictionary with key: (chrom,pos,start/stop), value: edge object
        return split_edges

    def score_SRVs(self,coverage_file, max_variant_CNV_distance):

        breakpoints_from_copy_number = self.segmented_coverage_to_CNV_calls_with_diff(coverage_file,cov_diff_threshold_to_split=0)

        SRV_score_breakpoint_1 = {}
        SRV_score_breakpoint_2 = {}

        for edge in self.edges:
            if edge.spansplit == "split":
                split = edge.weight

                port1 = edge.ports[0]
                port2 = edge.ports[1]

                side_port1 = self.select_spanning_port(port1)
                side_port2 = self.select_spanning_port(port2)

                span1 = port1.edges[side_port1].weight
                span2 = port2.edges[side_port2].weight

                chrom1 = port1.node.attributes["chrom"]
                chrom2 = port2.node.attributes["chrom"]

                port_name_1 = str(port1).split(":")
                port_name_2 = str(port2).split(":")
                
                pos1 = self.nodes[port_name_1[0]].attributes[port_name_1[1]]
                pos2 = self.nodes[port_name_2[0]].attributes[port_name_2[1]]

                #  For each variant, check if it has a CNV nearby and whether that CNV is in the correct direction on both sides
                SRV_score_breakpoint_1[edge.variant_name] = 0
                SRV_score_breakpoint_2[edge.variant_name] = 0

                num_CNV_in_right_direction = 0
                for CNV in breakpoints_from_copy_number.get(chrom1,[]):
                    if abs(pos1 - CNV) < max_variant_CNV_distance:
                        SRV_score_breakpoint_1[edge.variant_name] = max([1,SRV_score_breakpoint_1[edge.variant_name]])
                        diff = breakpoints_from_copy_number[chrom1][CNV][1] - breakpoints_from_copy_number[chrom1][CNV][0]
                        if port_name_1[1] =="start" and diff > 0:
                            # print "1", diff
                            SRV_score_breakpoint_1[edge.variant_name] = 2
                            num_CNV_in_right_direction += 1
                        elif port_name_1[1] =="stop" and diff < 0:
                            # print "1", diff
                            SRV_score_breakpoint_1[edge.variant_name] = 2
                            num_CNV_in_right_direction +=1
                if num_CNV_in_right_direction == 1:
                    SRV_score_breakpoint_1[edge.variant_name] = 3

                num_CNV_in_right_direction = 0
                for CNV in breakpoints_from_copy_number.get(chrom2,[]):
                    if abs(pos2 - CNV) < max_variant_CNV_distance:
                        SRV_score_breakpoint_2[edge.variant_name] = max([1,SRV_score_breakpoint_2[edge.variant_name]])
                        diff = breakpoints_from_copy_number[chrom2][CNV][1] - breakpoints_from_copy_number[chrom2][CNV][0]
                        if port_name_2[1] =="start" and diff > 0:
                            # print "2", diff
                            SRV_score_breakpoint_2[edge.variant_name] = 2
                            num_CNV_in_right_direction +=1
                        elif port_name_2[1] =="stop" and diff < 0:
                            # print "2", diff
                            SRV_score_breakpoint_2[edge.variant_name] = 2
                            num_CNV_in_right_direction +=1
                if num_CNV_in_right_direction == 1:
                    SRV_score_breakpoint_2[edge.variant_name] = 3


        reports = {}
        for variant_name in SRV_score_breakpoint_1:
            reports[variant_name] = (SRV_score_breakpoint_1[variant_name], SRV_score_breakpoint_2[variant_name])

        return reports


    def summarize_SRVs(self,SRVs):
        # SRVs is a dictionary with key = variant_name, value = tuple = (score for breakpoint 1, score for breakpoint 2)

        # summary is a dictionary with key = variant_name, value = string describing variant
        summary = {}
        for variant_name in SRVs:
            scores = SRVs[variant_name]
            s1 = scores[0]
            s2 = scores[1]
            if s1 == 3 and s2 == 3:
                summary[variant_name] = "Perfect,exactly 1 CNV on each side"
            elif s1 >= 2 and s2 >= 2:
                summary[variant_name] = "Great,1 or more explanatory CNVs"
            elif s1 < 2 and s2 < 2:
                summary[variant_name] = "Bad,CNVs missing or in wrong direction"
            else:
                summary[variant_name] = "Poor,explanatory CNVs only on one side"

        return summary


    def category_counts(self,reports):
        counts_by_category = {}
        for variant in reports:
            category = reports[variant]
            counts_by_category[category] = counts_by_category.get(category, 0) + 1

        return counts_by_category


    def calculate_flow_on_node(self,node_name):
        node = self.nodes[node_name]


        port1 = node.ports["start"]
        port2 = node.ports["stop"]

        sum_port1_edgeweights = 0
        for name in port1.edges:
            edge = port1.edges[name]
            sum_port1_edgeweights += edge.weight
        
        sum_port2_edgeweights = 0
        for name in port2.edges:
            edge = port2.edges[name]
            sum_port2_edgeweights += edge.weight

        return sum_port1_edgeweights, sum_port2_edgeweights



    def output_flow_on_nodes(self,output_filename):

        f = open(output_filename,'w')
        f.write("node_name,chromosome,start,stop,flow_in,flow_out,flow_balance,node_weight,in_matches_weight,out_matches_weight\n")

        for node_name in self.nodes:
            node_weight = self.nodes[node_name].weight
            sum_port1_edgeweights, sum_port2_edgeweights = self.calculate_flow_on_node(node_name)
            attributes = self.nodes[node_name].attributes
            flow_balance = self.approx_equal(sum_port1_edgeweights, sum_port2_edgeweights)
            in_matches_weight = self.approx_equal(sum_port1_edgeweights, node_weight)
            out_matches_weight = self.approx_equal(sum_port2_edgeweights, node_weight)
            f.write("%s,%s,%d,%d,%.2f,%.2f,%s,%.2f,in matches CN: %s,out matches CN: %s\n" % (node_name, attributes["chrom"],attributes["start"],attributes["stop"],sum_port1_edgeweights, sum_port2_edgeweights,flow_balance, node_weight,in_matches_weight,out_matches_weight))

        f.close()


    # def find_suitable_regions(self,CNV_breakpoints,genome_file,verbose=False):

    #     regions = self.nodes_from_breakpoints(CNV_breakpoints,genome_file)

    #     suitable_regions = []
    #     for region_name in regions:
    #         attributes = regions[region_name]
    #         chrom = attributes["chrom"]
    #         start = attributes["start"]
    #         stop = attributes["stop"]
    #         size = stop - start

    #         split_edges = self.split_edges_within_genome_interval(chrom,start,stop)
    #         num_splits = len(split_edges)


    #         if size > 1000000 and num_splits > 3 and num_splits < 1000:
    #             if verbose: print chrom,"\t", "Splits:", num_splits, "Size:", size
    #             profile = {}
    #             for key in attributes:
    #                 profile[key] = attributes[key]
    #             profile["num_splits"] = num_splits

    #             suitable_regions.append(profile)

    #     print "Found", len(suitable_regions), "suitable regions for Evolution"

    #     return suitable_regions


    def approx_equal(self,num1,num2):

        if abs(num1 - num2) < (num1+num2)/2:
            return True
        else:
            return False

    def significantly_smaller(self,num1,num2):

        if num1 < (num2 - absolute_difference):
            return True
        else:
            return False

    # def evaluate_regions(self,suitable_regions,output_filename):
    #     print 'UNFINISHED FUNCTION: evaluate_regions'
        
    #     for region in suitable_regions:
    #         node_names = self.nodes_within_genome_interval(region["chrom"],region["start"],region["stop"])
    #         node_balance = []
    #         for node_name in node_names:
                # node_weight = self.nodes[node_name].weight
                # sum_port1_edgeweights, sum_port2_edgeweights = self.calculate_flow_on_node(node_name)
                # attributes = self.nodes[node_name].attributes
                # flow_balance = self.approx_equal(sum_port1_edgeweights, sum_port2_edgeweights)
                # in_matches_weight = self.approx_equal(sum_port1_edgeweights, node_weight)
                # out_matches_weight = self.approx_equal(sum_port2_edgeweights, node_weight)

                # in_vs_out = self.approx_equal(sum_port1_edgeweights,sum_port2_edgeweights):
                # in_vs_node = self.approx_equal(sum_port1_edgeweights,node_weight)
                # out_vs_node = self.approx_equal(sum_port2_edgeweights,node_weight)

                # # If all 3 == True:
                #     # No rebalancing needed, but we can do a refinement to make the values closer
                # # If 1 == True:
                #     # Set the one that doesn't agree to an average of the other two
                # # If none agree:
                #     # Red flag the region and ignore it or cut the region into pieces at this point and start over



    # def calculate_flow_on_node(self,node_name):
    #     node = self.nodes[node_name]


    #     port1 = node.ports["start"]
    #     port2 = node.ports["stop"]

    #     sum_port1_edgeweights = 0
    #     for name in port1.edges:
    #         edge = port1.edges[name]
    #         sum_port1_edgeweights += edge.weight
        
    #     sum_port2_edgeweights = 0
    #     for name in port2.edges:
    #         edge = port2.edges[name]
    #         sum_port2_edgeweights += edge.weight

    #     return sum_port1_edgeweights, sum_port2_edgeweights


    def balance_port(self,port,increase):

        # if increase > 0: ################################################
        found = False
        for name in port.edges:
            edge = port.edges[name]
            if edge.spansplit == "split":
                edge.weight = max([edge.weight + increase,0])
                found = True
                break
        if found == False:
            for name in port.edges:
                edge = port.edges[name]
                # if edge.spansplit == "span":
                edge.weight = max([edge.weight + increase,0])
                found = True
                break

    def all_equal(self,values,allowance=0.001):
        for item in values:
            if abs(item - values[0])>allowance:
                return False
        return True

    def balance_flow_on_nodes(self,strategy = "average"):

        if not strategy in ["CN", "average"]:
            print "strategy must be 'CN' or 'average'"
            return []

        categories = {}
        missing_flow_values = {}

        unbalanceable_nodes = []

        for node_name in self.nodes:
            node = self.nodes[node_name]
            node_weight = self.nodes[node_name].weight

            sum_port1_edgeweights, sum_port2_edgeweights = self.calculate_flow_on_node(node_name)
            attributes = node.attributes
            flow_balance = self.approx_equal(sum_port1_edgeweights, sum_port2_edgeweights)
            in_matches_weight = self.approx_equal(sum_port1_edgeweights, node_weight)
            out_matches_weight = self.approx_equal(sum_port2_edgeweights, node_weight)

            in_vs_out = self.approx_equal(sum_port1_edgeweights,sum_port2_edgeweights)
            in_vs_node = self.approx_equal(sum_port1_edgeweights,node_weight)
            out_vs_node = self.approx_equal(sum_port2_edgeweights,node_weight)

            categories[node_name] = (in_vs_out,in_vs_node,out_vs_node)

            # Get ready to edit some of these values
            port1 = node.ports["start"]
            port2 = node.ports["stop"]

            # If all 3 == True:
                # No rebalancing needed, but we can do a refinement to make the values closer
            # If 1 == True:
                # Set the one that doesn't agree to an average of the other two
            # If none agree:
                # Red flag the region and ignore it or cut the region into pieces at this point and start over

            if sum([in_vs_out,in_vs_node,out_vs_node]) > 0:
                new_flow_level = node_weight
                if node_weight == sum_port1_edgeweights and node_weight == sum_port2_edgeweights:
                    categories[node_name] = "All exactly equal"
                    new_flow_level = node_weight
                elif (in_vs_out,in_vs_node,out_vs_node) == (True,True,True):
                    categories[node_name] = "All agree, refined"
                    new_flow_level = node_weight
                elif (in_vs_out,in_vs_node,out_vs_node) == (True,False,False):
                    categories[node_name] = "In and out agree"
                    new_flow_level = (sum_port1_edgeweights + sum_port2_edgeweights)/2
                elif (in_vs_out,in_vs_node,out_vs_node) == (False,True,False):
                    categories[node_name] = "Node and in agree"
                    if strategy == "average":
                        new_flow_level = (sum_port1_edgeweights + node_weight)/2
                    elif strategy == "CN":
                        new_flow_level = node_weight
                    else:
                        print "UNKNOWN STRATEGY"
                elif (in_vs_out,in_vs_node,out_vs_node) == (False,False,True):
                    categories[node_name] = "Node and out agree"
                    if strategy == "average":
                        new_flow_level = (sum_port2_edgeweights + node_weight)/2
                    elif strategy == "CN":
                        new_flow_level = node_weight
                    else:
                        print "UNKNOWN STRATEGY"
                elif (in_vs_out,in_vs_node,out_vs_node) == (True,True,False):
                    categories[node_name] = "In agrees with both"
                    new_flow_level = sum_port1_edgeweights
                elif (in_vs_out,in_vs_node,out_vs_node) == (True,False,True):
                    categories[node_name] = "Out agrees with both"
                    new_flow_level = sum_port2_edgeweights
                elif (in_vs_out,in_vs_node,out_vs_node) == (False,True,True):
                    categories[node_name] = "Node agrees with both"
                    new_flow_level = node_weight
                else:
                    print (in_vs_out,in_vs_node,out_vs_node), "MISSING"
                
                # Rebalance to the new_flow_level
                if new_flow_level < 0:
                    new_flow_level = 0
                self.balance_port(port1,new_flow_level - sum_port1_edgeweights)
                self.balance_port(port2,new_flow_level - sum_port2_edgeweights)
                self.nodes[node_name].weight = new_flow_level
            else:
                categories[node_name] = "All disagree, red flag this node"
                unbalanceable_nodes.append(node_name)

            # After all attempts at balancing, quantify exactly how much each node has a discrepancy in flow left over
            sum_port1_edgeweights, sum_port2_edgeweights = self.calculate_flow_on_node(node_name)
            biggest_diff = max([abs(node_weight - sum_port1_edgeweights), abs(node_weight - sum_port2_edgeweights), abs(sum_port1_edgeweights - sum_port2_edgeweights)])
            missing_flow_values[node_name] = biggest_diff

        print "sum missing flow:", sum(missing_flow_values.values())

        self.count_by_category(categories)

        return unbalanceable_nodes
        


    def count_balanced_nodes(self,verbose=False):

        num_balanced_nodes = 0
        num_total_nodes = 0

        for node_name in self.nodes:
            node = self.nodes[node_name]
            node_weight = self.nodes[node_name].weight

            sum_port1_edgeweights, sum_port2_edgeweights = self.calculate_flow_on_node(node_name)
            attributes = node.attributes
            flow_balance = self.approx_equal(sum_port1_edgeweights, sum_port2_edgeweights)
            in_matches_weight = self.approx_equal(sum_port1_edgeweights, node_weight)
            out_matches_weight = self.approx_equal(sum_port2_edgeweights, node_weight)

            in_vs_out = self.approx_equal(sum_port1_edgeweights,sum_port2_edgeweights)
            in_vs_node = self.approx_equal(sum_port1_edgeweights,node_weight)
            out_vs_node = self.approx_equal(sum_port2_edgeweights,node_weight)

            if sum([in_vs_out,in_vs_node,out_vs_node]) > 1:
                num_balanced_nodes += 1
            num_total_nodes += 1

        if verbose:
            print "Number of balanced nodes:", num_balanced_nodes
            print "Total number of nodes:", num_total_nodes

        return num_balanced_nodes, num_total_nodes



    def output_variants_from_graph(self,output_filename,flow_reports=None):
        fout = open(output_filename, 'w')
        fout.write("chrom1,pos1,chrom2,pos2,variant_name,strand1,strand2,split,span1,span2,flow_category,description\n")
        
        for edge in self.edges:
            if edge.spansplit == "split":
                split = edge.weight
                variant_name = edge.variant_name

                port1 = edge.ports[0]
                port2 = edge.ports[1]

                side_port1 = self.select_spanning_port(port1)
                side_port2 = self.select_spanning_port(port2)

                span1 = port1.edges[side_port1].weight
                span2 = port2.edges[side_port2].weight

                chrom1 = port1.node.attributes["chrom"]
                chrom2 = port2.node.attributes["chrom"]

                port_name_1 = str(port1).split(":")
                port_name_2 = str(port2).split(":")
                
                pos1 = self.nodes[port_name_1[0]].attributes[port_name_1[1]]
                pos2 = self.nodes[port_name_2[0]].attributes[port_name_2[1]]

                strand1 = "+"
                if port_name_1[1] == "start":
                    strand1 = "-"
                strand2 = "+"
                if port_name_2[1] == "start":
                    strand2 = "-"

                flow_info = "No flow info"
                if flow_reports != None:
                    flow_info = flow_reports.get(variant_name, "No flow info")
                
                fout.write("%s,%d,%s,%d,%s,%s,%s,%.2f,%.2f,%.2f,%s\n" % (chrom1,pos1,chrom2,pos2,variant_name,strand1,strand2,split,span1,span2,flow_info))

        fout.close()



    def output_nodes_as_boxes(self,output_filename):
        
        f = open(output_filename,'w')
        f.write("chromosome,start,end,y_start,height,path_ID\n")

        path_counter = 1
        for node_name in self.nodes:
            node = self.nodes[node_name]
            f.write("%s,%d,%d,%d,%d,path_%d\n" % (node.attributes["chrom"],node.attributes["start"], node.attributes["stop"],0,node.weight,path_counter))
            path_counter += 1
        f.close()


    def local_parsimony(self,chrom,start,end):
        primary_nodes = self.nodes_within_genome_interval(chrom,start,end)
        s1 = self.special_subgraph_from_genome_interval(chrom,start,end,degree=5)
        # print s1

        # s1.tree_parsimony()


    def node_flow_diff(self,node):
        sum_port1_edgeweights, sum_port2_edgeweights = self.calculate_flow_on_node(node.name)

        diff_flow = abs(sum_port1_edgeweights-sum_port2_edgeweights) + abs((sum_port1_edgeweights+sum_port2_edgeweights)/2 - node.weight)
        return diff_flow

    def score_variants_by_contribution_to_flow_on_nodes(self):

        variant_scores = {}

        for node_name in self.nodes:
            node = self.nodes[node_name]
            node_weight = self.nodes[node_name].weight

            sum_port1_edgeweights, sum_port2_edgeweights = self.calculate_flow_on_node(node_name)
            attributes = node.attributes
            flow_balance = self.approx_equal(sum_port1_edgeweights, sum_port2_edgeweights)
            in_matches_weight = self.approx_equal(sum_port1_edgeweights, node_weight)
            out_matches_weight = self.approx_equal(sum_port2_edgeweights, node_weight)

            in_vs_out = self.approx_equal(sum_port1_edgeweights,sum_port2_edgeweights)
            in_vs_node = self.approx_equal(sum_port1_edgeweights,node_weight)
            out_vs_node = self.approx_equal(sum_port2_edgeweights,node_weight)

            port1 = node.ports["start"]
            port2 = node.ports["stop"]

            original_node_flow_diff = self.node_flow_diff(node)


            if sum([in_vs_out,in_vs_node,out_vs_node]) > 1:
                for port in [port1,port2]:
                    for name in port.edges:
                        edge = port.edges[name]
                        if edge.spansplit == "split":
                            saved_edge_weight = edge.weight
                            edge.weight = 0
                            experimental_node_flow_diff = self.node_flow_diff(node)
                            edge.weight = saved_edge_weight
                            if original_node_flow_diff > experimental_node_flow_diff:
                                variant_scores[edge.variant_name] = max([variant_scores.get(edge.variant_name,0),100])
                            else:
                                variant_scores[edge.variant_name] = max([variant_scores.get(edge.variant_name,0),0])

        self.count_by_category(variant_scores)


    # def tree_parsimony(self):
    #     portal_name = "Portal"
    #     cycle_limit = 2
    #     depth_limit = 30


    #     self.add_portal(all_ports = False)

    #     print self.edges

        # before = time.time()
        # allpaths = self.depth_first_search(self.nodes[portal_name].ports["start"],self.nodes[portal_name].ports["start"],cycle_limit=cycle_limit,depth_limit=depth_limit)
        
        # print len(allpaths)
        
        # unique_paths = self.collapse_redundant_paths(allpaths)

        # print len(unique_paths)
        # for path in unique_paths:
        #     print path



    def show_flow_around_gene(self,gene,output_filename):
        annotation = self.annotation.get(gene)

        print annotation

        # Traverse the graph in both directions and go down all possible paths up to a certain distance 
        step = 200000
        checkpoint_locations = range(step,step*10,step)

        start_chrom = annotation["chrom"]
        start_pos = (annotation["start"]+annotation["stop"])/2

        node_name = self.find_nodename_by_position((start_chrom,start_pos))
        node = self.nodes[node_name]
        attr = node.attributes
        node_sequence_length = attr["stop"] - attr["start"]


        forward_and_reverse_flows = []
        forward_and_reverse_flows.append((self.nodes[node_name].attributes["chrom"],0))

        for start_direction in ["+","-"]:
            
            # Find all the positions of the threads after "distance" distance
            start_port = None
            current_distance = 0

            if start_direction == "+":
                start_port = self.nodes[node_name].ports["start"]
                current_distance = attr["start"]- start_pos # negative to offset the fact that we are gonna travel across the whole node
            elif start_direction == "-":
                start_port = self.nodes[node_name].ports["stop"]
                current_distance = start_pos - attr["stop"] # negative to offset the fact that we are gonna travel across the whole node
            else:
                print "ERROR: start_direction must be + or -"
                return

            print "start_port", start_port

            all_places = []
            self.depth_first_traverse(current_port=start_port, current_distance=current_distance, checkpoint_index=0,all_places=all_places,checkpoint_locations=checkpoint_locations)
            if start_direction == "+":
                forward_and_reverse_flows.extend(all_places)
            else:
                for place in all_places:
                    place = (place[0],-1*place[1])
                    forward_and_reverse_flows.append(place)


        f = open(output_filename,'w')
        f.write("chromosome,distance\n")
        for place in forward_and_reverse_flows:
            chrom = place[0]
            distance = place[1]

            f.write("%s,%d\n" % (chrom,distance))

        f.close()



    def depth_first_traverse(self,current_port,current_distance,checkpoint_index,all_places,checkpoint_locations,depth=0,max_recursion_depth = 20):
        target_distance = checkpoint_locations[checkpoint_index]

        attr = current_port.node.attributes
        node_sequence_length = attr["stop"] - attr["start"]

        current_distance += node_sequence_length

        jumped_port = current_port.jump()
        
        saveport = str(jumped_port)

        while current_distance >= target_distance:
            # all_places.append((current_port,target_distance))
            all_places.append((current_port.node.attributes["chrom"],target_distance))
            checkpoint_index += 1
        
            if checkpoint_index >= len(checkpoint_locations):
                return
            target_distance = checkpoint_locations[checkpoint_index]                
        
        if depth >= max_recursion_depth:
            # print "hit max recursion"
            all_places.append(("maxed out",target_distance))
            return
        # keep recursing
        edges = jumped_port.edges.values()
        for edge in edges:
            glide_port = edge.glide(jumped_port)
            self.depth_first_traverse(current_port=glide_port, current_distance=current_distance, checkpoint_index=checkpoint_index,all_places=all_places,checkpoint_locations=checkpoint_locations,depth=depth+1,max_recursion_depth=max_recursion_depth)


########################  Segment strategy for local evolution  ##########################

    def bottom_up_evolution(self,chromosome,start,end,output_prefix,min_weight_required):
        
        # 1. Find segments
        segments = self.find_segments_from_genome_interval(chromosome,start,end,min_weight_required=min_weight_required)
        print segments
        print "Found a set of", len(segments), "segments."
        
        self.boxes_from_parsimony(segments,output_filename=output_prefix+".boxes.csv")

        # 2. Link segments but remember to impose Parsimony!!
            # Also make sure inverted duplications work out correctly (currently may be counting only half and the looping edge was cut out of the subgraph when we found segments)

        # 3. That's it! Display the results


    def spanning_only_subgraph_from_genome_interval(self,chromosome,start,end,verbose=True):
        primary_nodes = self.nodes_within_genome_interval(chromosome,start,end)

        s = Graph()

        node_attributes = {}
        edges_to_add = {} # Create a dictionary of edges to avoid adding the same edge multiple times
        portals = {}

        for node_name in primary_nodes:
            # Copy the node itself
            node_attributes[node_name] = self.nodes[node_name].attributes
            # Find all the first-degree nodes and their edges
            for port_name in self.nodes[node_name].ports:
                this_port = self.nodes[node_name].ports[port_name]
                edges = this_port.edges
                if len(edges) == 0:
                    # print "NO EDGES:", port_name
                    portals[self.nodes[node_name].ports[port_name]] = self.nodes[node_name].weight
                for other_port in edges:
                    edge = edges[other_port]

                    # If this is a spanning edge within the region:
                    within_region = str(other_port.node.name) in map(str,primary_nodes)
                    
                    if (within_region and edge.spansplit in ["CNV","span"]) or this_port==other_port:
                        edges_to_add[str(edge)] = edge
                    else:
                        # Add a portal with the weight of the edge we are not including (sum these for the special case on the edge of the region where flow from both span and split edges need to escape to the portal )
                        portals[self.nodes[node_name].ports[port_name]] = portals.get(self.nodes[node_name].ports[port_name],0) + edge.weight

        # Add nodes to the graph
        s.create_nodes_with_attributes(node_attributes)

        for node_name in node_attributes:
            s.nodes[node_name].weight = self.nodes[node_name].weight

        edge_counter = 0

        # Add edges to the graph
        for edge_name in edges_to_add:
            edge = edges_to_add[edge_name]
            node1 = edge.ports[0].node.name
            port1 = edge.ports[0].name
            node2 = edge.ports[1].node.name
            port2 = edge.ports[1].name
            weight = edge.weight
            spansplit = edge.spansplit
            s.create_edge(node1,port1,node2,port2,weight=weight,spansplit=spansplit)
            edge_counter += 1


        # Add portal 
        s.create_nodes_with_attributes({"Portal":{"chrom":"0","start":0,"stop":0,"x":0,"y":0}})
        for dead_end_port in portals:
            s.create_edge(dead_end_port.node.name, dead_end_port.name, "Portal", "stop", portals[dead_end_port],spansplit="portal")

        return s


    def find_segments_from_genome_interval(self,chromosome,start,end,verbose = True,min_weight_required=0):

        # Build subgraph of primary nodes only (spanning edges only) and hook all split edges up to the Portal:
        s = self.spanning_only_subgraph_from_genome_interval(chromosome,start,end)   

        if verbose:
            print s
            s.print_edges()
        
        # Find all spanning paths
        # Do the usual iterative subtraction until all paths are gone
        
        reports = s.parsimony(depth_limit = 50,verbose=True,chop_end_nodes=0,min_weight_required=min_weight_required)
        if verbose:
            for report in reports:
                print report

        # Report these segments
        return reports


    def cover_uncovered_nodes(self):
        counter = 0
        for node_name in self.nodes:
            node = self.nodes[node_name]
            port1 = node.ports["start"]
            port2 = node.ports["stop"]

            weight_side_1 = 0
            weight_side_2 = 0

            CNV_on_side_1 = False
            for other_port in port1.edges:
                if port1.edges[other_port].spansplit == "CNV":
                    CNV_on_side_1 = True
                    weight_side_1 = other_port.node.weight
                    break

            CNV_on_side_2 = False
            for other_port in port2.edges:
                if port2.edges[other_port].spansplit == "CNV":
                    CNV_on_side_2 = True
                    weight_side_2 = other_port.node.weight
                    break

            if CNV_on_side_1 and CNV_on_side_2:
                # print node_name, "has CNV on both sides"
                # print weight_side_1, node.weight, weight_side_2
                if node.weight < min([weight_side_1,weight_side_2]):
                    node.weight = min([weight_side_1,weight_side_2])
                counter += 1
        print counter, "total nodes with CNVs on both sides"












    ########################################################################################################
    ################### All major functions needed for April 6th version below this line ###################
    ########################################################################################################
    # coverage, variants, and genome files are now all .csv files
    # span counts are now set by coverage only

    def bin_size(self,coverage_file):
        f = open(coverage_file)

        bin_sizes = []

        for line in f:
            fields = line.strip().split(",")
            chrom = fields[0]
            if not fields[1].isdigit(): # Header 
                # print line.strip()
                pass
            else:
                start = int(fields[1])
                end = int(fields[2])
                bin_sizes.append(end-start)

        f.close()

        return np.median(bin_sizes)


    def segmented_coverage_to_CNV_calls(self,coverage_file,cov_diff_threshold_to_split=0,csv=False):
        # MUST BE SORTED BY CHROMOSOME THEN BY START POSITION

        # sample:
        # chromosome      start   end     unsegmented_coverage    coverage
        # 1       0       10000   1e-05   7.83570386725362
        # 1       10000   20000   1e-05   7.83570386725362

        f = open(coverage_file)

        breakpoints_by_chromosome = {}

        current_chromosome = ""
        current_coverage = -1

        for line in f:
            fields = line.strip().split()
            if csv:
                fields = line.strip().split(",")
            chrom = fields[0]
            if not fields[1].isdigit(): # Header 
                continue
            start = int(fields[1])
            end = int(fields[2])
            # unsegmented_coverage = float(fields[3])
            segmented_coverage = float(fields[4])
            if chrom != current_chromosome:
                current_coverage = segmented_coverage
                current_chromosome = chrom
                # You can't have breakpoints on the first bin, so continue to the next bin
                continue
            if abs(segmented_coverage - current_coverage) > cov_diff_threshold_to_split:
                current_coverage = segmented_coverage
                breakpoints_by_chromosome[chrom] = breakpoints_by_chromosome.get(chrom, []) + [start]

        return breakpoints_by_chromosome

    def breakpoints_from_variants(self,variants_filename,csv=False):
        # Read file once to note all the breakpoints
        breakpoints_by_chromosome = {}
        f = open(variants_filename)
        for line in f:
            fields = line.strip().split()
            if csv:
                fields = line.strip().split(",")
            if not fields[1].isdigit(): # Header 
                continue
            chrom1 = fields[0]
            pos1 = (int(fields[1]) + int(fields[2]))/2
            breakpoints_by_chromosome[chrom1] = breakpoints_by_chromosome.get(chrom1,[]) + [pos1]
            chrom2 = fields[3]
            pos2 = (int(fields[4]) + int(fields[5]))/2
            breakpoints_by_chromosome[chrom2] = breakpoints_by_chromosome.get(chrom2,[]) + [pos2]

        f.close()
        return breakpoints_by_chromosome

    def binary_search_for_closest(self,x,a,start,end):
        # search for the closest position to x
        # O(log(m)) where m = number of CNVs

        if end-start <= 1:
            if abs(x-a[start]) <= abs(x-a[end]):
                return a[start]
            else:
                return a[end]
        else:
            mid = (start+end)/2
            if x > a[mid]:
                return self.binary_search_for_closest(x,a,start=mid,end=end)
            else:
                return self.binary_search_for_closest(x,a,start=start,end=mid)

    def snap_variants_to_CNVs(self, variant_breakpoints, copy_number_breakpoints):
        if len(copy_number_breakpoints)==0:
            variants_with_nearest_breakpoint = []
            for variant in variant_breakpoints:
                variants_with_nearest_breakpoint.append((variant,None))
            return variants_with_nearest_breakpoint
        # O(n*log(m)) where n = number of variants, m = number of CNVs
        # O(n*log(m)) = O(n) for looping through variants, O(log(m)) for binary search through CNVs

        # Sort these in case the original input coverage file was not sorted
        CNVs = np.array(copy_number_breakpoints)
        CNVs.sort()

        variants_with_nearest_breakpoint = []
        for variant in variant_breakpoints:
            variants_with_nearest_breakpoint.append((variant, self.binary_search_for_closest(variant,CNVs,0,len(CNVs)-1)))

        #  return tuple with (variant,nearest_CNV), just positions, no chromosomes
        return variants_with_nearest_breakpoint 


    def read_genome_file(self,genome_file,csv=False):
        f = open(genome_file)
        chromosome_lengths = {}
        for line in f:
            fields = line.strip().split()
            if csv:
                fields = line.strip().split(",")
            if not fields[1].isdigit(): # Header 
                continue
            chromosome_lengths[fields[0]] = int(fields[1])
        f.close()
        return chromosome_lengths

    def nodes_from_breakpoints(self,breakpoints_by_chromosome,genome_file,csv=False):
        # INPUT: breakpoints_by_chromosome is a dictionary with keys of chromosomes and entries as lists of positions to cut at

        # Use genome file to find ends of chromosomes
        chromosome_lengths = self.read_genome_file(genome_file,csv)

        nodes={}
        for chrom in chromosome_lengths:
            breakpoints = breakpoints_by_chromosome.get(chrom,[])
            breakpoints.sort()
            counter = 1
            node_start = 0
            for breakpoint in breakpoints:
                if breakpoint != node_start:
                    node_name = "%s.%d" % (chrom,counter)
                    nodes[node_name] = {"chrom":chrom,"start":node_start,"stop":breakpoint}
                    node_start = breakpoint
                    counter += 1
            node_name = "%s.%d" % (chrom,counter)
            nodes[node_name] = {"chrom":chrom,"start":node_start,"stop":chromosome_lengths[chrom]}

        return nodes


    def snap_breakpoints(self, breakpoints_from_copy_number, all_breakpoints_from_variants, resolution,verbose=False):

        variant_cluster_consensus_breakpoints = {}
        for chrom in all_breakpoints_from_variants:
            # Find which variants can snap onto CNVs and which ones cannot
            unsnapped_variants = []
            snapped_variants = []
            if verbose: print "chromosome:", chrom 
            variants_with_nearest_breakpoint = self.snap_variants_to_CNVs(variant_breakpoints=all_breakpoints_from_variants[chrom],copy_number_breakpoints=breakpoints_from_copy_number.get(chrom,[]))
            for var,CNV in variants_with_nearest_breakpoint:
                if CNV == None:
                    unsnapped_variants.append(var)
                elif abs(var-CNV) < resolution:
                    snapped_variants.append((var,CNV))
                else:
                    unsnapped_variants.append(var)
            if verbose: 
                print "snapped", len(snapped_variants)
                print "unsnapped", len(unsnapped_variants)
                print "___________________________"

            # For variants that cannot snap onto CNVs:
            # Cluster variants within our resolution and choose the average position as the breakpoint
            sorted_positions = np.array(unsnapped_variants)
            sorted_positions.sort()
            clusters = []
            current_cluster = []

            for i in xrange(len(sorted_positions)):
                if current_cluster == [] or abs(current_cluster[-1]-sorted_positions[i]) < resolution:
                    current_cluster.append(sorted_positions[i])
                else:
                    clusters.append(current_cluster)
                    current_cluster = [sorted_positions[i]]
            clusters.append(current_cluster)

            for cluster in clusters:
                variant_cluster_consensus_breakpoints[chrom] = variant_cluster_consensus_breakpoints.get(chrom,[]) + [int(np.mean(cluster))]

        # Cut the graph at all CNV and clustered variant breakpoints
        consolidated_breakpoints_by_chromosome = {}
        for chrom in set(variant_cluster_consensus_breakpoints.keys()).union(set(breakpoints_from_copy_number.keys())):
            consolidated_breakpoints_by_chromosome[chrom] = variant_cluster_consensus_breakpoints.get(chrom,[]) + breakpoints_from_copy_number.get(chrom,[])

        return consolidated_breakpoints_by_chromosome


    def breakpoints_from_graph(self):
        # output: dictionary of breakpoints[chrom]
        breakpoints = {}
        for node_name in self.nodes:
            attributes = self.nodes[node_name].attributes
            chrom = attributes["chrom"]
            breakpoints[chrom] = breakpoints.get(chrom,[]) + [attributes["start"]] + [attributes["stop"]]

        return breakpoints

    def create_split_edge(self, key1, key2, weight_split,split_variant_name,verbose=False):
        port1 = ""
        if key1[2] == "+":
            port1 = "stop"
        elif key1[2] == "-":
            port1 = "start"
        else:
            raise FileFormatError("strands must be + or - in columns 9 and 10")

        port2 = ""
        if key2[2] == "+":
            port2 = "stop"
        elif key2[2] == "-":
            port2 = "start"
        else:
            raise FileFormatError("strands must be + or - in columns 9 and 10")

        # split
        node1 = self.node_lookup.get(key1, "NA")
        node2 = self.node_lookup.get(key2, "NA")
        if node1 == "NA" or node2 == "NA":
            print "Keys attempted"
            print "key1:", key1, "found:", node1
            print "key2:", key2, "found:", node2
            print "___________________________________"
            
            if node1 == "NA":
                print "Dictionary for chrom of key1:"
                for key in self.node_lookup:
                    if key[0] == key1[0]:
                        print key
            if node2 == "NA":
                print "Dictionary for chrom of key2:"
                for key in self.node_lookup:
                    if key[0] == key2[0]:
                        print key
            raise FileFormatError("Attempted to create edge at breakpoint not represented among the nodes in the graph")
        self.create_edge(node1,port1,node2,port2,weight_split,spansplit="split",split_variant_name=split_variant_name)
        if verbose==True:
            print "nodes:", node1, node2

    def connect_nodes_at_variants(self,variants_filename,resolution,verbose=False):

        # get all breakpoints from graph
        graph_breakpoints = self.breakpoints_from_graph() # output: dictionary[chrom]
        for chrom in graph_breakpoints:
            graph_breakpoints[chrom] = np.array(graph_breakpoints[chrom])
            graph_breakpoints[chrom].sort()
            graph_breakpoints[chrom] = graph_breakpoints[chrom][1:-1] # avoid ends of chromosomes

        # read variants file
        # for variant in file:
        #       snap onto breakpoints
        #       create split edge

        f = open(variants_filename)
        # count_far_away = 0
        variant_names_used_already = set()

        variant_counter = 0
        for line in f:
            fields = line.strip().split(",")
            if not fields[1].isdigit(): # Header 
                continue
            chrom1 = fields[0]
            pos1 = (int(fields[1]) + int(fields[2]))/2
            closest1 = self.binary_search_for_closest(pos1,graph_breakpoints[chrom1],0,len(graph_breakpoints[chrom1])-1)
            strand1 = fields[8]
            # if abs(closest1 - pos1) > resolution:
            #     # if verbose: print "no breakpoint within resolution:        pos1:",pos1,"   closest1:",closest1
            #     count_far_away +=1 

            chrom2 = fields[3]
            pos2 = (int(fields[4]) + int(fields[5]))/2
            closest2 = self.binary_search_for_closest(pos2,graph_breakpoints[chrom2],0,len(graph_breakpoints[chrom2])-1)
            strand2 = fields[9]
            # if abs(closest2 - pos2) > resolution:
            #     # if verbose: print "no breakpoint within resolution:        pos2:",pos2,"   closest2:",closest2
            #     count_far_away +=1 

            variant_name = fields[6]
            if variant_name in variant_names_used_already:
                print "Variant names are not unique.", variant_name, "appeared twice"
                return
            variant_names_used_already.add(variant_name)

            weight_split = float(fields[11])

            # Find ports and connect them with an edge
            key1 = (chrom1,closest1,strand1)
            key2 = (chrom2,closest2,strand2)

            self.create_split_edge(key1, key2, weight_split=weight_split,split_variant_name = variant_name)
            variant_counter += 1
            # print "successfully created", variant_counter, "split edges"

        f.close()

    def connect_nodes_with_span_edges(self,breakpoints):

        for chrom in breakpoints:
            pos_list = breakpoints[chrom]
            for pos in pos_list:
                key1 = (chrom,pos,"+")
                node1 = self.node_lookup.get(key1, "NA")
                port1 = "stop"
                if node1 == "NA":
                    print "Edge cannot be created because no node has this start or end position as a port:"
                    print key1
                    raise FileFormatError("trying to create edge at position that is not a node start or end port")
                rev_key1 = (key1[0],key1[1],reverse(key1[2]))
                rev_node1 = self.node_lookup.get(rev_key1, "NA")
                rev_port1 = reverse(port1)

                weight_span = min([self.nodes[node1].weight, self.nodes[rev_node1].weight])
                self.create_edge(node1,port1,rev_node1,rev_port1,weight_span,spansplit="span")


    def create_graph(self, variants_filename, coverage_file, genome_file, snap_within_bin_resolution = 2, CN_difference_threshold_to_split=None, verbose=False):

        # Find bin size from coverage_file and set resolution from it
        bin_size = int(self.bin_size(coverage_file))
        print "bin_size:", bin_size

        resolution = bin_size*snap_within_bin_resolution

        ########################### Get all breakpoints from copy number and variants #################################
        # Record breakpoints of all CNVs 
        breakpoints_from_copy_number = self.segmented_coverage_to_CNV_calls(coverage_file,cov_diff_threshold_to_split=CN_difference_threshold_to_split,csv=True)

        # Record breakpoints of all split read variants
        all_breakpoints_from_variants = self.breakpoints_from_variants(variants_filename,csv=True)

        ########################### Snap breakpoints together #################################
        consolidated_breakpoints = self.snap_breakpoints(breakpoints_from_copy_number, all_breakpoints_from_variants, resolution,verbose=verbose)
        
        ########################### Create nodes #################################
        nodes = self.nodes_from_breakpoints(consolidated_breakpoints,genome_file,csv=True)
        self.create_nodes_with_attributes(nodes)


        # Connect the nodes with split read variants at the closest breakpoint
        self.connect_nodes_at_variants(variants_filename,resolution,verbose=verbose)


        # Create all spanning edges
        self.connect_nodes_with_span_edges(consolidated_breakpoints)


        if verbose:
            print "Nodes:", len(self.nodes)
            print "Edges:"
            print self.count_span_vs_split_edges()

    


