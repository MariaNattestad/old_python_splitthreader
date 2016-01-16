#! /usr/bin/env python

#######################################################################################################################################
################################################        SplitThreader.py         ######################################################
#######################################################################################################################################

# Node has a name, attributes like chromosome, pos_start, pos_stop, and 
# Node has 2 ports
    # Each port has an unrestricted number of edges to other ports, owned by other nodes
    # Each port has the name of its parent node, but does not own the parent node

class FileFormatError(Exception):
    pass

class Node(object):
    def __init__(self,name,attributes = {}):

        self.name = name
        self.ports = {}
        self.attributes = attributes

        # make two ports always
        self.add_port("start")
        self.add_port("stop")

    @property
    def x(self):
        import random
        return self.attributes.get("x",random.random())

    @property
    def y(self):
        import random
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
            print edge, edge.weight

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
            self.edges.append(e)
            # self.nodes[node1].ports[port1].edges.append(e)
            # self.nodes[node2].ports[port2].edges.append(e)
            self.nodes[node1].ports[port1].edges[self.nodes[node2].ports[port2]] = e
            self.nodes[node2].ports[port2].edges[self.nodes[node1].ports[port1]] = e


    def create_3_edges(self,key1,key2,weight_split,weight_span1,weight_span2,verbose=False):
        # if key1==key2:
        #     verbose = True


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
        self.create_edge(node1,port1,node2,port2,weight_split,spansplit="split")
        if verbose==True:
            print "nodes:", node1, node2
            
        # print "=================================================================="
        # print "Split:"
        # print key1,key2
        # print node1, port1
        # print node2, port2

        # span 1
        rev_key1 = (key1[0],key1[1],reverse(key1[2]))
        rev_node1 = self.node_lookup.get(rev_key1, "NA")
        rev_port1 = reverse(port1)

        # print "Span 1:"
        # print key1, rev_key1
        # print node1, port1
        # print rev_node1, rev_port1
        
        self.create_edge(node1,port1,rev_node1,rev_port1,weight_span1,spansplit="span")

        
        if key1 != key2:
            # span 2
            rev_key2 = (key2[0],key2[1],reverse(key2[2]))
            rev_node2 = self.node_lookup.get(rev_key2, "NA")
            rev_port2 = reverse(port2)

            # print "Span 1:"
            # print key2, rev_key2
            # print node2, port2
            # print rev_node2, rev_port2
            
            self.create_edge(node2,port2,rev_node2,rev_port2,weight_span2,spansplit="span")



    def create_edge(self,node1,port1,node2,port2,weight,spansplit="split"):
        # print node1
        # print port1
        # print self.nodes[node1].ports[port1]
        # print "_______________"
        # print node2
        # print port2
        # print self.nodes[node2].ports[port2]
        # print "_______________"
        e = Edge(self.nodes[node1].ports[port1], self.nodes[node2].ports[port2])
        e.weight = weight
        e.spansplit=spansplit
        self.edges.append(e)
        self.nodes[node1].ports[port1].edges[self.nodes[node2].ports[port2]] = e
        self.nodes[node2].ports[port2].edges[self.nodes[node1].ports[port1]] = e


    def read_spansplit(self,nodes_filename,edges_filename):
        f=open(nodes_filename)

        # Set up chromosomes in order first
        chromosome_locations = {}
        for i in xrange(23):
            chromosome_locations[str(i)] = i
        i+=1
        chromosome_locations["X"] = i
        i+=1 
        chromosome_locations["Y"] = i
        i+=1
        chromosome_locations["MT"] = i
        i+=1

        # Read in nodes from the spansplit.nodes file
        node_attributes = {}
        for line in f:
            fields = line.strip().split()
            y = 0
            if fields[0] in chromosome_locations.keys(): # Add other chromosomes to the dictionary if they aren't already in there
                y=chromosome_locations[fields[0]]
            else:
                chromosome_locations[fields[0]] = i
                i += 1
            node_name = "%s.%s" % (fields[0],fields[3])
            node_attributes[node_name] = {"chrom":fields[0],"start":int(fields[1]),"stop":int(fields[2]),"x":int(fields[1]),"y":y,"weight":float(fields[4])}
        self.create_nodes_with_attributes(node_attributes)

        f.close()

        # Read in edges from the spansplit file
        # counter = 0 # TESTING

        f=open(edges_filename)
        for line in f:
            fields = line.strip().split()
            # print fields
            chrom1 = fields[0]
            chrom2 = fields[3]
            pos1 = int(fields[1])
            pos2 = int(fields[4])
            strand1 = fields[8]
            strand2 = fields[9]
            weight_split = float(fields[11])
            weight_span1 = float(fields[13])
            weight_span2 = float(fields[16])


            key1 = (chrom1,pos1,strand1)
            key2 = (chrom2,pos2,strand2)

            # print key1
            # print key2

            self.create_3_edges(key1, key2, weight_split, weight_span1, weight_span2)

            # counter += 1        # TESTING
            # if counter > 10:    # TESTING
            #     break           # TESTING
        f.close()

    def to_csv(self,output_prefix):
        f=open(output_prefix + ".nodes.csv",'w')
        for node_name in self.nodes:
            node = self.nodes[node_name]
            f.write("%s,%s,%d,%d\n" % (node.name, node.attributes["chrom"], node.attributes["start"], node.attributes["stop"]));
        f.close()

        f=open(output_prefix + ".edges.csv",'w')
        for edge in self.edges:
            p1 = edge.ports[0]
            p2 = edge.ports[1]
            f.write("%s,%s,%s,%s,%f\n" % (p1.node.name,p1.name,p2.node.name,p2.name,edge.weight));
            
        f.close()

    def to_json(self,output_filename):
        node_id_dictionary = {}
        counter = 0
        f=open(output_filename,'w')
        f.write('{\n"nodes":[\n')
        for node_name in self.nodes:
            no_start = False
            no_stop = False
            if len(self.nodes[node_name].ports["start"].edges.values()) == 0:
                no_start = True
            if len(self.nodes[node_name].ports["stop"].edges.values()) == 0:
                no_stop = True

            node = self.nodes[node_name]
            prefix = ",\n"
            if counter == 0:
                prefix = ""
            length = 10
            if "start" in node.attributes and "stop" in node.attributes:
                length = node.attributes["stop"]-node.attributes["start"]
            chrom = "1"
            if "chrom" in node.attributes:
                chrom = node.attributes["chrom"]
            
            node_id_dictionary[node_name+":start"] = counter
            counter += 1
            significance = "None"
            if no_start == True:
                significance = "chrom_start"
            f.write(prefix + '\t{"name":"%s:start","chrom":"%s","seqlength":%d,"significance":"%s"}' % (node.name, chrom, length, significance));


            node_id_dictionary[node_name+":stop"] = counter
            counter += 1
            significance = "None"
            if no_stop == True:
                significance = "chrom_end"
            f.write(",\n" + '\t{"name":"%s:stop","chrom":"%s","seqlength":%d,"significance":"%s"}' % (node.name, chrom, length, significance));

        counter = 0
        f.write('\n],\n"links":[\n')
        for edge in self.edges:
            p1 = edge.ports[0]
            p2 = edge.ports[1]
            counter += 1
            prefix = ",\n"
            if counter == 1:
                prefix = ""
            chrom = "None"
            if edge.spansplit=="span":
                chrom = p1.node.attributes["chrom"]

            f.write(prefix + '\t{"source":%d,"target":%d,"value":%f,"attribute":"%s","chrom":"%s"}' % (node_id_dictionary[str(p1)],  node_id_dictionary[str(p2)],edge.weight,edge.spansplit,chrom));
        prefix = ",\n"

        for node_name in self.nodes:
            node = self.nodes[node_name]
            chrom = "1"
            if "chrom" in node.attributes:
                chrom = node.attributes["chrom"]
            length = 10
            if "start" in node.attributes and "stop" in node.attributes:
                length = node.attributes["stop"]-node.attributes["start"]
            f.write(prefix + '\t{"source":%d,"target":%d,"value":%f,"attribute":"sequence","chrom":"%s","seqlength":%d}' % (node_id_dictionary[node_name+":start"],  node_id_dictionary[node_name+":stop"],500,chrom,length));
            
        f.write("\n]\n}\n")
        f.close()

    # def paths_to_json(self,paths,output_filename):
    #     f = open(output_filename,"w")
    #     f.write('{\n"paths":[')
    #     for path in paths:
    #         f.write(path)

    #     f.close()


    def add_portal(self,portal_name="Portal"):
        # Set up portals at the unconnected nodes
        portals = set()
        for node_name in self.nodes:
            no_start = False
            no_stop = False
            if len(self.nodes[node_name].ports["start"].edges.values()) == 0:
                no_start = True
            if len(self.nodes[node_name].ports["stop"].edges.values()) == 0:
                no_stop = True
            if no_start and no_stop:
                pass # Don't make a node into a portal if it has no edges
            elif no_start:
                portals.add(self.nodes[node_name].ports["start"])
            elif no_stop:
                portals.add(self.nodes[node_name].ports["stop"])

        self.create_nodes_with_attributes({portal_name:{"chrom":"0","start":0,"stop":0,"x":0,"y":0}})
        for dead_end_port in portals:
            self.create_edge(dead_end_port.node.name, dead_end_port.name, portal_name, "stop", float("inf"))
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

    def min_weight(self,path):
        edge_list = self.edges_from_path(path)
        edge_dict = {}
        for edge in edge_list:
            edge_dict[edge] = edge_dict.get(edge,0) + 1

        minweight = -1
        for edge in edge_dict:
            corrected_weight = edge.weight*1.0/edge_dict[edge]
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
            if intermediate_node.name in count_only_nodes_in_set:
                distance += intermediate_node.attributes["stop"] - intermediate_node.attributes["start"]
        return distance

    def find_paths_by_directions(self,node1,node2,pos1,pos2,strand1,strand2,depth_limit=5):
        node1_start_port = "start" if strand1=="+" else "stop"
        node2_end_port = "stop" if strand2=="+" else "start"

        allpaths = self.breadth_first_search(node1.ports[node1_start_port], node2.ports[node2_end_port], depth_limit = depth_limit, stop_when_found = False)
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
        print "sample gene names used:"
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
            if counter < 5:
                print gene_name
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

    def gene_fusion_distance(self,gene_name1,gene_name2,additional_info = None, depth_limit=20, verbose=False):
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
                        paths,distances = self.find_paths_by_directions(node1,node2,pos1,pos2,strand1,strand2,depth_limit=depth_limit)  
                        if len(paths)>0:
                            direction_gene1 = "Forward"
                            if (annot1["strand"] == "-" and gene_end1 == "stop") or (annot1["strand"]=="+" and gene_end1 == "start"):
                                # print "Forward through", gene_name1
                                pass
                            else:
                                # print "Reverse through", gene_name1
                                direction_gene1 = "Reverse"
                            direction_gene2 = "Forward"
                            if (annot2["strand"] == "-" and gene_end2 == "start") or (annot2["strand"]=="+" and gene_end2 == "stop"):
                                # print "Forward through",gene_name2
                                pass
                            else:
                                # print "Reverse through",gene_name2
                                direction_gene2 = "Reverse"
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

        ################# This is the old way before I started encoding spansplit as an edge property
        # num_splits = -1

        # current_chromosome = ""
        # position_where_we_left_off = 0
        # for item in path:
        #     node,port = item.split(":")
        #     this_chromosome = self.nodes[node].attributes["chrom"]
        #     this_position = self.nodes[node].attributes[reverse(port)] # port refers to after the jump, so we reverse that to get the entry point into this node
        #     if this_chromosome == current_chromosome and this_position == position_where_we_left_off: 
        #         pass # reference spanning, not a split
        #     else:
        #         num_splits += 1
        #         current_chromosome = this_chromosome
        #     position_where_we_left_off = self.nodes[node].attributes[port] # port refers to after the jump, so that reflects the exit port out of this node
        return num_splits

    def split_weights_in_path(self,path):
        split_weights = []
        edges = self.edges_from_path(path)
        for edge in edges:
            if edge.spansplit == "split":
                split_weights.append(edge.weight)

        ################# This is the old way before I started encoding spansplit as an edge property
        # num_splits = -1

        # current_chromosome = ""
        # position_where_we_left_off = 0
        # for item in path:
        #     node,port = item.split(":")
        #     this_chromosome = self.nodes[node].attributes["chrom"]
        #     this_position = self.nodes[node].attributes[reverse(port)] # port refers to after the jump, so we reverse that to get the entry point into this node
        #     if this_chromosome == current_chromosome and this_position == position_where_we_left_off: 
        #         pass # reference spanning, not a split
        #     else:
        #         num_splits += 1
        #         current_chromosome = this_chromosome
        #     position_where_we_left_off = self.nodes[node].attributes[port] # port refers to after the jump, so that reflects the exit port out of this node
        return split_weights

    def gene_fusion_report(self,gene_name1,gene_name2,additional_info = None, depth_limit=15,verbose = False):
        if verbose:
            print gene_name1,"-",gene_name2
        # print gene_name1, self.annotation.get(gene_name1)
        # print gene_name2, self.annotation.get(gene_name2)
        reports = self.gene_fusion_distance(gene_name1,gene_name2,additional_info = additional_info, depth_limit=depth_limit,verbose=False)

        if len(reports) == 0:
            if verbose:
                print "No gene fusion detected"
        # elif len(reports) == 1:
        #     report = reports[0]
        #     print report["Gene1_direction"],"-",report["Gene2_direction"]
        #     # print gene_name1,report["Gene1_direction"]
        #     # print gene_name2,report["Gene2_direction"]
        #     print report["distance"]/1000., "kb"
        #     print len(report["path"])-1, "translocations"
        #     print report["path"]
            return None
        else:
            if verbose:
                print len(reports), "gene fusion(s) detected"
            scores = []
            for report in reports:
                score = 0
                num_splits = self.count_splits_in_path(report["path"])
                if report["Gene1_direction"] == report["Gene2_direction"] and num_splits==1 and report["distance"]<100000:
                    score = 150
                elif report["Gene1_direction"] == report["Gene2_direction"] and num_splits==1 and report["distance"]<1000000:
                    score = 100
                elif num_splits==1 and report["distance"]<1000000:
                    score = 70
                elif report["Gene1_direction"] == report["Gene2_direction"] and num_splits==2 and report["distance"]<1000000:
                    score = 50
                elif num_splits==2 and report["distance"]<1000000:
                    score = 40
                else:
                    score = 20
                if report["distance"] < 1000000:
                    score = score - 20*report["distance"]*1.0/1000000.
                scores.append(score)
            import numpy as np
            scores = np.array(scores)
            if max(scores) <= 20:
                return None
            indices = np.argsort(scores)[::-1]
            to_return = None
            for index in indices:
                if scores[index] == max(scores):
                    report = reports[index]
                    report["number_of_splits"] = self.count_splits_in_path(report["path"])
                    if verbose:
                        print scores[index], report["Gene1_direction"],"-",report["Gene2_direction"],"|", report["distance"]/1000., "kb","|", self.count_splits_in_path(report["path"]), "translocation(s)", report["path"]
                    to_return = report
            if verbose:
                print '__________________________'

            return to_return
        # maybe put in a safety so genes that are already close to each other aren't reported as fusions unless the variants bring them closer: so it's not read-through transcription

    ################# Still works but no longer used by Parsimony or in any command-line program ############################
    def find_longest_path(self,use_breadth_first_search=False,depth_limit=50,portal_name="Portal",required_minimum_edge_weight=1,cycle_limit=2):
        # Two methods for finding all the paths:
        allpaths = []
        if use_breadth_first_search:
            allpaths = self.breadth_first_search(self.nodes[portal_name].ports["start"],self.nodes[portal_name].ports["start"],depth_limit=depth_limit)
        else:
            allpaths = self.depth_first_search(self.nodes[portal_name].ports["start"],self.nodes[portal_name].ports["start"],cycle_limit=cycle_limit,depth_limit=depth_limit)

        # print allpaths
        
        longest_uninterrupted_path_so_far = []
        longest_uninterrupted_length_so_far = 0
        
        for path in allpaths:
            current_uninterrupted_length = 0
            current_chromosome = ""
            position_where_we_left_off = 0
            for item in path:
                node,port = item.split(":")
                # print "Stop:",self.nodes[node].attributes["stop"]
                # print "Start:",self.nodes[node].attributes["start"]
                # # Attributes: # {"chrom":fields[0],"start":int(fields[1]),"stop":int(fields[2]),"x":int(fields[1]),"y":y}
                seq_length = self.nodes[node].attributes["stop"]-self.nodes[node].attributes["start"]
                this_chromosome = self.nodes[node].attributes["chrom"]
                this_position = self.nodes[node].attributes[reverse(port)] # port refers to after the jump, so we reverse that to get the entry point into this node
                if this_chromosome == current_chromosome and this_position == position_where_we_left_off: 
                    current_uninterrupted_length += seq_length
                else:
                    # Save if this path is the best so far
                    if current_uninterrupted_length > longest_uninterrupted_length_so_far and self.min_weight(path) >= required_minimum_edge_weight:
                        longest_uninterrupted_length_so_far = current_uninterrupted_length
                        longest_uninterrupted_path_so_far = path
                    # Reset chromosome and length
                    current_uninterrupted_length = seq_length
                    current_chromosome = this_chromosome
                position_where_we_left_off = self.nodes[node].attributes[port] # port refers to after the jump, so that reflects the exit port out of this node
        return longest_uninterrupted_path_so_far,longest_uninterrupted_length_so_far

    def find_intact_length(self,path,count_only_nodes_in_set=[]):
        longest_uninterrupted_length_so_far = 0
        current_uninterrupted_length = 0
        current_chromosome = ""
        position_where_we_left_off = 0
        for item in path:
            node,port = item.split(":")
            # print "Stop:",self.nodes[node].attributes["stop"]
            # print "Start:",self.nodes[node].attributes["start"]
            # # Attributes: # {"chrom":fields[0],"start":int(fields[1]),"stop":int(fields[2]),"x":int(fields[1]),"y":y}
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

    def subtract(self,path,weight):
        edges = self.edges_from_path(path)
        for edge in edges:
            edge.weight = edge.weight - weight

    def breadth_first_search(self, current_port, destination_port, depth_limit = -1, stop_when_found = False):
        
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
        # jump
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





    def cycle_depth_first_search_recurse(self,current_port,destination_port,cycles_found,depth_limit,path_so_far=[],depth=0):
        # saving the ports after jumping, so if the path contains A:start, it means you went through A in the reverse direction. A:stop means forward direction. 

        ############# Basic steps: ################
        # jump
        # if current_port == destination_port:
        #     return allpaths + [path_so_far]
        # find edges of port
        # for edge in edges:
            # glide
            # recurse
        ###########################################
        
        jumped_port = current_port.jump()
        
        saveport = str(jumped_port)
        if str(jumped_port) == str(destination_port) or depth > depth_limit :
            pass
            # print "MATCH"
            # allpaths.append(path_so_far + [saveport]) # new
        else:
            edges = jumped_port.edges.values()
            for edge in edges:
                glide_port = edge.glide(jumped_port)
                if str(jumped_port) in path_so_far:
                    cycles_found.append(path_so_far+[saveport])
                else: # keep recursing
                    self.cycle_depth_first_search_recurse(current_port=glide_port, destination_port=destination_port, cycles_found=cycles_found, path_so_far=path_so_far+[saveport],depth_limit=depth_limit,depth=depth+1) # new

    def find_cycles(self,depth_limit=20):
        portal_name="Portal"
        self.add_portal()

        cycles_found = []
        self.cycle_depth_first_search_recurse(current_port=self.nodes[portal_name].ports["start"], destination_port=self.nodes[portal_name].ports["start"], cycles_found=cycles_found,depth_limit=depth_limit)
        
        cycles_only = []
        for cycle in cycles_found:
            index = cycle.index(cycle[-1])
            cycles_only.append(cycle[index:])

        return self.collapse_redundant_paths(cycles_only)

    def group_redundant_gene_fusions(self,paths_dict_by_fusion_name,split_edges_only = False):
        paths_used = []
        edge_lists_used = []
        names_used = []

        groups = {}

        # i = 0
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
                        # print names_used[j], "<----", names[i], "\t(same path, possibly overlapping genes)"
                        # print "\t%s\t(same fusion path as %s)" % (names[i], names_used[j])
                        groups[names_used[j]].append(fusion_name)
                        break
                if redundant == False:
                    paths_used.append(path)
                    edge_lists_used.append(edge_list)
                    names_used.append(fusion_name)
                    groups[fusion_name] = [fusion_name]
            # i += 1

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

    def subgraph_from_genome_interval(self,chromosome,start,end,degree):
        nodes = self.nodes_within_genome_interval(chromosome,start,end)
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
        
        # if len(matching_nodes) == 0:
            # print "Warning: No node matches: ",chromosome, ":",start,"-",end
            # raise FileFormatError("No nodes match the point")


        # Returns node names
        return matching_nodes

    def first_degree_nodes(self,node_names):
        first_degree_nodes = node_names + []
        for node_name in node_names:
            node = self.nodes[node_name]
            for port_name in self.nodes[node_name].ports:
                edges = self.nodes[node_name].ports[port_name].edges
                for other_port in edges:
                    first_degree_nodes.append(other_port.node.name)

        return list(set(first_degree_nodes))


    def subgraph_from_nodes(self,node_names_given,degree_given=2): 
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
        
        # Add nodes to the graph
        s.create_nodes_with_attributes(node_attributes)

        # Add edges to the graph
        for edge_name in edges_to_add:
            edge = edges_to_add[edge_name]
            node1 = edge.ports[0].node.name
            port1 = edge.ports[0].name
            node2 = edge.ports[1].node.name
            port2 = edge.ports[1].name
            weight = edge.weight
            spansplit = edge.spansplit
            s.create_edge(node1,port1,node2,port2,weight,spansplit=spansplit)

        return s

    def parsimony(self,use_breadth_first_search=False,portal_name="Portal",verbose=False,depth_limit=20,cycle_limit=0,min_weight_required = 10,chop_end_nodes=0,count_only_nodes_in_set=[]):
        import time
        self.add_portal()
        # recording = []

        before = time.time()
        allpaths = self.depth_first_search(self.nodes[portal_name].ports["start"],self.nodes[portal_name].ports["start"],cycle_limit=cycle_limit,depth_limit=depth_limit)
        if verbose==True:
            print "DFS:  %.2f seconds" % (time.time()-before)
            print "Number of paths", len(allpaths)
        
        # Chop off the end nodes on each path for the subgraph/local region case as these are inflated and can differentiate paths that are actually the same within the specific region
        if chop_end_nodes != 0:
            newpaths = []
            for path in allpaths:
                if len(path) > chop_end_nodes*2:
                    newpaths.append(path[chop_end_nodes:-chop_end_nodes])
            allpaths = newpaths

        intact_lengths = self.find_path_intact_lengths(allpaths,count_only_nodes_in_set)
        total_lengths = self.find_path_total_lengths(allpaths,count_only_nodes_in_set)
        sorting_parameters = []

        # Sorting by 1) Longest intact path within region, 2) Longest total path within region, 3) Fewest splits
        for i in xrange(len(allpaths)):
            sorting_parameters.append((i,intact_lengths[i],total_lengths[i],-1*self.count_splits_in_path(allpaths[i])))
        
        import operator
        ordering = sorted(sorting_parameters,key=operator.itemgetter(1,2,3))[::-1]

        recordings = []
        # for i in xrange(len(indices)): 
        for item in ordering:
            index = item[0]
            # index = indices[i]
            path = allpaths[index]
            weight = self.min_weight(path)
            if weight < min_weight_required:
                continue
            else:
                # print path
                # print weight
                recordings.append([path,intact_lengths[index],weight])
                self.subtract(path=path,weight=weight)

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

    def local_parsimony(self,chromosome,start,end,output_prefix,degree=3,min_weight_required=10,depth_limit=30):
        self.calculate_edge_weights_from_node_coverages()




        s = self.subgraph_from_genome_interval(chromosome,start,end,degree=degree)


        # s.calculate_edge_weights_from_node_coverages()
        




        print "Created subgraph of this region with %d nodes and %d edges\n" % (len(s.nodes),len(s.edges))
        # s.to_json(output_prefix+".subgraph.json")
        region_nodes = set(s.nodes_within_genome_interval(chromosome,start,end))

        print region_nodes

        reports = s.parsimony(depth_limit = depth_limit,chop_end_nodes=degree,count_only_nodes_in_set=region_nodes,min_weight_required=min_weight_required) 

        # print "Beautiful art showing the paths in this region:"
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
                # print report
                output += "\n"
        print output
        filtered_reports = filtered_reports[::-1] # flip back around, to undo the ordering that is best for drawing, and return to the order that the paths were found in: longest intact sequence first

        print "Found a parsimonious set of", len(filtered_reports), "paths."
        
        # s.bed_files_from_parsimony(filtered_reports,output_filename=output_prefix)
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

    def bed_files_from_parsimony(self,recordings,output_filename):
        path_counter = 1

        # print "Header for bed file:\nchromosome\tstart\tend\tpath_ID\ttotal_length\tstrand\tminimum_read_depth_on_breakpoints\tlongest_intact_sequence_length\n"
        for record in recordings:
            f=open(output_filename+".path_%03d.bed"%(path_counter),'w')
            path = record[0][1:-1] # Cut off Portal nodes on either side of the path
            minweight = record[2]
            longest_intact_length = record[1]
            total_length = self.find_total_length(path)
            # previous_chromosome = "0"
            # chrom_length = 0
            # f_records.write("path_%03d" % (record_counter))
            for item in path:
                node_name,port = item.split(":")
                node = self.nodes[node_name]
                strand = "+"
                if port == "start":
                    strand = "-"
                f.write("%s\t%d\t%d\tpath_%03d\t%d\t%s\t%d\t%d\n" % (node.attributes["chrom"],node.attributes["start"], node.attributes["stop"],path_counter, total_length, strand,minweight,longest_intact_length))
            path_counter += 1   
            f.close()

    def karyotype_from_parsimony(self,recordings,output_filename):
        # record_counter = 1
        # f_karyotype=open(output_filename + ".karyotype.txt",'w')
        # f_records = open(output_filename + ".parsimony.txt",'w')
        allpaths = []
        path_counter = 1
        f_bed = open(output_filename + ".paths.bed",'w')
        f_bed.write("chromosome\tstart\tend\tpath_ID\ttotal_length\tstrand\tminimum_read_depth_on_breakpoints\tlongest_intact_sequence_length\n")
        for record in recordings:
            path = record[0][1:-1] # Cut off Portal nodes on either side of the path
            minweight = record[2]
            longest_intact_length = record[1]
            total_length = self.find_total_length(path)
            # previous_chromosome = "0"
            # chrom_length = 0
            # f_records.write("path_%03d" % (record_counter))
            for item in path:
                node_name,port = item.split(":")
                node = self.nodes[node_name]
                strand = "+"
                if port == "start":
                    strand = "-"
                f_bed.write("%s\t%d\t%d\tpath_%03d\t%d\t%s\t%d\t%d\n" % (node.attributes["chrom"],node.attributes["start"], node.attributes["stop"],path_counter, total_length, strand,minweight,longest_intact_length))
            path_counter += 1

            #     node,port = item.split(":")
            #     this_chromosome = self.nodes[node].attributes["chrom"]
            #     if node != "Portal":
            #         f_records.write("\t%s:%d-%d;%s" % (this_chromosome,self.nodes[node].attributes["start"],self.nodes[node].attributes["stop"], "Forward" if port == "stop" else "Reverse"))
            #         # f_records.write("\t%s" % (item))
            #     node_length = self.nodes[node].attributes["stop"]-self.nodes[node].attributes["start"]
            #     if this_chromosome == previous_chromosome:
            #         chrom_length += node_length
            #     else:
            #         if previous_chromosome != "0":
            #             f_karyotype.write("path_%03d\t%s\t%d\t%.2f\n" % (record_counter,previous_chromosome,chrom_length,record[2]))
            #             # f_karyotype.write("path_"+str(record_counter)+"\t"+previous_chromosome + "\t" + str(chrom_length) + "\n")
            #         previous_chromosome = this_chromosome
            #         chrom_length = node_length
            # f_records.write("\n")
            # record_counter += 1

        # f_karyotype.close()
        # f_records.close()
        f_bed.close()

    def calculate_edge_weights_from_node_coverages(self):
        print "//////////////////////////////////////////////////////"

        flagged_splits = []
        adjusted_splits = []

        counter = 0
        for edge in self.edges:

            if edge.spansplit == "split":
            
                # print edge
                # print edge.spansplit

                port1 = edge.ports[1]
                port2 = edge.ports[0]

                node1 = port1.node
                node2 = port2.node


                n1 = node1.attributes["weight"]
                n2 = node2.attributes["weight"]

                # print node1, "weight:", n1
                # print node2, "weight:", n2

                # print edge.ports[0].edges.keys()
                # print edge.ports[1].edges.keys()
                
                
                # print port1
                # print port2

                span_port1 = ""
                span_edge1 = ""

                for key in port1.edges.keys():
                    if key != port2:
                        if span_port1 == "":
                            span_port1 = port1.edges[key].glide(port1)
                            span_edge1 = port1.edges[key]
                        else:
                            # print "span_port1 already set. Old:", span_port1, "   New:", port1.edges[key].glide(key)
                            # print port1.edges
                            span_port1 = ""
                            span_edge1 = ""
                            break

                span_port2 = ""
                span_edge2 = ""

                for key in port2.edges.keys():
                    if key != port1:
                        if span_port2 == "":
                            span_port2 = port2.edges[key].glide(port2)
                            span_edge2 = port2.edges[key]
                        else:
                            # print "span_port2 already set. Old:", span_port2, "   New:", port2.edges[key].glide(key)
                            # print port2.edges
                            span_port2 = ""
                            span_edge2 = ""
                            break

                # print "Opposite ports:"
                # print span_port1
                # print span_port2


                if span_port1 != "" and span_port2 != "": # Only if the spanning nodes are present do we continue calculations
                    x1 = span_port1.node.attributes["weight"]
                    # print "span weight 1:", x1
                    
                    x2 = span_port2.node.attributes["weight"]
                    # print "span weight 2:", x2

                    est_split_weight_1 = n1 - x1
                    est_split_weight_2 = n2 - x2

                    # print "estimated weights:", est_split_weight_1, est_split_weight_2
                    percent_diff = 100*abs(est_split_weight_1-est_split_weight_2)/((est_split_weight_1+est_split_weight_2)/2)
                    # print "percentage difference:", percent_diff

                    if est_split_weight_1 < 0 or est_split_weight_2 < 0:
                        # print edge
                        # print node1, "weight:", n1
                        # print node2, "weight:", n2

                        # print span_port1
                        # print span_port2
                        # print "n1:", n1
                        # print "n2:", n2
                        # print "x1:", x1
                        # print "x2:", x2

                        flagged_splits.append(edge)
                        # if port1 == port2:
                        #     print "Inverted duplication with negatives"
                        #     print edge
                        #     print span_port1
                        #     print span_port2
                        #     print "n1:", n1
                        #     print "n2:", n2
                        #     print "x1:", x1
                        #     print "x2:", x2
                        #     print "est_split_weight_1:",est_split_weight_1
                        #     print "est_split_weight_2:",est_split_weight_2
                        #     print "-----------"

                    else: #if percent_diff < 100:
                        # Then it's reasonable and we will use the average of the est_split_weight_1/2 instead of the split read count
                        # old_edge_weight = edge.weight
                        edge.weight = (est_split_weight_1+est_split_weight_2)/2

                        # print edge, edge.weight-old_edge_weight
                        span_edge1.weight = x1
                        span_edge2.weight = x2
                        adjusted_splits.append(edge)

                        # if port1 == port2:
                        #     print "Inverted duplication"
                        #     print edge
                        #     print span_port1
                        #     print span_port2
                        #     print "n1:", n1
                        #     print "n2:", n2
                        #     print "x1:", x1
                        #     print "x2:", x2
                        #     print "est_split_weight_1:",est_split_weight_1
                        #     print "est_split_weight_2:",est_split_weight_2
                        #     print "-----------"



                        
                    # print "-----------"

                    # if counter > 6:
                    #     break
                    counter += 1
        print "total:",counter

        print "flagged:",len(flagged_splits)
        print "adjusted:",len(adjusted_splits)


        print "//////////////////////////////////////////////////////"

    def search_for_fusions(self,gene_pair_list_file,output_fusion_report_file):
        f=open(gene_pair_list_file)
        
        gene_fusion_reports = []
        gene_names = []
        for line in f:
            fields = line.strip().split()
            gene_names.append((fields[0],fields[1]))
            gene_fusion_reports.append(self.gene_fusion_report(fields[0],fields[1],fields[2:]))

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
                # f_output_fusion_report.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (gene_names[i][0],gene_names[i][1], "none","none","none","none","none"))
            else:
                total_supported_fusions += 1
                if report["Gene1_direction"]=="Reverse" and report["Gene2_direction"]=="Reverse":
                    report["Gene1_direction"]="Forward"
                    report["Gene2_direction"]="Forward"
                    tmp = report["Gene1"]
                    report["Gene1"] = report["Gene2"]
                    report["Gene2"] = tmp
                key = (report["Gene1"],report["Gene2"])
                paths_dict_by_fusion_name[key] = report["path"]
                reports_dict_by_fusion_name[key] = report
                
        f.close()
        print "Found evidence of %d gene fusions out of the list of %d given:" % (total_supported_fusions,total_putative_fusions)


        groups = self.group_redundant_gene_fusions(paths_dict_by_fusion_name, split_edges_only = True)
        final_reports = []
        num_split_RNA_reads = []
        redundant_fusions = {}
        for group in groups:
            gene_lengths = {}
            for fusion_name in group:
                # print reports_dict_by_fusion_name[fusion_name]
                # print reports_dict_by_fusion_name[fusion_name]["info"]
                annot1 = self.annotation[fusion_name[0]]
                # print annot1
                annot2 = self.annotation[fusion_name[1]]
                # print annot2
                length1 = annot1["stop"]-annot1["start"]
                length2 = annot2["stop"]-annot2["start"]
                gene_lengths[fusion_name] = length1 + length2
            # print gene_lengths
            max_fusion_name = ""
            max_gene_length = 0
            for fusion_name in gene_lengths:
                if gene_lengths[fusion_name] > max_gene_length:
                    max_fusion_name = fusion_name
                    max_gene_length = gene_lengths[fusion_name]
            # print "MAX:", max_fusion_name,max_gene_length
            redundant_fusions[max_fusion_name] = set(group)-set([max_fusion_name])

            final_reports.append(reports_dict_by_fusion_name[max_fusion_name])
            num_split_RNA_reads.append(int(reports_dict_by_fusion_name[max_fusion_name]["info"][0]))


        f_output_fusion_report = open(output_fusion_report_file,"w")
        f_output_fusion_report.write("gene_1\tgene2\tRNA_split_read_count\tDNA_split_reads_at_variants\ttotal_fusion_gene_length\tstrand_gene_1\tstrand_gene_2\tnumber_of_splits_to_thread_through\tminimum_spansplit_weights_along_path\tother_fusions_through_same_variant(s)\n")

        import numpy as np
        sorted_indices = np.argsort(num_split_RNA_reads)[::-1]
        for index in sorted_indices:
            report = final_reports[index]
            all_split_weights = map(str,self.split_weights_in_path(report["path"]))
            alternate_names = []
            for item in redundant_fusions[(report["Gene1"], report["Gene2"])]:
                alternate_names.append(item[0] + "=" + item[1])
            rna_split_read_count = int(report["info"][0])
            f_output_fusion_report.write("%s\t%s\t%d\t%s\t%d\t%s\t%s\t%d\t%f\t%s\n" % (report["Gene1"], report["Gene2"], rna_split_read_count, ",".join(all_split_weights),report["distance"],report["Gene1_direction"],report["Gene2_direction"],report["number_of_splits"], self.min_weight(report["path"]), ",".join(alternate_names)   ))

        f_output_fusion_report.close()

        # print "%d of these are through unique paths" % (len())

        # print "Wrote output to %s" % (output_fusion_report_file)



############### edges are saved in the graph but are objects that are accessible everywhere, and when changed from one node's perspective change everywhere ###############
# print g.edges
# print "before:",g.nodes["A"].ports["start"].edges[0].weight

# g.nodes["A"].ports["start"].edges[0].weight = 0
# print "after:",g.nodes["A"].ports["start"].edges[0].weight
###########################################################################################################################################################################
