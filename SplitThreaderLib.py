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


    def create_3_edges(self,key1,key2,weight_split,weight_span1,weight_span2):
        
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
            node_attributes[node_name] = {"chrom":fields[0],"start":int(fields[1]),"stop":int(fields[2]),"x":int(fields[1]),"y":y}
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

    def edges_from_path(self,path):
        second_port = self.port_from_path_item(path[0])
        edge_list = []
        for item in path[1:]:
            first_port = self.opposite_port_from_path_item(item)
            edge_list.append(second_port.edges[first_port])
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
            raise FileFormatError("More than one node contains the point")
        if len(matching_nodes) == 0:
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

    def path_sequence_length(self,path):
        distance = 0
        for item in path:
            intermediate_node = self.port_from_path_item(item).node
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
            # print "Intermediate nodes:", distance
            # print "First node:", node1.attributes["stop"]-pos1
            # print "Last node:", pos2 - node2.attributes["start"]
            distance += abs(node1.attributes[reverse(node1_start_port)]-pos1)
            distance += abs(pos2 - node2.attributes[reverse(node2_end_port)])
            distances.append(distance)
            
            # print node1.attributes["stop"]
            # total length of intermediate nodes + partial lengths of first and last nodes
        return allpaths,distances


    def read_annotation(self,annotation_file,remove_chr=True,name_field=4,strand_field=6):
        f = open(annotation_file)
        # c = 0
        for line in f:
            fields = line.strip().split()
            chrom = fields[0]
            if remove_chr:
                chrom = chrom.split("chr")[1]
            start = float(fields[1])
            stop = float(fields[2])
            name = fields[name_field-1]
            strand = fields[strand_field-1]
            self.annotation[name] = {"chrom":chrom,"start":start,"stop":stop,"strand":strand}
            # c+=1
            # if c > 1000:
            #     break
        f.close()

    def gene_fusion_distance(self,gene_name1,gene_name2,depth_limit,verbose=False):
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
                            for i in xrange(len(distances)):
                                reports.append({"Gene1_direction":direction_gene1, "Gene2_direction":direction_gene2, "path":paths[i],"distance":distances[i]})
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

    def gene_fusion_report(self,gene_name1,gene_name2,depth_limit=15,verbose = False):
        if verbose:
            print gene_name1,"-",gene_name2
        # print gene_name1, self.annotation.get(gene_name1)
        # print gene_name2, self.annotation.get(gene_name2)
        reports = self.gene_fusion_distance(gene_name1,gene_name2,depth_limit=depth_limit,verbose=False)

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
        # maybe put in a safety so genes that are already close to each other aren't reported as fusions unless the variants bring them closer so it's not read-through transcription


    def karyotype_from_parsimony(self,recordings,output_filename):
        record_counter = 1
        f_karyotype=open(output_filename + ".karyotype",'w')
        f_records = open(output_filename + ".parsimony",'w')
        for record in recordings:
            path = record[0]
            previous_chromosome = "0"
            chrom_length = 0
            f_records.write("%d" % (record_counter))
            for item in path:
                node,port = item.split(":")
                this_chromosome = self.nodes[node].attributes["chrom"]
                if node != "Portal":
                    f_records.write("\t%s:%d-%d;%s" % (this_chromosome,self.nodes[node].attributes["start"],self.nodes[node].attributes["stop"], "Forward" if port == "stop" else "Reverse"))
                    # f_records.write("\t%s" % (item))
                node_length = self.nodes[node].attributes["stop"]-self.nodes[node].attributes["start"]
                if this_chromosome == previous_chromosome:
                    chrom_length += node_length
                else:
                    if previous_chromosome != "0":
                        f_karyotype.write("%d\t%s\t%d\t%.2f\n" % (record_counter,previous_chromosome,chrom_length,record[2]))
                        # f_karyotype.write("path_"+str(record_counter)+"\t"+previous_chromosome + "\t" + str(chrom_length) + "\n")
                    previous_chromosome = this_chromosome
                    chrom_length = node_length
            f_records.write("\n")
            record_counter += 1

        f_karyotype.close()
        f_records.close()

        
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

    def find_path_lengths(self,allpaths):
        intact_lengths = []

        for path in allpaths:
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
                if this_chromosome == current_chromosome and this_position == position_where_we_left_off: 
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
            intact_lengths.append(longest_uninterrupted_length_so_far)
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

    def collapse_redundant_paths(self,paths):
        paths_used = []
        edge_lists_used = []
        for path in paths:
            edge_list = set(self.edges_from_path(path))
            if paths_used == []:
                paths_used.append(path)
                edge_lists_used.append(edge_list)
            else:
                redundant = False
                for previous_edge_list in edge_lists_used:
                    if edge_list == previous_edge_list:
                        redundant = True
                if redundant == False:
                    paths_used.append(path)
                    edge_lists_used.append(edge_list)

        return paths_used


    def parsimony(self,use_breadth_first_search=False,portal_name="Portal",verbose=True,depth_limit=20,cycle_limit=2,min_weight_required = 5):
        import time
        self.add_portal()
        # recording = []

        before = time.time()
        allpaths = self.depth_first_search(self.nodes[portal_name].ports["start"],self.nodes[portal_name].ports["start"],cycle_limit=cycle_limit,depth_limit=depth_limit)
        if verbose==True:
            print "DFS:  %.2f seconds" % (time.time()-before)
            print "Number of paths", len(allpaths)
        
        before = time.time()
        lengths = self.find_path_lengths(allpaths)
        if verbose==True: print "Intact length finding:  %.2f seconds" % (time.time()-before)

        import numpy as np
        lengths = np.array(lengths)
        indices = np.argsort(lengths)[::-1]
        before = time.time()
        
        if verbose==True: print "Sorting lengths:  %.2f seconds" % (time.time()-before)

        recordings = []
        for i in xrange(len(indices)): 
            index = indices[i]
            path = allpaths[index]
            weight = self.min_weight(path)
            if weight < min_weight_required:
                continue
            else:
                recordings.append([path,lengths[index],weight])
                self.subtract(path=path,weight=weight)

        return recordings

    def franken_path(self,path,output_filename = None):
        print "Franken-path:"
        print path
        counter = 1
        output = ""
        for item in path:
            node_name,port = item.split(":")
            node = self.nodes[node_name]
            strand = "+"
            if port == "start":
                strand = "-"
            output += "%s\t%d\t%d\t%03d\t0\t%s\n" % (node.attributes["chrom"],node.attributes["start"], node.attributes["stop"],counter, strand) 
            # output += "%03d\t%s:%d-%d\t%s\n" % (counter, node.attributes["chrom"],node.attributes["start"], node.attributes["stop"], strand)
            counter += 1
        if output_filename == None:
            print output
        else:
            f=open(output_filename,"w")
            f.write(output)
            f.close()
            print "bedtools getfasta -s -fi $REFERENCE -bed %s -fo %s.fasta" % (output_filename,output_filename)
        # print "bedtools getfasta -s -fi $REFERENCE -bed $FILE -fo $FILE.fasta"











############### edges are saved in the graph but are objects that are accessible everywhere, and when changed from one node's perspective change everywhere ###############
# print g.edges
# print "before:",g.nodes["A"].ports["start"].edges[0].weight

# g.nodes["A"].ports["start"].edges[0].weight = 0
# print "after:",g.nodes["A"].ports["start"].edges[0].weight
###########################################################################################################################################################################
