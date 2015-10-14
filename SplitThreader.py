#! /usr/bin/env python

#######################################################################################################################################
################################################        SplitThreader.py         ######################################################
#######################################################################################################################################

# Node has a name, attributes like chromosome, pos_start, pos_stop, and 
# Node has 2 ports
    # Each port has an unrestricted number of edges to other ports, owned by other nodes
    # Each port has the name of its parent node, but does not own the parent node

separator = "_____________________________________________________"

# Only necessary for draw function:
neato = "/usr/local/bin/neato" 
import random


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
        return self.attributes.get("x",random.random())

    @property
    def y(self):
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
        self.edges = [] # keep it a list instead of a dictionary so have can have multiple edges between the same ports

    def add_edge(self,edge_instance):
        self.edges.append(edge_instance)

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
        raise Exception("reverse only takes +, -, start, end")

class Graph(object):
    FileFormatError=FileFormatError

    def __init__(self):
        self.nodes = {} # key: node_name; value: Node object
        self.edges = [] # list of Edge objects
        self.node_lookup = {} # key: (chrom, pos, strand) tuple; #value: node_name

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

    def create_edges(self,edge_list,weight_list = []):
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
            self.edges.append(e)
            self.nodes[node1].ports[port1].edges.append(e)
            self.nodes[node2].ports[port2].edges.append(e)

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
        self.create_edge(node1,port1,node2,port2,weight_split)

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
        
        self.create_edge(node1,port1,rev_node1,rev_port1,weight_span1)

        # span 2
        rev_key2 = (key2[0],key2[1],reverse(key2[2]))
        rev_node2 = self.node_lookup.get(rev_key2, "NA")
        rev_port2 = reverse(port2)

        # print "Span 1:"
        # print key2, rev_key2
        # print node2, port2
        # print rev_node2, rev_port2
        
        self.create_edge(node2,port2,rev_node2,rev_port2,weight_span2)



    def create_edge(self,node1,port1,node2,port2,weight):
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
        self.edges.append(e)
        self.nodes[node1].ports[port1].edges.append(e)
        self.nodes[node2].ports[port2].edges.append(e)


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


    def depth_first_search_recurse(self,current_port,destination_port,all_paths,path_so_far=[],stop_when_found = False):
        # saving the ports after jumping, so if the path contains A:start, it means you went through A in the reverse direction. A:stop means forward direction. 

        ############# Basic steps: ################
        # jump
        # if current_port == destination_port:
        #     return all_paths + [path_so_far]
        # find edges of port
        # for edge in edges:
            # glide
            # recurse
        ###########################################
        
        jumped_port = current_port.jump()
        
        saveport = str(jumped_port)
        if str(jumped_port) == str(destination_port):
            # print "MATCH"
            all_paths.append(path_so_far + [saveport]) # new
        else:
            edges = jumped_port.edges
            for edge in edges:
                glide_port = edge.glide(jumped_port)
                if stop_when_found and len(all_paths)>0:
                    return
                else: # keep recursing
                    self.depth_first_search_recurse(glide_port, destination_port, all_paths, path_so_far+[saveport],stop_when_found=stop_when_found) # new

    def depth_first_search(self,current_port,destination_port,stop_when_found = False):
        all_paths = []
        self.depth_first_search_recurse(current_port=current_port,destination_port=destination_port,all_paths=all_paths,stop_when_found=stop_when_found)
        return all_paths

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
            for edge in edges:
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

    def draw(self,output_filename,call_neato=True,max_linewidth=10,max_x_pixels = 600, max_y_pixels = 600):
        f=open(output_filename,'w')
        f.write("graph structs\n")
        f.write('{')
        # f.write('label = ""')
        f.write("\n")
        f.write('\tgraph [rankdir = "LR" bgcolor = "white" style = "filled"];\n')
        f.write('\tnode [shape = "record" style = "filled"];\n')
        f.write('\tedge [label="Edge" penwidth=12 fontcolor="red"];\n')

        # x_max = 0.0000001
        # y_max = 0.0000001
        # for node in self.nodes.values():
        #     if node.x > x_max:
        #         x_max = node.x
        #     if node.y > y_max:
        #         y_max = node.y
        # print x_max
        # print y_max

        #####################
        x_maxes = {}
        x_mins = {}
        y_max = 0.0000001
        y_min = 9999999999999
        for node in self.nodes.values():
            if node.y > y_max:
                y_max = node.y
            if node.y < y_min:
                y_min = node.y
            if node.x > x_maxes.get(node.y,0):
                x_maxes[node.y] = node.x
            if node.x < x_mins.get(node.y,99999999999):
                x_mins[node.y] = node.x


        #####################
        for node in self.nodes.values():
            # f.write('\t\tstruct' + str(counter) + '[label = "' + node.name + ' |{<f1> start |<f2> end}" shape = record fillcolor = "white" ];\n')
            f.write('\t\t' + node.name + ' [label = "' + node.name + ' |{<f1> start |<f2> end}" shape = record fillcolor = "white" pos = "' + str((float(node.x)-x_mins[node.y])/(x_maxes[node.y]-x_mins[node.y])*max_x_pixels) + "," + str(float(node.y)/y_max*max_y_pixels) + '"];\n')
        
        max_edgeweight = 0.0000001
        for edge in self.edges:
            if edge.weight > max_edgeweight:
                max_edgeweight = edge.weight

        for edge in self.edges:
            print edge
            p1 = 'f1' if edge.p1 == 'start' else 'f2'
            p2 = 'f1' if edge.p2 == 'start' else 'f2'
            label = edge.weight
            linewidth = (edge.weight/max_edgeweight) * max_linewidth
            f.write('\t' + edge.n1 + ':' + p1 + ' -- ' + edge.n2 + ":" + p2 + ' [label="' +  str(label)  + '", color="black", penwidth=' + str(float(linewidth)) + ', fontcolor=gray];\n')

        f.write("}")
        f.close()

        if call_neato == True:
            import subprocess
            text_to_run = neato + " -n2 -Tpng " + output_filename + " > " + output_filename + ".png" # use flag -n2 when position is set
            print text_to_run
            #subprocess.call([text_to_run],shell=True) # for this to work, the output file must be f.close()d first 
            subprocess.check_output([text_to_run],shell=True) # for this to work, the output file must be f.close()d first 
            text_to_run = "open " + output_filename + ".png"
            subprocess.call([text_to_run],shell=True)

    def parsimony(self,use_breadth_first_search=False):
        portal_name = "Portal"
        self.add_portal(portal_name=portal_name)
        self.find_longest_path(use_breadth_first_search=use_breadth_first_search,portal_name=portal_name)

    def add_portal(self,portal_name="Portal"):
        # Set up portals at the unconnected nodes
        portals = set()
        for node_name in self.nodes:
            no_start = False
            no_stop = False
            if self.nodes[node_name].ports["start"].edges == []:
                no_start = True
            if self.nodes[node_name].ports["stop"].edges == []:
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


        
    def find_longest_path(self,use_breadth_first_search=False,portal_name="Portal"):
        # Two methods for finding all the paths:

        allpaths = []
        if use_breadth_first_search:
            allpaths = self.breadth_first_search(self.nodes[portal_name].ports["start"],self.nodes[portal_name].ports["start"])
        else:
            allpaths = self.depth_first_search(self.nodes[portal_name].ports["start"],self.nodes[portal_name].ports["start"])

                        # Checking that DFS and BFS produce the same results:

                        # print len(bfs)
                        # print len(dfs)
                        # bfs_lengths = set()
                        # for item in bfs:
                        #     bfs_lengths.add(len(item))
                        # dfs_lengths = set()
                        # for item in dfs:
                        #     dfs_lengths.add(len(item))
                        # print bfs_lengths-dfs_lengths
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
                    if current_uninterrupted_length > longest_uninterrupted_length_so_far:
                        longest_uninterrupted_length_so_far = current_uninterrupted_length
                        longest_uninterrupted_path_so_far = path
                    # Reset chromosome and length
                    current_uninterrupted_length = seq_length
                    current_chromosome = this_chromosome
                position_where_we_left_off = self.nodes[node].attributes[port] # port refers to after the jump, so that reflects the exit port out of this node
        return longest_uninterrupted_path_so_far,longest_uninterrupted_length_so_far




############### edges are saved in the graph but are objects that are accessible everywhere, and when changed from one node's perspective change everywhere ###############
# print g.edges
# print "before:",g.nodes["A"].ports["start"].edges[0].weight

# g.nodes["A"].ports["start"].edges[0].weight = 0
# print "after:",g.nodes["A"].ports["start"].edges[0].weight
###########################################################################################################################################################################
