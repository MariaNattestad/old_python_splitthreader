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
            assert("glide(): current port not in edge")


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


class Graph(object):
    def __init__(self):
        self.nodes = {}
        self.edges = []

    def __str__(self):
        return ", ".join(self.nodes.keys())

    def __repr__(self):
        return "Graph: [" +", ".join(self.nodes.keys()) + "]"

    def print_nodes(self):
        for node in g.nodes:
            print node, g.nodes[node].attributes
    def print_edges(self):        
        for edge in g.edges:
            print edge, edge.weight

    def create_nodes_with_attributes(self,node_attributes):
        for node_name in node_attributes:
            if node_name in self.nodes.keys():
                print "Node already exists in graph, not overwriting"
            else:
                self.nodes[node_name] = Node(node_name,node_attributes[node_name])

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

    def depth_first_search(self,current_port,destination_port,path_so_far=[],all_paths=[]):
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
        current_port = current_port.jump()
        saveport = str(current_port)
        if str(current_port) == str(destination_port):
            return all_paths + [path_so_far + [saveport]]
        edges = current_port.edges
        for edge in edges:
            current_port = edge.glide(current_port)
            return self.depth_first_search(current_port,destination_port,path_so_far+[saveport],all_paths)

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
        current_port = current_port.jump()
        # if match, then return:
        if str(current_port) == str(destination_port):
            return [[str(destination_port)]]

        queue = [(current_port,[current_port])]
        while queue:
            (port, path) = queue.pop(0)
            if depth_limit != -1 and len(path) > depth_limit:
                return None
            edges = port.edges
            for edge in edges:
                # ignore if already in path
                # glide
                current_port = edge.glide(current_port)
                # jump
                current_port = current_port.jump()
                # if match, then return
                if str(current_port) == str(destination_port):
                    return [map(str,path) + [str(destination_port)]]
                # else append to queue
                else:
                    queue.append((current_port, path+[current_port]))

    def bfs(graph, start, goal):
        queue = [(start, [start])]
        while queue:
            (vertex, path) = queue.pop(0)
            for next in graph[vertex] - set(path):
                if next == goal:
                    yield path + [next]
                else:
                    queue.append((next, path + [next]))

    def draw(self,output_filename,call_neato=True,max_linewidth=10,max_x_pixels = 600, max_y_pixels = 600):
        f=open(output_filename,'w')
        f.write("graph structs\n")
        f.write('{')
        # f.write('label = ""')
        f.write("\n")
        f.write('\tgraph [rankdir = "LR" bgcolor = "white" style = "filled"];\n')
        f.write('\tnode [shape = "record" style = "filled"];\n')
        f.write('\tedge [label="Edge" penwidth=12 fontcolor="red"];\n')

        x_max = 0.0000001
        y_max = 0.0000001
        for node in self.nodes.values():
            if node.x > x_max:
                x_max = node.x
            if node.y > y_max:
                y_max = node.y
        for node in self.nodes.values():
            # f.write('\t\tstruct' + str(counter) + '[label = "' + node.name + ' |{<f1> start |<f2> end}" shape = record fillcolor = "white" ];\n')
            f.write('\t\t' + node.name + ' [label = "' + node.name + ' |{<f1> start |<f2> end}" shape = record fillcolor = "white" pos = "' + str(float(node.x)/x_max*max_x_pixels) + "," + str(float(node.y)/y_max*max_y_pixels) + '"];\n')
        
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
            f.write('\t' + edge.n1 + ':' + p1 + ' -- ' + edge.n2 + ":" + p2 + ' [label="' +  str(label)  + '", color="black", penwidth=' + str(float(linewidth)) + ', fontcolor=yellow];\n')

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






# g = Graph()
# g.create_nodes(["A","B","C"])
# # g.create_nodes_with_attributes({"A":{"x":0,"y":0}, "B":{"x":2,"y":0}, "C":{"x":1,"y":1}})
# g.create_edges([   (("A","start"),("B","stop")),    (("B","start"),("C","start"))   ])

# g.print_nodes()


# # # print g

# # # g.print_edges()


# g.draw("/Users/mnattest/Desktop/SplitThreader_testcases/test.dot",call_neato=True)




############### edges are saved in the graph but are objects that are accessible everywhere, and when changed from one node's perspective change everywhere ###############
# print g.edges
# print "before:",g.nodes["A"].ports["start"].edges[0].weight

# g.nodes["A"].ports["start"].edges[0].weight = 0
# print "after:",g.nodes["A"].ports["start"].edges[0].weight
###########################################################################################################################################################################









