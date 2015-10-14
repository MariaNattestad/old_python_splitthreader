#!/usr/bin/env python

import unittest
from SplitThreader import *




class TestNode(unittest.TestCase):
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.node = Node("A")
        self.nodeatts = Node("A",{"x":10,"y":50})
    def test_node_created(self):
        self.assertIsInstance(self.node,Node)
    def test_node_has_name(self):
        self.assertEqual(self.node.name, "A")
        self.assertEqual(self.nodeatts.name, "A")
    def test_node_has_attributes(self):
        self.assertEqual(self.node.attributes, {})
        self.assertEqual(self.nodeatts.attributes, {"x":10,"y":50})
    def test_node_has_ports(self):
        self.assertTrue(len(self.node.ports)==2)
        self.assertIsInstance(self.node.ports["start"],Port)
        self.assertIsInstance(self.node.ports["stop"],Port)

class TestGraph(unittest.TestCase):
    def setUp(self):
        self.g = Graph()
        self.g.create_nodes(["A","B","C"])
        self.g.create_edges([   (("A","start"),("B","stop")),    (("B","start"),("C","start"))  ])
        
        self.g2 = Graph()
        self.g2.create_nodes_with_attributes({"A":{"x":0,"y":0}, "B":{"x":2,"y":0}, "C":{"x":1,"y":1}})
        self.g2.create_edges([   (("A","start"),("B","stop")),    (("B","start"),("C","start"))   ], [40,100])

    def test_graph_created(self):
        self.assertIsInstance(self.g, Graph)
        self.assertIsInstance(self.g2, Graph)

    def test_graph_has_nodes(self):
        self.assertEqual(len(self.g.nodes), 3)
        self.assertEqual(len(self.g2.nodes), 3)

    def test_graph_nodes_have_attributes(self):
        for node in self.g.nodes.values():
            self.assertEqual(node.attributes, {})
        for node in self.g2.nodes.values():
            self.assertNotEqual(node.attributes, {})

    def test_graph_has_edges(self):
        self.assertEqual(len(self.g.edges),2)
        self.assertEqual(len(self.g2.edges),2)

    def test_port_has_edges(self):
        self.assertEqual(self.g.nodes["A"].ports["start"].edges[self.g.nodes["B"].ports["stop"]].glide(self.g.nodes["A"].ports["start"]),self.g.nodes["B"].ports["stop"])

    def test_graph_edges_have_weights(self):
        for edge in self.g.edges:
            self.assertEqual(edge.weight, 0)
        for edge in self.g2.edges:
            self.assertNotEqual(edge.weight, 0)

    def test_edge_weights_changed_everywhere(self):
        # (("A","start"),("B","stop")),
        first_port = self.g.nodes["A"].ports["start"]
        second_port = self.g.nodes["B"].ports["stop"]


        self.g.nodes["A"].ports["start"].edges[second_port].weight = 100
        self.assertEqual(self.g.nodes["A"].ports["start"].edges[second_port].weight,100)
        self.assertEqual(self.g.nodes["B"].ports["stop"].edges[first_port].weight,100)
        self.assertEqual(self.g.edges[0].weight,100)

    def test_glide_along_edge(self):
        # (("A","start"),("B","stop"))
        this_port = self.g.nodes["A"].ports["start"]
        other_port = self.g.nodes["B"].ports["stop"]
        self.assertEqual(  this_port.edges[other_port].glide(this_port)  ,  other_port  )
        self.assertEqual(  other_port.edges[this_port].glide(other_port)  ,  this_port  )

        #(("B","start"),("C","start"))
        this_port = self.g.nodes["B"].ports["start"]
        other_port = self.g.nodes["C"].ports["start"]
        self.assertEqual(  this_port.edges[other_port].glide(this_port)  ,  other_port  )
        self.assertEqual(  other_port.edges[this_port].glide(other_port)  ,  this_port  )

    def test_jump_across_node(self):
        this_port = self.g.nodes["A"].ports["start"]
        other_port = self.g.nodes["A"].ports["stop"]
        self.assertEqual(  this_port.jump()  ,  other_port  )
        self.assertEqual(  other_port.jump()  ,  this_port  )

    def test_depth_first_search(self):
        # only need to jump
        this_port = self.g.nodes["A"].ports["stop"]
        other_port = self.g.nodes["A"].ports["start"]
        self.assertEqual(  self.g.depth_first_search(this_port,other_port),  [["A:start"]])
        self.assertEqual(  self.g.depth_first_search(this_port,other_port,stop_when_found=True),  [["A:start"]])
        self.assertEqual(  str(this_port), "A:stop" )
        
        # travel one edge
        this_port = self.g.nodes["A"].ports["stop"]
        other_port = self.g.nodes["B"].ports["start"]
        self.assertEqual(  self.g.depth_first_search(this_port,other_port),  [["A:start","B:start"]])
        self.assertEqual(  self.g.depth_first_search(this_port,other_port,stop_when_found=True),  [["A:start","B:start"]])
        self.assertEqual(  str(this_port), "A:stop" )

        # travel two edges
        this_port = self.g.nodes["A"].ports["stop"]
        other_port = self.g.nodes["C"].ports["stop"]
        self.assertEqual(  self.g.depth_first_search(this_port,other_port),  [["A:start","B:start","C:stop"]])
        self.assertEqual(  self.g.depth_first_search(this_port,other_port,stop_when_found=True),  [["A:start","B:start","C:stop"]])
        self.assertEqual(  str(this_port), "A:stop" )

        # no connection
        this_port = self.g.nodes["A"].ports["stop"]
        other_port = self.g.nodes["C"].ports["start"]
        self.assertEqual(  self.g.depth_first_search(this_port,other_port), [])
        self.assertEqual(  self.g.depth_first_search(this_port,other_port,stop_when_found=True), [])
        self.assertEqual(  str(this_port), "A:stop" )


    def test_breadth_first_search(self):
        # only need to jump
        this_port = self.g.nodes["A"].ports["stop"]
        other_port = self.g.nodes["A"].ports["start"]
        self.assertEqual(  self.g.breadth_first_search(this_port,other_port),  [["A:start"]])
        self.assertEqual(  str(this_port), "A:stop" )
        
        # travel one edge
        this_port = self.g.nodes["A"].ports["stop"]
        other_port = self.g.nodes["B"].ports["start"]
        self.assertEqual(  self.g.breadth_first_search(this_port,other_port),  [["A:start","B:start"]])
        self.assertEqual(  str(this_port), "A:stop" )

        # travel two edges
        this_port = self.g.nodes["A"].ports["stop"]
        other_port = self.g.nodes["C"].ports["stop"]
        self.assertEqual(  self.g.breadth_first_search(this_port,other_port),  [["A:start","B:start","C:stop"]])
        self.assertEqual(  str(this_port), "A:stop" )

        # no connection
        this_port = self.g.nodes["A"].ports["stop"]
        other_port = self.g.nodes["C"].ports["start"]
        self.assertEqual(  self.g.breadth_first_search(this_port,other_port), [])
        self.assertEqual(  str(this_port), "A:stop" )
    
    def test_breadth_first_search_depth_limit(self):

        # jump, depth limit 0: FOUND
        this_port = self.g.nodes["A"].ports["stop"]
        other_port = self.g.nodes["A"].ports["start"]
        self.assertEqual(  self.g.breadth_first_search(this_port,other_port,depth_limit=0),  [["A:start"]])
        self.assertEqual(  str(this_port), "A:stop" )

        # travel one edge, depth limit 1: FOUND
        this_port = self.g.nodes["A"].ports["stop"]
        other_port = self.g.nodes["B"].ports["start"]
        self.assertEqual(  self.g.breadth_first_search(this_port,other_port,depth_limit=1),  [["A:start","B:start"]])
        self.assertEqual(  str(this_port), "A:stop" )

        # travel one edge, depth limit 0: NOT found
        this_port = self.g.nodes["A"].ports["stop"]
        other_port = self.g.nodes["B"].ports["start"]
        self.assertEqual(  self.g.breadth_first_search(this_port,other_port,depth_limit=0), [])
        self.assertEqual(  str(this_port), "A:stop" )

        # travel 2 edges, depth limit 1: NOT found
        this_port = self.g.nodes["A"].ports["stop"]
        other_port = self.g.nodes["C"].ports["stop"]
        self.assertEqual(  self.g.breadth_first_search(this_port,other_port,depth_limit=1), [])
        self.assertEqual(  str(this_port), "A:stop" )

        # travel 2 edges, depth limit 2: FOUND
        this_port = self.g.nodes["A"].ports["stop"]
        other_port = self.g.nodes["C"].ports["stop"]
        self.assertEqual(  self.g.breadth_first_search(this_port,other_port,depth_limit=2), [["A:start","B:start","C:stop"]])
        self.assertEqual(  str(this_port), "A:stop" )

        # travel 2 edges, depth limit -1 (means don't use depth limit): FOUND
        this_port = self.g.nodes["A"].ports["stop"]
        other_port = self.g.nodes["C"].ports["stop"]
        self.assertEqual(  self.g.breadth_first_search(this_port,other_port,depth_limit=-1), [["A:start","B:start","C:stop"]])
        self.assertEqual(  str(this_port), "A:stop" )

    def test_subtract_edge_twice(self):
        self.g3 = Graph()
        self.g3.create_nodes_with_attributes({"A":{"x":0,"y":0}, "B":{"x":2,"y":0}, "C":{"x":1,"y":1}})
        self.g3.create_edges([   (("A","start"),("B","stop"))  ,  (("B","start"),("C","start"))  ,  (("C","stop"),("C","stop"))   ], [40,100,80])

        this_port = self.g3.nodes["A"].ports["stop"]

        self.g3.print_edges()

        allpaths = self.g3.depth_first_search(this_port,this_port)
        path = allpaths[0]

        minweight = self.g3.min_weight(path)
        self.g3.subtract(path,minweight)
        self.assertEqual(self.g3.min_weight(path),0)


class TestReadingSpansplit(unittest.TestCase):

    def setUp(self):
        self.g = Graph()
    def test_read_spansplit(self):
        testcase_dir = "/Users/mnattest/Desktop/SplitThreader_testcases/"
        nodes_filename = testcase_dir + "Her2.spansplit.nodes.bed"
        edges_filename = testcase_dir + "Her2.spansplit.bedpe"
        
        # nodes_filename = testcase_dir + "bwamem.hg19.readname_sorted.mq60.bd200.mw10.primary_chr.over10kb.spansplit.nodes.bed"
        # edges_filename = testcase_dir + "bwamem.hg19.readname_sorted.mq60.bd200.mw10.primary_chr.over10kb.spansplit.bedpe"
        
        self.g.read_spansplit(nodes_filename,edges_filename)
        # self.g.print_edges()
        
        self.assertEqual(len(self.g.nodes),33)
        self.assertEqual(len(self.g.edges),15)
        

        ############## Draw function ###############
        # self.g.draw("/Users/mnattest/Desktop/SplitThreader_testcases/test.dot")
        


class TestParsimony(unittest.TestCase):
    def setUp(self):
        self.g = Graph()
        testcase_dir = "/Users/mnattest/Desktop/SplitThreader_testcases/"
        nodes_filename = testcase_dir + "Her2.spansplit.nodes.bed"
        edges_filename = testcase_dir + "Her2.spansplit.bedpe"
        self.g.read_spansplit(nodes_filename,edges_filename)

    def test_portal_and_that_BFS_and_DFS_agree(self):
        portal_name="Portal"
        self.g.add_portal(portal_name=portal_name)
        ## Checking that DFS and BFS produce the same results:
        dfs = self.g.depth_first_search(self.g.nodes[portal_name].ports["start"],self.g.nodes[portal_name].ports["start"])
        bfs = self.g.breadth_first_search(self.g.nodes[portal_name].ports["start"],self.g.nodes[portal_name].ports["start"])
        self.assertEqual(len(bfs),30)
        self.assertEqual(len(dfs),30)
        self.assertEqual(set(map(tuple,bfs)),set(map(tuple,dfs)))
        
    def test_find_longest_path(self):
        portal_name="Portal"
        self.g.add_portal()
        path,length = self.g.find_longest_path()

        self.assertEqual(len(path),8)
        self.assertEqual(length,4995158)

    def test_min_weight_and_subtract(self):
        portal_name="Portal"
        self.g.add_portal()
        path,length = self.g.find_longest_path()
        self.assertEqual(self.g.min_weight(path),14.0)

        # Subtract 13, check new minweight is 1
        self.g.subtract(path,13)
        path,length = self.g.find_longest_path()
        self.assertEqual(self.g.min_weight(path),1.0)

        # Subtract the last 1, check new longest path is different
        self.g.subtract(path,1)
        new_path,length = self.g.find_longest_path()
        self.assertNotEqual(set(map(tuple,path)),set(map(tuple,new_path)))


    # def test_parsimony(self):
        # self.g.parsimony()


def main():
    unittest.main()
    
    
if __name__ == '__main__':
    main()
    