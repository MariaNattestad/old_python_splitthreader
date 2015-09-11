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

    def test_graph_edges_have_weights(self):
        for edge in self.g.edges:
            self.assertEqual(edge.weight, 0)
        for edge in self.g2.edges:
            self.assertNotEqual(edge.weight, 0)

    def test_edge_weights_changed_everywhere(self):
        # (("A","start"),("B","stop")),
        self.g.nodes["A"].ports["start"].edges[0].weight = 100
        self.assertEqual(self.g.nodes["A"].ports["start"].edges[0].weight,100)
        self.assertEqual(self.g.nodes["B"].ports["stop"].edges[0].weight,100)
        self.assertEqual(self.g.edges[0].weight,100)

    def test_glide_along_edge(self):
        # (("A","start"),("B","stop"))
        this_port = self.g.nodes["A"].ports["start"]
        other_port = self.g.nodes["B"].ports["stop"]
        self.assertEqual(  this_port.edges[0].glide(this_port)  ,  other_port  )
        self.assertEqual(  other_port.edges[0].glide(other_port)  ,  this_port  )

        #(("B","start"),("C","start"))
        this_port = self.g.nodes["B"].ports["start"]
        other_port = self.g.nodes["C"].ports["start"]
        self.assertEqual(  this_port.edges[0].glide(this_port)  ,  other_port  )
        self.assertEqual(  other_port.edges[0].glide(other_port)  ,  this_port  )

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
        self.assertEqual(  str(this_port), "A:stop" )
        
        # travel one edge
        this_port = self.g.nodes["A"].ports["stop"]
        other_port = self.g.nodes["B"].ports["start"]
        self.assertEqual(  self.g.depth_first_search(this_port,other_port),  [["A:start","B:start"]])
        self.assertEqual(  str(this_port), "A:stop" )

        # travel two edges
        this_port = self.g.nodes["A"].ports["stop"]
        other_port = self.g.nodes["C"].ports["stop"]
        self.assertEqual(  self.g.depth_first_search(this_port,other_port),  [["A:start","B:start","C:stop"]])
        self.assertEqual(  str(this_port), "A:stop" )

        # no connection
        this_port = self.g.nodes["A"].ports["stop"]
        other_port = self.g.nodes["C"].ports["start"]
        self.assertEqual(  self.g.depth_first_search(this_port,other_port), None)
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
        self.assertEqual(  self.g.breadth_first_search(this_port,other_port), None)
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
        self.assertEqual(  self.g.breadth_first_search(this_port,other_port,depth_limit=0), None)
        self.assertEqual(  str(this_port), "A:stop" )

        # travel 2 edges, depth limit 1: NOT found
        this_port = self.g.nodes["A"].ports["stop"]
        other_port = self.g.nodes["C"].ports["stop"]
        self.assertEqual(  self.g.breadth_first_search(this_port,other_port,depth_limit=1), None)
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


def main():
    unittest.main()
    
    
if __name__ == '__main__':
    main()
    