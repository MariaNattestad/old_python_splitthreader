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


    def test_parsimony(self):
        self.assertEqual(len(self.g.parsimony(depth_limit=50)),11)


class TestParsimonyGenomeWide(unittest.TestCase):
    def test_parsimony_on_whole_genome(self):
        self.g = Graph()

        testcase_dir = "/Users/mnattest/Desktop/SplitThreader_testcases/"
        nodes_filename = testcase_dir + "bwamem.hg19.readname_sorted.mq60.bd200.mw10.primary_chr.over10kb.spansplit.nodes.bed"
        edges_filename = testcase_dir + "bwamem.hg19.readname_sorted.mq60.bd200.mw10.primary_chr.over10kb.spansplit.bedpe"
        
        self.g.read_spansplit(nodes_filename,edges_filename)

        # print "Number of nodes:", len(self.g.nodes)
        # print "Number of edges:", len(self.g.edges)

        recordings = self.g.parsimony(use_breadth_first_search=True,verbose=False,depth_limit=10)
        # print len(recordings)
        # 3.7-3.8 seconds per longest path subtraction (finding all paths once), when depth_limit is 20
        # when depth limit is 21, it takes 7.5 seconds per build
        # That means we probably have O(2^N) complexity

        # for item in recordings:
        #     print item
        #     print "_____________________________________________________________"

class TestAnnotations(unittest.TestCase):
    def setUp(self):
        self.g = Graph()
        testcase_dir = "/Users/mnattest/Desktop/SplitThreader_testcases/"
        nodes_filename = testcase_dir + "Her2.spansplit.nodes.bed"
        edges_filename = testcase_dir + "Her2.spansplit.bedpe"
        nodes_filename = testcase_dir + "bwamem.hg19.readname_sorted.mq60.bd200.mw10.primary_chr.over10kb.spansplit.nodes.bed"
        edges_filename = testcase_dir + "bwamem.hg19.readname_sorted.mq60.bd200.mw10.primary_chr.over10kb.spansplit.bedpe"
        
        self.g.read_spansplit(nodes_filename,edges_filename)

    def test_calculate_distance(self):
        point1 = ("17",37278042) # within node 17.9
        point2 = ("17",38066674) # within node 17.9
        self.assertEqual(self.g.find_nodename_by_position(point1),"17.9")
        path,distance = self.g.calculate_distance(point1,point2)
        self.assertEqual(distance,38066674-37278042) # 2 points in same node

        point1 = ("17",37278042) # within node 17.9
        point2 = ("17",38473392) # within node 17.10
        self.assertEqual(self.g.find_nodename_by_position(point1),"17.9")
        self.assertEqual(self.g.find_nodename_by_position(point2),"17.10")
        path,distance = self.g.calculate_distance(point1,point2)
        self.assertEqual(distance,38473392-37278042) # 2 points on neighboring nodes

        point1 = ("17",37278042) # within node 17.9
        point2 = ("17",39541958) # within node 17.10
        self.assertEqual(self.g.find_nodename_by_position(point1),"17.9")
        self.assertEqual(self.g.find_nodename_by_position(point2),"17.11")
        path,distance = self.g.calculate_distance(point1,point2,depth_limit=3) # with whole genome there is now a shortcut, so set the depth_limit to 3 
        self.assertEqual(distance,39541958-37278042) # 2 points on nodes with one node in the middle

        point1 = ("17",37278042) # within node 17.9
        point2 = ("17",41401528) # end of node 17.13
        self.assertEqual(self.g.find_nodename_by_position(point1),"17.9")
        self.assertEqual(self.g.find_nodename_by_position(point2),"17.13")
        path,distance = self.g.calculate_distance(point1,point2)
        self.assertEqual(distance,41401528-37278042) # 2 points on nodes with three nodes in the middle

        point1 = ("17",37278042) # within node 17.9
        point2 = ("8",121560937) # within node 8.165
        self.assertEqual(self.g.find_nodename_by_position(point1),"17.9")
        self.assertEqual(self.g.find_nodename_by_position(point2),"8.165")
        # travel along edge: 8.165:start--17.9:start 
        path,distance = self.g.calculate_distance(point1,point2)
        self.assertEqual(distance,11000)

    def test_read_annotation(self):
        testcase_dir = "/Users/mnattest/Desktop/SplitThreader_testcases/"
        annot_filename = testcase_dir + "gencode.v19.annotation.gtf.genes.bed"
        self.g.read_annotation(annot_filename,name_field=8)
        self.assertEqual(len(self.g.annotation),55765)
        
        self.assertEqual(self.g.count_splits_in_path(self.g.gene_fusion_report("CYTH1","EIF3H")["path"]),2)
        self.assertEqual(self.g.count_splits_in_path(self.g.gene_fusion_report("TATDN1","GSDMB")["path"]),1)
        self.assertEqual(self.g.count_splits_in_path(self.g.gene_fusion_report("SUMF1","LRRFIP2")["path"]),1)
        self.assertEqual(self.g.count_splits_in_path(self.g.gene_fusion_report("SUMF1","LRRFIP2")["path"]),1)
        self.assertEqual(self.g.count_splits_in_path(self.g.gene_fusion_report("PHF20","RP4-723E3.1")["path"]),1)
        self.assertEqual(self.g.count_splits_in_path(self.g.gene_fusion_report("RARA","PKIA")["path"]),1)
        self.assertEqual(self.g.count_splits_in_path(self.g.gene_fusion_report("MTBP","SAMD12")["path"]),1)
        self.assertEqual(self.g.count_splits_in_path(self.g.gene_fusion_report("CYTH1","MTBP")["path"]),1)
        self.assertEqual(self.g.count_splits_in_path(self.g.gene_fusion_report("TOX2","STAU1")["path"]),1)
        self.assertEqual(self.g.count_splits_in_path(self.g.gene_fusion_report("LINC00536","PVT1")["path"]),1)
        
        self.assertEqual(self.g.count_splits_in_path(self.g.gene_fusion_report("SNTB1","KLHDC2")["path"]),2) # Maybe only in Lumpy calls

        # self.g.gene_fusion_report("CYTH1","EIF3H") # good
        # self.g.gene_fusion_report("TATDN1","GSDMB") # good
        # self.g.gene_fusion_report("SUMF1","LRRFIP2") # good
        # self.g.gene_fusion_report("TAF2","COLEC10") # good
        # self.g.gene_fusion_report("TRIO","FBXL7") # good
        # self.g.gene_fusion_report("ATAD5","TLK2") # good
        # self.g.gene_fusion_report("DHX35","ITCH") # good
        # self.g.gene_fusion_report("RARA","PKIA") # good
        # self.g.gene_fusion_report("PHF20","RP4-723E3.1") # good

        # These go through multiple spanning edges but still have only 1 real split, so I implemented the count_splits_in_path() function to show that they have only 1 translocation
        # self.g.gene_fusion_report("MTBP","SAMD12") # good
        # self.g.gene_fusion_report("LINC00536","PVT1") # good
        # self.g.gene_fusion_report("TBC1D31","ZNF704") # good
        # self.g.gene_fusion_report("TOX2","STAU1") # good
        # self.g.gene_fusion_report("CYTH1","MTBP") # good
        # self.g.gene_fusion_report("SAMD12","EXT1") # good # Forward-Reverse

        ###########################
        # self.g.gene_fusion_report("PREX1","PHF20") # DNA can link up in two different ways. Which one is right can be seen from IsoSeq to see which ends of the genes are transcribed in the fusion. 
        # PHF20-PREX1 with 2 splits (Beginning of PHF20 into beginning of PREX1)
        # PREX1-PHF20 with 1 split (End of PREX1 into end of PHF20)
        ###########################

        # Should appear with Sniffles but not with Lumpy:
        # self.g.gene_fusion_report("AMZ2","CASC8") 

        # # Putative gene fusion through 2 translocations, possibly novel
        # self.g.gene_fusion_report("SNTB1","KLHDC2") # YAY

def main():
    unittest.main()



if __name__ == '__main__':
    main()


