import networkx as nx
class Graph:
    def __init__(self,node_list,edge_list):
        self.node_list=node_list
        self.edge_list=edge_list
        self.init_node(node_list)
        self.init_edge(edge_list)

    def read_node_file(self,path):
        nodes=[]
        with open(path,"r",encoding='UTF-8-sig') as fp:
            lines=fp.readlines()
            for line in lines:
                line=line[:-1]
                nodes.append(line)
        return nodes

    def test_find_pair(self):

        nodes = self.read_node_file("zyd_network/node/node_gene.csv")
        for node in nodes:
            for node2 in nodes:
                if self.getGraph().has_edge(node, node2):
                    print(node + " and " + node2 + " connected")

    def read_edge_file(self,path):
        edges=[]
        with open(path,"r",encoding='UTF-8-sig') as fp:
            lines=fp.readlines()
            for line in lines:
                pairs=line.split(',')
                edge=(pairs[0],pairs[1][:-1])
                edges.append(edge)
        return edges

    def init_node(self,node_list):
        self.nodes=[]
        for node_path in node_list:
            node=self.read_node_file(node_path)
            self.nodes=self.nodes+node

    def init_edge(self,edge_list):
        self.edges=[]
        for edge_path in edge_list:
            edge=self.read_edge_file(edge_path)
            self.edges=self.edges+edge

    def getGraph(self):
        graph=nx.Graph()
        #slice=random.sample(self.edges,186469)
        graph.add_nodes_from(self.nodes)
        graph.add_edges_from(self.edges)
        return graph


def init_graph():
    node_root = "zyd_network/node/"
    edge_root = "zyd_network/edge/"
    node_paths = [node_root + "node_disease.csv",
                  node_root + "node_drugs.csv",
                  node_root + "node_gene.csv",
                  node_root + "node_go.csv",
                  node_root + "node_miRNA.csv",
                  node_root + "node_phenotype(text).csv",
                  node_root + "node_protein.csv", ]
    edge_paths = [edge_root + "edge_disease2phenotype(text).csv",
                  edge_root + "edge_drugs2protein.csv",
                  edge_root + "edge_gene2disease.csv",
                  edge_root + "edge_gene2protein.csv",
                  edge_root + "edge_go2gene.csv",
                  edge_root + "edge_miRNA2disease.csv",
                  edge_root + "edge_miRNA2gene.csv",
                  ]
    g = Graph(node_paths, edge_paths)
    #nodes=g.read_node_file("zyd_network/node/node_gene.csv")
    my_graph = g.getGraph()

    return my_graph

def read_node_file(path):
    nodes = []
    with open(path, "r",encoding='UTF-8-sig') as fp:
        lines = fp.readlines()
        for line in lines:
            line = line[:-1]
            nodes.append(line)
    return nodes
def get_node_type(node):
    """
    Get node type
    :param node:
    :return: type
    """
    type=node[0]
    return type
