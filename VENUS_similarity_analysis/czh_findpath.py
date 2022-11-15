import networkx as nx
import math

class CSGH_algo:
    """
    1.改变搜索策略
    """

    def __init__(self,graph):
        self.map={}
        self.map1={}
        self.map2={}
        self.graph=graph


    def start(self,gene1,gene2,max_length):
        self.ccc=0
        """
        :param gene1:
        :param gene2:
        :param max_length: max length of meta path.
        :return:meta_paths
        """
        self.map = {}
        self.map1 = {}
        self.map2 = {}

        graph = self.graph

    def is_connected(self,gene1,gene2,graph,max_step):

        try:
            sp = nx.astar_path_length(graph, gene1, gene2)
        except:
            return False

        if sp>max_step:
            return False
        else:
            return True

    def is_symmetry(self, str):
        for index in range(len(str)):
            if str[index] != str[len(str) - index - 1]:
                return False
        return True
    def find_meta_path(self,gene):
        visited=[]
        current_layer=[]
        pathway_item={}
        pathway=[]
        pathway_item[gene]="G"
        current_layer.append(gene)
        visited.append(gene)
        neighbor_num=0
        while len(current_layer)>0:
            candidate_nodes=current_layer
            print(current_layer)
            current_layer=[]
            for g in candidate_nodes:
                for neighbors in self.graph.neighbors(g):
                    if neighbors in visited:
                        continue
                    visited.append(neighbors)
                    pathway_item[neighbors]=pathway_item[g]+self.get_node_type(neighbors)
                    if self.get_node_type(neighbors) == "G":
                        neighbor_num = neighbor_num + 1
                        continue
                    current_layer.append(neighbors)
        for keys in pathway_item.keys():
            if(self.is_symmetry(pathway_item[keys])):
                pathway.append(pathway_item[keys])
        return neighbor_num, list(set(pathway))

    def find_neighbor(self,gene,path):
        visited = []
        current_layer = []
        current_layer.append(gene)
        visited.append(gene)
        layer=1
        while len(current_layer) > 0 and layer <len(path):
            candidate_nodes = current_layer
            print(current_layer)
            current_layer = []
            for g in candidate_nodes:
                for neighbors in self.graph.neighbors(g):
                    if neighbors in visited:
                        continue
                    visited.append(neighbors)
                    if self.get_node_type(neighbors) == path[layer]:
                        current_layer.append(neighbors)
            layer=layer+1
        return list(set(current_layer))


    def get_node_type(self,node):
        type=node[0]
        return type



