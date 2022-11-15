import networkx as nx
import math

class FMP_algo:
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
        max_step=max_length-1

        if not self.is_connected(gene1,gene2,graph,max_step*2+1):
            print( "there is no meta paths between "+gene1+" and "+gene2)
            return []
        meta_paths=self.find_meta_paths(gene1,gene2,graph,max_step)
        if len(meta_paths)==0:
            return []
        return meta_paths

    def is_connected(self,gene1,gene2,graph,max_step):
        """
        Check if the two gene node is connected or not.
        Find a shortest path between two nodes.
        Return False if there is no path or the length of the path is
        longer than the max step.
        :param gene1:
        :param gene2:
        :param graph:
        :param max_step:
        :return:True,False
        """

        try:
            sp = nx.astar_path_length(graph, gene1, gene2)
        except:
            return False

        if sp>max_step:
            return False
        else:
            return True

    def find_meta_paths(self,gene1,gene2,graph,max_step):
        """
        Find all meta paths between two nodes within max_step.
        Meanwhile record all the meta paths which connect each nodes
        to the center node, in order to compute the score after.
        :param gene1:
        :param gene2:
        :param graph:
        :param max_step:
        :return:mPaths: meta paths between gene1 and gene2
                hmPaths: meta paths between gene1 and center nodes
                tmPaths: meta paths between gene2 and center nodes
        """
        center=max_step

        self.find_meta_paths_g1_center(gene1, graph, center, self.get_node_type(gene1),gene1)
        hmPaths = self.map1
        self.find_meta_paths_g2_center(gene2, graph, center, self.get_node_type(gene2),gene2)
        tmPaths = self.map2
        mPaths=self.match_meta_path(hmPaths,tmPaths)

        return mPaths

    def match_meta_path(self,head_map,tail_map):
        meta_path_list=[]
        for i,k in head_map.items():
            #print(i)
            if tail_map.__contains__(i):
                meta_path_name=i+i[0:-1][::-1]
                meta_count=0
                found=tail_map[i]
                for node,count in k.items():
                    if found.__contains__(node):
                        count2=found[node]['ins']
                        meta_count+=count['ins']*count2

                        if len(i)>3:
                            check_len=len(i)-3
                            list1 = count['list']
                            list2 = found[node]['list']
                            check_nodes=[]
                            check_nodes2=[]
                            for path1 in list1:
                                nodes=path1.split("|")
                                check_nodes.append(nodes[1:1+check_len])

                            for path2 in list2:
                                nodes2=path2.split('|')
                                check_nodes2.append(nodes2[1:1+check_len])
                            for c1 in check_nodes:
                                for c2 in check_nodes2:
                                    if len(set(c1).intersection(set(c2)))>0:
                                        meta_count-=1



                                        #{'mPath': 'GMDpDMG', 'score': 0.8, 'ins': 96, 'gene1': 'G:HGNC:9884','gene2': 'G:HGNC:24086'}
                                        #{'mPath': 'GMDpDMG', 'score': 0.5333333333333333, 'ins': 64, 'gene1': 'G:HGNC:9884', 'gene2': 'G:HGNC:24086'}]

                if meta_count ==0:
                    continue
                meta_path_list.append({
                    "meta_path_name":meta_path_name,
                    "ins":meta_count,

                })
        return meta_path_list



    def find_meta_paths_g1_center(self, gene, graph, k, mPathName,path_name_detail, old_node=None ):
        """
        Find all meta paths between two nodes within max_step.
        The loop stops only when k = 0.
        And the paths it find restore in self.map1.
        :param gene:
        :param graph:
        :param k:
        :param mPathName:
        :param old_node:
        :return: None
        """
        if k==0:
            return 0
        else:
            k-=1

        for neighbors in graph.neighbors(gene):
            # avoid creating circle in the meta path
            if neighbors == old_node:
                continue
            new_path_name_detail=path_name_detail+"|"+neighbors
            type2=self.get_node_type(neighbors)
            m_path_name2=mPathName+type2
            # use pruning strategy
            if self.pruning(m_path_name2):
                continue
            self.save_map1(m_path_name2, neighbors,new_path_name_detail)
            if self.pruning2(m_path_name2):
                continue
            self.find_meta_paths_g1_center(neighbors, graph, k, m_path_name2, new_path_name_detail,old_node=gene)

    def find_meta_paths_g2_center(self, gene, graph, k, mPathName, path_name_detail, old_node=None):
        """
        Find all meta paths between two nodes within max_step.
        The loop stops only when k = 0.
        And the paths it find restore in self.map2.
        :param gene:
        :param graph:
        :param k:
        :param mPathName:
        :param old_node:
        :return: None
        """
        if k == 0:
            return 0
        else:
            k -= 1
        for neighbors in graph.neighbors(gene):
            # avoid creating circle in the meta path
            if neighbors == old_node:
                continue
            new_path_name_detail = path_name_detail + "|" + neighbors
            type2 = self.get_node_type(neighbors)
            m_path_name2 = mPathName + type2
            # use pruning strategy
            if self.pruning(m_path_name2):
                continue
            self.save_map2(m_path_name2, neighbors, new_path_name_detail)
            if self.pruning2(m_path_name2):
                continue
            self.find_meta_paths_g2_center(neighbors, graph, k, m_path_name2,new_path_name_detail,old_node=gene)

    def save_map1(self, m_path_name2, neighbors,path_name_detail):
        if self.map1.__contains__(m_path_name2):
            if self.map1[m_path_name2].__contains__(neighbors):
                self.map1[m_path_name2][neighbors]['ins'] += 1
                self.map1[m_path_name2][neighbors]['list'].append(path_name_detail)

            else:
                self.map1[m_path_name2][neighbors]={
                'ins': 1,
                'list': [path_name_detail]
            }
        else:
            neig_dict = {
                'ins': 1,
                'list': [path_name_detail]
            }
            self.map1[m_path_name2] = {
                neighbors:neig_dict
            }

    def save_map2(self, m_path_name2, neighbors, path_name_detail):
        if self.map2.__contains__(m_path_name2):
            if self.map2[m_path_name2].__contains__(neighbors):
                self.map2[m_path_name2][neighbors]['ins'] += 1
                self.map2[m_path_name2][neighbors]['list'].append(path_name_detail)

            else:
                self.map2[m_path_name2][neighbors] = {
                    'ins': 1,
                    'list': [path_name_detail]
                }
        else:
            neig_dict = {
                'ins': 1,
                'list': [path_name_detail]
            }
            self.map2[m_path_name2] = {
                neighbors: neig_dict
            }



    def pruning(self,mPathName):
        """
        pruning strategy
        use the pruning strategy to reduce the searching space
        :param mPathName: String
        :return: True or False
        """
        if mPathName[-1]=="G":
            return True
        elif mPathName=="GMDM":
            return True
        elif mPathName=="GPdP":
            return True
        else:
            return False

    def pruning2(self,mPathName):
        """
        pruning strategy
        use the pruning strategy to reduce the searching space
        :param mPathName: String
        :return: True or False
        """
        if mPathName[-1]=="G":
            return True
        # elif mPathName=="GMD":
        #     return True
        # elif mPathName=="GPd":
        #     return True
        else:
            return False
    def get_node_type(self,node):
        """
        Get node type
        :param node:
        :return: type
        """
        type=node[0]
        return type

    def get_ins_of_meta_path(self,mPath,mPathDict):
        meta_path=mPathDict[mPath]
        ins=0
        if meta_path.__len__() ==0:
            return 0
        for i,k in meta_path.items():
            ins+=(k['ins']*k['ins'])
        return ins

    def cut_meta_path(self,mPath):
        """
        Cut meta_path string into two part.
        One is a meta_path from the started gene node to the center node.
        Another is from the ended gene node.
        For example, input mPath like 'GDPdG' and the output are 'GDP' and 'Gdp"
        :param mPath:
        :return: head2center,tail2center
        """
        le=len(mPath)
        if le%2==0:
            # 'GDDG'
            center=int(le/2)
            head2center=mPath[0:center]
            tail2center=mPath[center:le][::-1]
            return head2center,tail2center
        else:
            #'GDG' 'GDPG
            center = math.floor(le / 2)
            head2center=mPath[0:center+1]
            tail2center=mPath[center:le][::-1]
            return head2center,tail2center




