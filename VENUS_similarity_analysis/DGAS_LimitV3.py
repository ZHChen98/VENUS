import networkx as nx
import math
import datetime
import copy

class DGAS_limit_core:

    def __init__(self,limits,graph):
        self.map={}
        self.map1={}
        self.map2={}
        self.limit=[]
        for i in limits:
            flag = True
            for k in limits:
                if i==k:
                    continue
                elif i in k:
                    flag=False
            if flag:
                self.limit.append(i)
        print(self.limit)
        self.graph=graph

    def start(self,gene1,gene2,max_length):
        """
        :param gene1:
        :param gene2:
        :param max_length: max length of meta path.
        :return:score_list
                like [
                {
                'mPath': 'GMDMG',
                'score': 0.4,
                'ins': 3,   instance of meta Path
                'gene1': 'G:HGNC:17256',
                'gene2': 'G:HGNC:21497'
                }
                ]
        """
        self.map = {}
        self.map1 = {}
        self.map2 = {}
        graph = self.graph
        max_step=max_length-1

        # if not self.is_connected(gene1,gene2,graph,max_step*2+1):
        #     print( "there is no meta paths between "+gene1+" and "+gene2)
        #     return []
        meta_paths, meta_paths_from_gene1, meta_path_from_gene2=self.find_meta_paths(gene1,gene2,graph,max_step)
        if len(meta_paths)==0:
            return []
        meta_paths_list=[]
        for i in meta_paths:
            meta_paths_list.append(i["meta_path_name"])
        meta_paths_combine_list=self.combine_whole(meta_paths_list)

        score_list=self.score(meta_paths,meta_paths_from_gene1,meta_path_from_gene2,meta_paths_combine_list)
        for i in score_list:
            i['gene1']=gene1
            i['gene2']=gene2
        return score_list

    def combine_whole(self,l):
        lenth=len(l)
        list=[]
        if lenth==0:
            return []
        else:
            for i in range(1,lenth+1):
                list=list+self.combine(l,i)
        return list

    def combine(self,l, n):
        answers = []
        one = [0] * n

        def next_c(li=0, ni=0):
            if ni == n:
                answers.append(copy.copy(one))
                return
            for lj in range(li, len(l)):
                one[ni] = l[lj]
                next_c(lj + 1, ni + 1)

        next_c()
        return answers

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

        # print(time2-time1)
        # print(time3-time2)

        return mPaths,hmPaths,tmPaths

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

                        if len(i)>=3:
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
                    "half":i,
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
            type2 = self.get_node_type(neighbors)
            m_path_name2 = mPathName + type2
            # use pruning strategy
            if self.pruning(m_path_name2):
                continue
            new_path_name_detail=path_name_detail+"|"+neighbors
            self.save_map1(m_path_name2, neighbors,new_path_name_detail)
            if self.pruning2(m_path_name2):
                continue
            self.find_meta_paths_g1_center(neighbors, graph, k, m_path_name2, new_path_name_detail,old_node=gene)

    def find_meta_paths_g2_center(self, gene, graph, k, mPathName, path_name_detail,old_node=None):
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
            type2 = self.get_node_type(neighbors)
            m_path_name2 = mPathName + type2
            # use pruning strategy
            if self.pruning(m_path_name2):
                continue
            new_path_name_detail = path_name_detail + "|" + neighbors
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
        elif not any(mPathName in k for k in self.limit):
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

        if any(mPathName == k for k in self.limit):

            return True
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

    def score(self,mPaths,mPathDict1,mPathDict2,meta_paths_combine_list):
        """

        :param mPaths: meta paths between gene1 and gene2
        :param mPathDict1: meta paths between gene1 and center nodes
        :param mPathDict2: meta paths between gene2 and center nodes
        :return: score_list:
        """
        score_list=[]
        for i in meta_paths_combine_list:
            total_name=str(i)
            list_ins = []
            list_hmIns = []
            list_tmIns = []
            for name in i:
                for a in mPaths:
                    mPath=a['meta_path_name']
                    if mPath==name:
                        ins=a['ins']
                        list_ins.append(ins)
                        meta_path_half=a['half']
                        hmIns=self.get_ins_of_meta_path(meta_path_half,mPathDict1)
                        list_hmIns.append(hmIns)
                        tmIns=self.get_ins_of_meta_path(meta_path_half,mPathDict2)
                        list_tmIns.append(tmIns)
            score=self.score_kernal_function(list_ins,list_hmIns,list_tmIns)
            score_list.append({
                'mPath':total_name,
                'score':score,
            })

        return score_list

    def score_kernal_function(self,ins,hmIns,tmIns):
        """
        The score function
        :param ins: instance of meta path
        :param hmIns: instance of meta path from head to center node
        :param tmIns: instance of meta path from tail to center node
        :return: score
        """
        ins_total=1
        hmIns_total=1
        tmIns_total=1
        for i in ins:
            ins_total*=i
        for j in hmIns:
            hmIns_total*=j
        for k in tmIns:
            tmIns_total*=k

        score=2*ins_total/(hmIns_total+tmIns_total)
        return score


