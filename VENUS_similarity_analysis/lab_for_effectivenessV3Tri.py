# -*- coding: utf-8 -*-
from DGAS_Limit_Tri import DGAS_limit_core
from DGAS_FMPV2_Tri import FMP_algo
from build_graph import init_graph, read_node_file
import pandas as pd
import datetime
import math
import multiprocessing
import sys
import os
import copy


def mutil_prcoessing(dict_parameter):
    gene1 = dict_parameter['gene1']
    piece_gene_nodes = dict_parameter['gene_nodes']
    full_gene_nodes = dict_parameter['full_gene']
    max_length = dict_parameter['max_length']
    limit = dict_parameter['limit']
    graph = dict_parameter['graph']
    now_time = datetime.datetime.now().timestamp()
    algo = DGAS_limit_core(limit, graph)
    count = 0
    lens = len(piece_gene_nodes)
    final_len = lens
    print("totol nodes:"+final_len)
    score_list = []
    for gene2 in piece_gene_nodes:
        count += 1
        if count % 100 ==0:
            print(count)
        for gene3 in full_gene_nodes:
            # count += 1
            # if count % 500 == 0:
            #     time = datetime.datetime.now().timestamp() - now_time
            #     array_time = time / (count)
            #     left = array_time * (final_len - count) / (60)
            #     print("the estimated time of processing:" + str(os.getpid()) + " is:" + str(left)[0:5] + " min")
            if gene1 == gene2 or gene2 == gene3 or gene1 == gene3:
                continue
            score = algo.start(gene1, gene2, gene3, max_length)
            # if score.__len__() > 0:
            #     print(score)
            score_list = score_list + score
    print(" ############Processing:" + str(os.getpid()) + " is done ###################")
    return score_list

def get_tuple(gr):
    Tnode = []
    with open("zyd_network/node/node_miRNA.csv", "r", encoding='UTF-8-sig') as fp:
        lines = fp.readlines()
        for line in lines:
            line = line[:-1]
            Tnode.append(line)
    f_c = {}
    f_l={}
    count = 0
    for node in Tnode:
        f_n=[]
        for G_neighbor in gr.neighbors(node):
            if G_neighbor[0]!='G':
                continue
            f_n.append(G_neighbor)
            if f_l.__contains__(G_neighbor):
                f_l[G_neighbor] += 1
            else:
                f_l[G_neighbor] = 1
        print(len(f_n))
        for i in range(0,len(f_n)):
            for j in range(i+1, len(f_n)):
                for k in range(j+1, len(f_n)):
                    f_set=[]
                    f_set.append(f_n[i])
                    f_set.append(f_n[j])
                    f_set.append(f_n[k])
                    f_set.sort()
                    f_string = "|".join(f_set)
                    if f_c.__contains__(f_string):
                        f_c[G_neighbor] += 1
                    else:
                        f_c[G_neighbor] = 1
    print("-----")
def combine_whole(l):
    lenth = len(l)
    list = []
    if lenth == 0:
        return []
    else:
        for i in range(1, lenth + 1):
            list = list + combine(l, i)
    return list


def combine(l, n):
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


def create_multi_task_data(gene_nodes, cores, max_length, gene1, limit, num_of_processings, graphs):
    list_len = len(gene_nodes)
    cut_count = math.floor(list_len / cores)
    cut_list = []
    print("Split data(" + str(list_len) + ") into " + str(cores) + " set:")

    for i in range(0, cores - 1):
        print(str(i * cut_count) + "--" + str(i * cut_count + cut_count - 1))
        piece = gene_nodes[i * cut_count:i * cut_count + cut_count - 1]
        cut_list.append({
            'gene1': gene1,
            'gene_nodes': piece,
            'full_gene':gene_nodes,
            'max_length': max_length,
            'limit': limit,
            'graph': graphs1[i]
        })
    i = cores - 1
    final_piece = gene_nodes[i * cut_count:list_len - 1]
    print(str(i * cut_count) + "--" + str(list_len - 1))
    cut_list.append({
        'gene1': gene1,
        'gene_nodes': final_piece,
        'full_gene': gene_nodes,
        'max_length': max_length,
        'limit': limit,
        'graph': graphs[i]
    })
    return cut_list


def get_parameter():
    if len(sys.argv) < 4:
        print(" Use G:HGNC:6932 and G:HGNC:9236 as input")
        # exit()
        return None, None, None
        # gene1="G:HGNC:6932"
        # fileName="HGNC6932"
    else:
        gene1 = "G:HGNC:" + str(sys.argv[1])
        gene2 = "G:HGNC:" + str(sys.argv[2])
        gene3 = "G:HGNC:" + str(sys.argv[3])
        fileName = "HGNC" + str(sys.argv[1]) + "_HGNC" + str(sys.argv[2]) + "_HGNC" + str(sys.argv[3])

        return gene1, gene2, gene3, fileName


def cut_meta_path(mPath):
    """
    Cut meta_path string into two part.
    One is a meta_path from the started gene node to the center node.
    Another is from the ended gene node.
    For example, input mPath like 'GDPdG' and the output are 'GDP' and 'Gdp"
    :param mPath:
    :return: head2center,tail2center
    """
    le = len(mPath)
    if le % 2 == 0:
        print("This meta path:" + mPath + " is not even")
        return None
    else:
        # 'GDG' 'GDPG
        center = math.floor(le / 2)
        head2center = mPath[0:center + 1]
        tail2center = mPath[center:le][::-1]
        if head2center == tail2center:
            return head2center
        else:
            print("This meta path:" + mPath + " is not symmetrical")
            return None


def read_rank(pd, meta_path_candidate, gene1, gene2, gene3):
    rank_list = []
    for meta_path in meta_path_candidate:
        pd2 = pd[pd['mPath'] == meta_path]
        pd2 = pd2.sort_values('score', ascending=False)
        pd2['score_rank'] = pd2['score'].rank(method='dense', ascending=False)
        max_rank = pd2['score_rank'].max()
        pd3 = pd2[
            (pd2['gene1'] == gene1) &
            (pd2['gene2'] == gene2) &
            (pd2['gene3'] == gene3)
            ]
        if len(pd3) == 0:
            continue
        meta_path_rank = pd3.iloc[0]['score_rank']
        rank_list.append({
            "meta_path_rank": meta_path_rank,
            "meta_path": meta_path,
            "gene1": gene1,
            "gene2": gene2,
            "gene3": gene3,
            "total_rank": max_rank,
            "rank_percent": round(1 - meta_path_rank / max_rank, 3)
        })

    return rank_list


if __name__ == '__main__':


    start_time = datetime.datetime.now().timestamp()
    gr = init_graph()
    get_tuple(gr)
    print(len(gr.nodes))
    print(len(gr.edges))
    gene_pair = pd.read_csv("label_result816/postive_pair_50_2.csv")
    gene_pair2 = pd.read_csv("label_result816/negative_pair_100.csv")
    max_length = 4
    gene_pair1 = gene_pair.to_records(index=None)
    gene_pair2 = gene_pair2.to_records(index=None)
    gene_pair = []
    for i in gene_pair1:
        gene_pair.append(i)
    for i in gene_pair2:
        gene_pair.append(i)
    # print(gene_pair)
    total = len(gene_pair)
    re_list = []
    re_strs = []
    re_rank_percent = []
    cores = 18
    graphs1 = []
    for i in range(0, math.floor(cores / 2) + 1):
        graphs1.append(init_graph())
    graphs2 = []
    for i in range(0, math.floor(cores / 2) + 1):
        graphs2.append(init_graph())
    cores = 16
    pool = multiprocessing.Pool(processes=cores)
    count_pair = 0

    time1 = datetime.datetime.now().timestamp()
    g1 = 'G:HGNC:4236'
    g2 = 'G:HGNC:30249'
    g3 = 'G:HGNC:9756'
    zyd = FMP_algo(init_graph())
    re = zyd.start(g1, g2, g3, 4)
    meta_path_candidate = []
    for r in re:
        meta_path_candidate.append(r["meta_path_name"])
    print(str(len(meta_path_candidate)) + " meta paths were found")

    meta_path_limit = []
    meta_path_chosen = []
    for candidate in meta_path_candidate:
        meta_path = cut_meta_path(candidate)
        if not meta_path == None:
            meta_path_chosen.append(candidate)
            print(candidate + "---" + meta_path)
            meta_path_limit.append(meta_path)
    print("The number of cpu cores is " + str(cores))
    gene_nodes = read_node_file("zyd_network/node/node_gene.csv")

    cores_for_gene1 = math.floor(cores / 3)
    cores_for_gene2 = math.floor(cores / 3)
    cores_for_gene3 = cores - cores_for_gene1 - cores_for_gene2
    multi_prcoessing_data1 = create_multi_task_data(gene_nodes, cores_for_gene1, 4, g1, meta_path_limit, cores,
                                                    graphs1)
    multi_prcoessing_data2 = create_multi_task_data(gene_nodes, cores_for_gene2, 4, g2, meta_path_limit, cores,
                                                    graphs2)
    multi_prcoessing_data3 = create_multi_task_data(gene_nodes, cores_for_gene3, 4, g3, meta_path_limit, cores,
                                                    graphs2)
    multi_prcoessing_data = multi_prcoessing_data1 + multi_prcoessing_data2 + multi_prcoessing_data3

    final_score_list = []

    for y in pool.imap_unordered(mutil_prcoessing, multi_prcoessing_data):
        final_score_list = final_score_list + y

    score_pd = pd.DataFrame(final_score_list)

    meta_path_candidate = combine_whole(meta_path_candidate)
    meta_path_candidate2 = []
    for i in meta_path_candidate:
        meta_path_candidate2.append(str(i))

    print(meta_path_candidate2)
    rank_list = read_rank(score_pd, meta_path_candidate2, g1, g2, g3)
    re_list = re_list + rank_list
    rank_pd = pd.DataFrame(rank_list)
    print(rank_pd)
    for it in rank_pd:
        print(it)
    rank_pd = rank_pd.sort_values('meta_path_rank')

    first_one = rank_pd.iloc[0]
    print("The result is:")
    print("Meta path:" + first_one['meta_path'] + ", rank is " + str(
        int(first_one['meta_path_rank'])) + " and total rank is " + str(int(first_one['total_rank'])))
    re_strs.append({
        "result": "Meta path:" + first_one['meta_path'] + ", rank is " + str(
            int(first_one['meta_path_rank'])) + " and total rank is " + str(int(first_one['total_rank'])),
        "gene1": g1,
        "gene2": g2
    })
    # rank_percent
    time2 = datetime.datetime.now().timestamp()
    print("Cost time:")
    cost_time = (time2 - time1) / 60
    print(str(cost_time)[0:4])
    rank_pd.to_csv("lab_result816/lab_effectiveness_150.csv", mode='a')
    rank_pd = rank_pd.sort_values('rank_percent', ascending=False)
    first_one = rank_pd.iloc[0]
    re_rank_percent.append({
        "rank_percent": first_one['rank_percent'],
        "gene1": g1,
        "gene2": g2,
        "meta_path": first_one['meta_path'],
    })
    print(first_one['rank_percent'])

    re_pd = pd.DataFrame(re_list)
    re_pd.to_csv("lab_result816/lab_result_effectiveness_150.csv")
    re_strs = pd.DataFrame(re_strs)
    re_strs.to_csv("lab_result816/lab_result_str_150.csv")
    re_rank_percent_pd = pd.DataFrame(re_rank_percent)
    re_rank_percent_pd.to_csv("lab_result816/lab_result_percent_150.csv", mode="a")
    end_time = datetime.datetime.now().timestamp()
    during_time = (end_time - start_time) / 60
    print("during time is:")
    print(during_time)






