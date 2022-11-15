# -*- coding: utf-8 -*-
from DGAS_FMPV2 import FMP_algo
from czh_findpath import CSGH_algo
from build_graph import init_graph, read_node_file
import pandas as pd
import datetime
import math
import multiprocessing
import sys
import os
import copy

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
    print(len(gr.nodes))
    print(len(gr.edges))
    g1 = 'G:HGNC:4236'
    g2 = 'G:HGNC:30249'
    g3 = 'G:HGNC:9756'
    czh = CSGH_algo(init_graph())
    #num, path=czh.find_meta_path(g4)
    candidate_meta_path=['GMDMDMG', 'GMDMG', 'G', 'GMG']
    for path in candidate_meta_path:
        nei1=czh.find_neighbor(g1,path)
        print("-------")
        print(path)
        print(len(nei1))
        print(nei1)

    meta_path_candidate = []








