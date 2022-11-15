# -*- coding: utf-8 -*-
# @Time : 2020/12/24 13:26
# @Author : Jclian91
# @File : util.py
# @Place : Yangpu, Shanghai

# 数据相关的配置
event_type = "label"

train_file_path = "./data/%s.train" % event_type
test_file_path = "./data/%s.test" % event_type
valid_file_path = "./data/%s.valid" % event_type

# 模型相关的配置
BERT_MODEL_DIR = "./PubMedBERT-base-uncased-abstract-fulltext"
MAX_LEN = 256
TRAIN_BATCH_SIZE = 32
VALID_BATCH_SIZE = 16
EPOCHS = 100
LEARNING_RATE = 2e-05
