# -*- coding: utf-8 -*-
# @Time : 2021/1/31 21:09
# @Author : Jclian91
# @File : torch_model_predict.py
# @Place : Yangpu, Shanghai
import torch
import numpy as np
from transformers import BertTokenizer
from torch.utils.data import DataLoader
import os
import re
from load_data import bert_encode
from util import event_type, BERT_MODEL_DIR, MAX_LEN, VALID_BATCH_SIZE
from torch_model_train import BERTClass, CustomDataset, idx2tag

dev = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
model = BERTClass()
model.to(dev)
model.load_state_dict(torch.load("{}_ner.pth".format(event_type)))
tokenizer = BertTokenizer.from_pretrained('./{}'.format(BERT_MODEL_DIR))
path = "content"
outpath="predict_content"
mat="{:<15}"

def get_text_predict(text):
    test_set = CustomDataset(tokenizer, [text], [[]], MAX_LEN)
    test_params = {'batch_size': VALID_BATCH_SIZE, 'shuffle': True, 'num_workers': 0}
    test_loader = DataLoader(test_set, **test_params)

    model.eval()
    with torch.no_grad():
        for _, data in enumerate(test_loader):
            ids = data['ids'].to(dev, dtype=torch.long)
            mask = data['mask'].to(dev, dtype=torch.long)
            targets = data['tags'].to(dev, dtype=torch.long)

            output = model(ids, mask, labels=targets)
            loss, logits = output[:2]
            logits = logits.detach().cpu().numpy()
            tags = [idx2tag[_] for _ in [list(p) for p in np.argmax(logits, axis=2)][0]]

            # 输出预测结果
            real_tag = []
            i = 0
            for word in text.split():
                new_word = bert_encode(word)
                if i < len(tags):
                    real_tag.append(tags[i])
                    i += len(new_word)

    return real_tag
if __name__ == '__main__':
    files = os.listdir(path)
    for file in files:
        print(file)
        f1 = open(os.path.join(outpath, file), "w", encoding='utf-8')
        with open(os.path.join(path, file)) as f:
            test_text = f.read()
            for text_i in test_text.split(". "):
                text = re.sub(r'[^A-Za-z0-9 ]+', '', text_i)
                print(text)
                test_set = CustomDataset(tokenizer, [text], [[]], MAX_LEN)
                test_params = {'batch_size': VALID_BATCH_SIZE, 'shuffle': True, 'num_workers': 0}
                test_loader = DataLoader(test_set, **test_params)
                model.eval()
                with torch.no_grad():
                    for _, data in enumerate(test_loader):
                        ids = data['ids'].to(dev, dtype=torch.long)
                        mask = data['mask'].to(dev, dtype=torch.long)
                        targets = data['tags'].to(dev, dtype=torch.long)

                        output = model(ids, mask, labels=targets)
                        loss, logits = output[:2]
                        logits = logits.detach().cpu().numpy()
                        tags = [idx2tag[_] for _ in [list(p) for p in np.argmax(logits, axis=2)][0]]

                        # 输出预测结果
                        real_tag = []
                        i = 0
                        for word in text.split():
                            print(word)
                            new_word = bert_encode(word)
                            if i < len(tags):
                                # print(mat.format(word)+tags[i])
                                f1.write(word+ "\t" + tags[i] + "\n")
                                real_tag.append(tags[i])
                                i += len(new_word)
                f1.write("\n")
        f1.close()
