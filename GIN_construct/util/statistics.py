with open("./data/label.test", 'r') as f:
    token_list = f.readlines()

sample_num = 1
labeled_sample_num = 0
flag = False

for t in token_list:
    if t != "\n":
        if t.strip("\n").split("\t")[1] != "O":
            if not flag:
                labeled_sample_num += 1
            flag = True
        continue
    else:
        sample_num += 1
        flag = False
print("sample number: {}".format(sample_num))
print("labeled sample number: {}".format(labeled_sample_num))