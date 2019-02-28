import json

i=0
with open("enumerate_disordered_whole_seq_2716.dat") as f:
    lst = json.load(f)

    print len(lst)

    for m in lst['models']:
        if m['score'] > 10000:
            i+=1


print i