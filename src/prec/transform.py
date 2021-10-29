#!/usr/bin/python3

import csv


# transform list
# {3,0,2,0} corresponds to {1,1,1,3,3}
# add num of components to first coordinate
# {5,1,1,1,3,3}

def tr(li):
    out = []
    out.append(sum(li))
    for i in range(len(li)):
        for j in range(li[i]):
            out.append(i+1)
    return out

def stringify(vec):
    return "{" + ",".join([str(x) for x in vec]) + "}"

arr = []

stringify([1,2,3])

# read csv
with open('cacti_l.txt', newline='') as f:
    reader = csv.reader(f)
    for row in reader:
        line = "".join(row)             # make string out of list
        line = line[1:]                 # remove first {
        line = line.replace('}','')     # remove all }
        li = line.split('{')            # split substrings 
        # create array of int arrays and add it to list
        arr.append([[int(y) for y in x.split(' ')] for x in li])

outarr = [[tr(y)  for y in x] for x in arr]

linenumber=0
for line in outarr:
    linenumber = linenumber + 1
    vecindex = -1
    for vec in line:
        vecindex = vecindex + 1
        print("cactipart[{}][{}] = (INT *) malloc({}*sizeof(INT));".format(linenumber, vecindex, len(vec)))
        print("memcpy(cactipart[{}][{}], (INT []){}, {}*sizeof(INT));".format(linenumber, vecindex, stringify(vec), len(vec)))
    print("")


