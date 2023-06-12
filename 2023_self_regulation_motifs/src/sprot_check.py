#!/usr/bin/env python3

import csv

def extact_uniprot_from_AF(af_id):
    s = af_id.split('-')
    if s[0] == 'AF':
        return s[1]
    return af_id

def overlap(a, b):
    if a is None or b is None:
        return False
    return a[0] <= b[0] <= a[1] or b[0] <= a[0] <= b[1]

table1 = open('SwissProt/SwissProt_results.csv', 'r')
reader1 = csv.reader(table1)
next(reader1)

table2 = open('autoinhibitory_sprot.csv', 'r')
reader2 = csv.reader(table2)
next(reader2)

table3 = open('sprot_check.csv', 'w')
writer = csv.writer(table3)
writer.writerow(["Accession", "Organism", "Type", "Location", "Sequence", "Note", "Found"])

sprot_dic = {}
for line in reader2:
    if line[0] not in sprot_dic.keys():
        sprot_dic[line[0]] = [line[1:]]
    else:
        sprot_dic[line[0]].append(line[1:])

found = {}
not_found = {}

c = 0
k = 0

for line in reader1:
    id = extact_uniprot_from_AF(line[1])
    if id in sprot_dic.keys():
        for data in sprot_dic[id]:
            # if the positions overlap we will condsider that we found the motif/domain
            if overlap(tuple((line[7], line[8])), data[2].strip('[]').split(':')):
                writer.writerow([id, data[0], data[1], data[2], data[3], data[4], True])
                print('\n', id, ':')
                print(data)
                print('-------------')
                print(line)
                if id not in found.keys():
                    found[id] = [data]
                else:
                    found[id].append(data)
                c+=1
            else:
                writer.writerow([id, data[0], data[1], data[2], data[3], data[4], False])
                if id not in not_found.keys():
                    not_found[id] = [data]
                else:
                    not_found[id].append(data)
                k+=1
print(c, 'found out of', c+k, '(', k, 'not found)')

table1.close()
table2.close()
table3.close()

