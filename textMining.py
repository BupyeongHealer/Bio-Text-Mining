import copy
import math
from collections import deque

kegg_names = {}
name_kegg = {}
f1 = open('C:\\Users\Administrator\Desktop\genes.txt', 'r')
for line in f1.readlines():
    
    t1 = line.split(';')[0]
    t2 = t1.split('\t')
    kegg_id = t2[0]
    kegg_names[kegg_id] = []
    for name in t2[1].split(','):
        name = name.strip()
        kegg_names[kegg_id].append(name)
        name_kegg[name] = kegg_id
f1.close()
#print(name_kegg['CLCN5'])

f2 = open('C:\\Users\Administrator\Desktop\\abstracts1.txt', 'r')
data = f2.readlines()
f2.close()

print("Data",len(data))     #207248
temp = []
abstract = []
for i in range(0, 207248):
  temp.append(data[i])
  if "PMID" in data[i]:
    x = "".join(temp)
    abstract.append(x)
    #print("abstract", abstract)
    temp = []


file_data = "abstract-dentdiseas-set.txt"
handle_file = open(file_data, "r")
data = handle_file.readlines()
handle_file.close()
#print(len(data))  총 207248줄
#print(data[27])

documents = deque()
temp = []
abstract = []

for i in range(0, 207248):
  temp.append(data[i])
  if "MEDLINE" in data[i]:
    x = "".join(temp)
    abstract.append(x)
    #print("abstract", abstract)
    temp = []

#print(abstract[5])
#print(len(abstract)) #9234


def probability(data, target):
  length = len(data)
  count = 0
  for i in range(0, length):
    if target in data[i]:
      count += 1
  result = count/length
  return result

def probability2(data, target):
      length = len(data)
  count = 0
  for i in range(0, length):
    if target[0] in data[i]:
      count += 1
    elif target[1] in data[i]:
      count += 1
    elif target [2] in data[i]:
      count += 1
  result = count/length
  return result

def probability3(data, target1, target2):
      length = len(data)
  count = 0 
  for i in range(0, length):
    if target1 and target2[0] in data[i]:
      count += 1
    elif target1 and target2[1] in data[i]:
      count += 1
    elif target1 and target2[2] in data[i]:
      count += 1
  result = count/length
  return result

def main():
      gene_list = ["CLCN5", "OCRL", "OCRL1"]
  disease_probability = probability(abstract, "dent")
  gene_probability = probability2(abstract, gene_list)
  disease_gene_probability = probability3(abstract, "dent", gene_list)
  p_value = math.log2(disease_gene_probability / (disease_probability * gene_probability))
  
  print("disease_probability :", disease_probability)
  print("gene_probability :", gene_probability)
  print("==================================================")  
  print(gene_list[0], ":", gene_probability)
  print(gene_list[1], ":", gene_probability)
  print(gene_list[2], ":", gene_probability)
  print("==================================================")
  print("disease_gene_probability :", disease_gene_probability)
  print("p_value :", p_value)

if __name__ == "__main__":
  main()

"""
